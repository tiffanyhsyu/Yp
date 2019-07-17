import os
import emcee
import time
import numpy as np
import model_flux_ratio as mfr
import galaxy
import pdb


# Allowed galaxy names:
class MCMCgal:
    def __init__(self, galaxyname):
        self.galaxyname = galaxyname

        #galdict = galaxy.load_AOS2015(self.galaxyname) # optical+NIR
        #galdict = galaxy.load_AOS2012(self.galaxyname) # optical only
        #galdict = galaxy.load_synthetic(self.galaxyname) # synthetic runs
        #galdict = galaxy.load_ours_noHaHb(self.galaxyname)
        #galdict = galaxy.load_ours(self.galaxyname)
        galdict = galaxy.load_SDSS(self.galaxyname)

        self._full_tbl = galdict['full_tbl']
        self._T_OIII = galdict['T_OIII']

        # Read in measured data (wavelength, flux ratios, and EWs)
        #self._flux_ratios = self._full_tbl[:-1]  # Ignore the entry for P-gamma for MCMC'
        self._flux_ratios = self._full_tbl # Use full table for only optical systems..!
        #### Might want to remove this [:-1] in the future and instead, in the loop over emission lines, add an elif self._emis_lines[w] == 10941.082: continue, or something like that!

        #print (self._flux_ratios)
        # Names of wavelenghts of interest for MCMC
        # self._y_names = ['HeI+H83890', 'HeI4027', 'Hd', 'Hg', 'HeI4472', 'Hb', 'HeI5017', 'HeI5877', 'Ha', 'HeI6679', 'HeI7067', 'HeI10830']

        # Balmer and Helium lines of interest for MCMC
        self._hydrogen_lines = np.array([10941.082, 6564.612, 4862.721, 4341.684, 4102.891, 3890.166])  # Pa-g, Ha, Hb, Hg, Hd, H8
        self._helium_lines = np.array([10833.306, 7067.198, 6679.994, 5877.299, 5017.079, 4472.755, 4027.328, 3890.151])

        self._allowed_lines = np.sort(np.concatenate((self._hydrogen_lines, self._helium_lines)))[1:]  # [1:] to remove the duplicate ~3890 wavelength

        # Wavelengths we care about for MCMC, based on what is given in the input flux file (concatenating self._hydrogen_lines and self._helium_lines means some emlines could be mistakenly modeled even though they are not measured
        # Not sorting anymore because input could potentially not be in increasing Wavelength, and want to make sure we grab the right EW, Flux Ratios for the corresponding Wavelength
        self._emis_lines = self._flux_ratios['Wavelength']

        # Measured data from spectra
        self._EWs_meas = np.array(self._flux_ratios['EW'])

        self._y = np.array(self._flux_ratios['Flux Ratio'])  # F(lambda) / F(H-beta)
        try:
            self._y_error = np.array(self._flux_ratios['Flux Ratio Errors'])
        except:
            self._y_error = np.array(self._flux_ratios['Flux Ratio'] * 0.002)
        self._x = np.zeros(self._y.size)

        # Range of values for 8 parameters: y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, xi/n_HI
        self._min_y_plus, self._max_y_plus = 0.05, 0.1  # fraction of singly ionized He; y+ = He+/H+
        self._min_temp, self._max_temp = 10000, 22000  # electron temperature (K)
        self._min_log_dens, self._max_log_dens = 0, 3  # log10(electron density) (cm^-3)
        self._min_c_Hb, self._max_c_Hb = 0, 0.5  # reddening
        self._min_a_H, self._max_a_H = 0, 10  # underlying stellar H absorption (Angstroms)
        self._min_a_He, self._max_a_He = 0, 4  # underlying stellar HeI absorption (Angstroms)
        self._min_tau_He, self._max_tau_He = 0, 5  # optical depth; range of values from Izotov & Thuan 2010
        # min_n_HI, max_n_HI = 1e-4, 1e-1  # neutral hydrogen density (cm^-3)
        self._min_log_xi, self._max_log_xi = -6, -0.0969  # equals to xi=0-0.8; ratio of neutral hydrogen to singly ionized hydrogen densities; xi = n(HI)/n(HII)
        self.mcmc_steps()

    def __call__(self, theta, x, y, yerr):
        return self.lnprob(theta, x, y, yerr)

    # Set up MCMC
    def get_model(self, theta):
        #    y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, n_HI = theta
        y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi = theta

        # Reparameterize density and neutral hydrogen input
        dens = 10 ** log_dens
        xi = 10 ** log_xi

        # Take into account error on EW(Hb) by perturbing EW(Hb) by its measured EW error
        try:
            EW_Hb = np.random.normal(self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0],
                                 self._flux_ratios['EW Errors'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0])
        except:
            EW_Hb = np.random.normal(self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0],
                                 0.1*self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0])

        # Continuum level ratio; Eq. 2.4 of AOS2012
        h = self._y * EW_Hb / self._EWs_meas  # relative to H-beta; i.e., h(lambda) / h(Hbeta)
#        EW_meas = np.random.normal(self._EWs_meas, self._flux_ratios['EW Errors'])
#        EW_Hb = EW_meas[np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]]
#        h = self._y * EW_Hb / EW_meas  # relative to H-beta; i.e., h(lambda) / h(Hbeta)

        # Model flux
        model_flux = np.zeros(self._y.size)  # emission lines we want to model

        # Some values, calculated at Hbeta, for later use in model flux
        collisional_to_recomb_Hbeta = mfr.hydrogen_collision_to_recomb(xi, self._hydrogen_lines[2], temp, method='A2002')
        f_lambda_at_Hbeta = mfr.f_lambda_avg_interp(self._hydrogen_lines[2])
        #AHbeta_Av = mfr.reddening_coefficient(hydrogen_lines[2]) # CCM 1989 reddening curve

        done_3889 = False
        for w in range(len(self._emis_lines)):
            # Determine if working with hydrogen or helium line; within 3.5 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
            nearest_wave = self._allowed_lines[np.where(np.abs(self._allowed_lines - self._emis_lines[w]) < 3.5)[0]][0]
            # The above line is redundant for my input waves, but allows for cases where emis_lines[w] is some other array, say waves_of_interest[w],
            # and not exactly at the wavelengths given in the allowed_lines array (which is concatenated from arrays hydrogen_lines and helium_lines)

            # Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166) and P-gamma
            if nearest_wave in self._hydrogen_lines and nearest_wave != 3890.166 and nearest_wave != 10941.082:# and nearest_wave != 4862.721:
                line_species = 'hydrogen'

                emissivity_ratio = mfr.hydrogen_emissivity_S2018(self._emis_lines[w], temp, dens)
                a_H_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_H, ion=line_species)
                collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, self._emis_lines[w], temp, method='A2002')
                reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.
                #reddening_function = ( mfr.reddening_coefficient(self._emis_lines[w]) / AHbeta_Av ) - 1. # CCM 1989 reddening curve

                #			flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) *
                #				( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                #				10**-(reddening_function * c_Hb)
                # Reparameterization of flux to use the continuum level
                flux = (emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) *
                        10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_H_at_wave / EW_Hb) * (h[w]))

            # Any HeI line besides the blended HeI+H8 line or the NIR HeI10830 line
            elif nearest_wave in self._helium_lines and nearest_wave != 3890.151 and nearest_wave != 10833.306:
                line_species = 'helium'

                emissivity_ratio = mfr.helium_emissivity_PFSD2012(self._emis_lines[w], temp, dens)
                a_He_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_He, ion=line_species)
                optical_depth_at_wave = mfr.optical_depth_function(self._emis_lines[w], temp, dens, tau_He)
                reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.
                #reddening_function = ( mfr.reddening_coefficient(self._emis_lines[w]) / AHbeta_Av ) - 1. # CCM 1989 reddening curve

                flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                        10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_He_at_wave / EW_Hb) * (h[w]))

            # The blended HeI+H8 line
            elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
                if done_3889: continue
                done_3889 = True
                reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.
                #reddening_function = ( mfr.reddening_coefficient(self._emis_lines[w]) / AHbeta_Av ) - 1. # CCM 1989 reddening curve

                # HeI 3890.151 contribution:
                line_species = 'helium'

                emissivity_ratio = mfr.helium_emissivity_PFSD2012(self._emis_lines[w], temp, dens)
                a_He_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_He, ion=line_species)
                optical_depth_at_wave = mfr.optical_depth_function(self._emis_lines[w], temp, dens, tau_He)

                flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                        10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / EW_Hb)) - ((a_He_at_wave / EW_Hb) * (h[w]))

                # H8 contribution:
                line_species = 'hydrogen'

                emissivity_ratio = mfr.hydrogen_emissivity_S2018(self._emis_lines[w], temp, dens)
                a_H_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_H, ion=line_species)
                # C/R scaling for H8 from Hg already done in the function
                collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, self._emis_lines[w], temp, method='A2002')
#                collisional_to_recomb_factor = np.exp(( -13.6 * ((1/5**2)-(1/8**2)) ) / (8.6173303e-5 * temp)) # scale factor for going from C/R(Hg) to C/R(H8)
#                collisional_to_recomb_ratio = collisional_to_recomb_factor * mfr.hydrogen_collision_to_recomb(xi, 4341.684, temp) # Calculate C/R(Hg) and multiply by above scale factor

                flux += (emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) *
                         10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_H_at_wave / EW_Hb) * (h[w]))

            # Infrared HeI10830 line
            elif nearest_wave == 10833.306:
                # Theoretical F(Pg)/F(Hb) ratio, aka the 'model-dependent scaling ratio'
                line_species = 'hydrogen'

                emissivity_ratio = mfr.hydrogen_emissivity_S2018(10941.082, temp, dens)  # hard-coded Pg wavelength; could also be hydrogen_lines[0]
                a_H_at_wave = mfr.stellar_absorption(10941.082, a_H, ion=line_species)
                # C/R scaling for Pg from Pb already done in the function
                collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, self._hydrogen_lines[0], temp, method='A2002')
                reddening_function = (mfr.f_lambda_avg_interp(self._hydrogen_lines[0]) / f_lambda_at_Hbeta) - 1.  # hard-coded Pg wavelength; could also be hydrogen_lines[0]
                #reddening_function = ( mfr.reddening_coefficient(self._emis_lines[w]) / AHbeta_Av ) - 1. # CCM 1989 reddening curve

                EW_Pg = self._full_tbl[np.where(self._full_tbl['Wavelength'] == 10941.082)[0][0]]['EW']
                Pg_to_Hb_flux = emissivity_ratio *  ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) * \
                                ((EW_Hb + a_H) / (EW_Hb)) / ((EW_Pg + a_H_at_wave) / (EW_Pg)) * 10 ** -(reddening_function * c_Hb)

                # Theoretical F(HeI10830)/F(Hbeta) ratio
                line_species = 'helium'

                emissivity_ratio = mfr.helium_emissivity_PFSD2012(self._emis_lines[w], temp, dens)
                a_He_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_He, ion=line_species)
                optical_depth_at_wave = mfr.optical_depth_function(self._emis_lines[w], temp, dens, tau_He)
                reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.
                #reddening_function = ( mfr.reddening_coefficient(self._emis_lines[w]) / AHbeta_Av ) - 1. # CCM 1989 reddening curve

                # The way h is defined above and given the format of the input fluxes gives ( F(HeI10830)/F(Pg) ) * ( EW(Hb) / EW(HeI10830) ) here; must be multiplied by the calculated
                # theoretical F(Pg)/F(Hb) ratio from above to get the HeI10830 to Hbeta continuum level ratio, which is the definition of h, from Eq. 2.4 of AOS2012
                HeI10830_to_Hb_flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) * \
                                       10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ( (a_He_at_wave / EW_Hb) * (h[w] * Pg_to_Hb_flux) )

                # Want to get theoretical F(HeI10830)/F(Pg) to match that in the input table of flux_ratios -- can do this by ( F(HeI10830)/F(Hbeta) ) / ( F(Hbeta)/F(Pg) )!
                flux = HeI10830_to_Hb_flux / Pg_to_Hb_flux

            else:
                print ('Check your input wavelength -- not a recognized hydrogen or helium line for MCMC analysis')
                pdb.set_trace()

            model_flux[w] = flux

        return model_flux

    # Define the probability function as likelihood * prior.
    def lnprior(self, theta):
        #    y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, n_HI = theta
        y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi = theta

        if self._min_y_plus <= y_plus <= self._max_y_plus and \
                self._min_temp <= temp <= self._max_temp and \
                self._min_log_dens <= log_dens <= self._max_log_dens and \
                self._min_c_Hb <= c_Hb <= self._max_c_Hb and \
                self._min_a_H <= a_H <= self._max_a_H and \
                self._min_a_He <= a_He <= self._max_a_He and \
                self._min_tau_He <= tau_He <= self._max_tau_He and \
                self._min_log_xi <= log_xi <= self._max_log_xi:
            return -(temp - self._T_OIII) ** 2 / (0.2 * self._T_OIII) ** 2

        return -np.inf

    def lnlike(self, theta, x, y, yerr):
        model = self.get_model(theta)
        inv_sigma2 = 1.0 / yerr ** 2
        return -0.5 * (np.sum((y - model) ** 2 * inv_sigma2 - np.log(inv_sigma2)))

    def lnprob(self, theta, x, y, yerr):
        lp = self.lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.lnlike(theta, x, y, yerr)

    def mcmc_steps(self):
        # Set up sampler
        ndim, nwalkers = 8, 500

        pos = [np.array([np.random.uniform(self._min_y_plus, self._max_y_plus),
                         np.random.uniform(self._min_temp, self._max_temp),
                         np.random.uniform(self._min_log_dens, self._max_log_dens),
                         np.random.uniform(self._min_c_Hb, self._max_c_Hb),
                         np.random.uniform(self._min_a_H, self._max_a_H),
                         np.random.uniform(self._min_a_He, self._max_a_He),
                         np.random.uniform(self._min_tau_He, self._max_tau_He),
                         # np.random.uniform(min_n_HI, max_n_HI),
                         np.random.uniform(self._min_log_xi, self._max_log_xi)]) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self, args=(self._x, self._y, self._y_error), threads=ndim)

        print('Running MCMC...')
        nmbr = 1000
        a = time.time()
        for i, result in enumerate(sampler.run_mcmc(pos, nmbr, rstate0=np.random.get_state())):
            if True:  # (i+1) % 100 == 0:
                print('{0:5.1%}'.format(float(i) / nmbr))
        print('Done!')
        print((time.time() - a) / 60.0, 'mins')

        print('Saving samples')
        np.save('{0:s}_{1:d}walkers_{2:d}steps'.format(self.galaxyname, nwalkers, nmbr), sampler.chain)

        print('Making plots')
        burnin = int(0.8 * nmbr)

        samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
        # Names of 8 parameters and input 'true' parameter values
        prenams = ['y+', 'temperature', '$log(n_{e})$', 'c(H\\beta)', '$a_{H}$', '$a_{He}$', '$\\tau_{He}', '$log(\\xi)$']
        input_vals = np.array([0.08, 18000, 2, 0.1, 1.0, 1.0, 1.0, -2])  # Input parameters for fake spectra
        # input_vals = np.array([0.08634, 12979, 1.987, 0.15, 2.31, 0.37, 2.27, -1.767]) # AOS 2015's solved parameters for Mrk450 No.1

        print ('Best parameter values:')
        # y_plus_mcmc, temp_mcmc, dens_mcmc, c_Hb_mcmc, a_H_mcmc, a_He_mcmc, tau_He_mcmc, n_HI_mcmc = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        y_plus_mcmc, temp_mcmc, log_dens_mcmc, c_Hb_mcmc, a_H_mcmc, a_He_mcmc, tau_He_mcmc, log_xi_mcmc = map(
            lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))

        print ('y+', y_plus_mcmc)
        print ('T', temp_mcmc)
        print ('log(n_e)', log_dens_mcmc)
        print ('c(Hb)', c_Hb_mcmc)
        print ('a_H', a_H_mcmc)
        print ('a_He', a_He_mcmc)
        print ('tau_He', tau_He_mcmc)
        # print ('n_HI', n_HI_mcmc)
        print ('log(xi)', log_xi_mcmc)
        #print ('\n Input parameter values:')
        #print (input_vals)

        dens = 10.0 ** (samples[:, 2])
        v = np.percentile(dens, [16, 50, 84])
        dens_mcmc = (v[1], v[2] - v[1], v[1] - v[0],)
        # Save some output
        allpars = np.hstack((y_plus_mcmc, dens_mcmc, a_He_mcmc, tau_He_mcmc, temp_mcmc, c_Hb_mcmc, a_H_mcmc, log_xi_mcmc))
        outdat = open('all_output', 'r').readlines()
        sendout = open('all_output', 'w')
        for ii in outdat:
            sendout.write(ii)
        sendout.write('{0:s} '.format(self.galaxyname))
        for ii in allpars:
            sendout.write('{0:f} '.format(ii))
        sendout.write('\n')
        sendout.close()

        '''
        fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, 12))
        for i in range(ndim):
            axes[i].plot(sampler.chain[:, :, i].T, color='k', alpha=0.4)
            axes[i].yaxis.set_major_locator(MaxNLocator(5))
            axes[i].axvline(burnin, color='red')
            axes[i].set_ylabel(prenams[i])
        axes[7].set_xlabel('Steps')
        fig.tight_layout(h_pad=0.0)
        fig.savefig('{0:s}_{1:d}_walkers{2:d}_steps.pdf'.format('test_MCMC_time_evol', nwalkers, nmbr), overwrite=True)
        #fig.savefig('{0:s}_{1:d}_walkers{2:d}_steps.pdf'.format('LeoP_time_evol', nwalkers, nmbr), overwrite=True)
        
        fig = corner.corner(samples, labels=prenams, truths=input_vals)
        fig.savefig('{0:s}_{1:d}_walkers{2:d}_steps.pdf'.format('test_MCMC_params', nwalkers, nmbr), overwrite=True)
        #fig.savefig('{0:s}_{1:d}_walkers{2:d}_steps.pdf'.format('LeoP_MCMC_params', nwalkers, nmbr), overwrite=True)
        '''


if __name__ == '__main__':
    # The allowed names
    AOS2015 = ['IZw18SE1', 'SBS0335-052E1', 'SBS0335-052E3', 'J0519+0007', 'SBS0940+5442', 'Tol65', 'SBS1415+437No13',
             'SBS1415+437No2', 'CGCG007-025No2', 'Mrk209', 'SBS1030+583', 'Mrk71No1', 'SBS1152+579', 'Mrk59',
             'SBS1135+581', 'Mrk450No1']
    ours = ['J0000p3052A', 'J0000p3052B', 'J0018p2345', 'J0118p3512', 'J0140p2951', 'J0214m0835', 'J0220p2044A',
            'J0220p2044B', 'J0452m0541', 'J0743p4807', 'J0757p4750', 'J0943p3326', 'J1044p6306', 'J1204p5259',
            'J1214p1245', 'J1322p5425', 'J1414m0208', 'J1425p4441', 'J1655p6337', 'J1705p3527', 'J1732p4452', 'J1757p6454',
            'J2030m1343', 'J2213p1722', 'J2319p1616', 'J2230m0531', 'J2339p3230', 'KJ2', 'KJ29', 'KJ5', 'KJ5B', 'KJ97', 'LeoP']
    SDSS = ['spec-0266-51630-0407', 'spec-0284-51943-0408', 'spec-0283-51959-0572', 'spec-0284-51943-0007', 'spec-0278-51900-0392',
            'spec-0299-51671-0083', 'spec-0299-51671-0311', 'spec-0289-51990-0369', 'spec-0285-51930-0154', 'spec-0267-51608-0421',
            'spec-0279-51984-0293', 'spec-0279-51984-0520', 'spec-0301-51942-0531', 'spec-0327-52294-0042', 'spec-0287-52023-0230',
            'spec-0278-51900-0081', 'spec-0304-51609-0020', 'spec-0304-51609-0034', 'spec-0380-51792-0253', 'spec-0394-51913-0493',
            'spec-0394-51913-0599', 'spec-0280-51612-0192', 'spec-0281-51614-0129', 'spec-0358-51818-0311', 'spec-0358-51818-0364',
            'spec-0358-51818-0403', 'spec-0358-51818-0504', 'spec-0328-52282-0491', 'spec-0292-51609-0566', 'spec-0287-52023-0236',
            'spec-0304-51609-0583', 'spec-0339-51692-0437', 'spec-0281-51614-0499', 'spec-0282-51658-0296', 'spec-0296-51984-0416',
            'spec-0306-51637-0515', 'spec-0381-51811-0370', 'spec-0305-51613-0594', 'spec-0339-51692-0622', 'spec-0282-51658-0543',
            'spec-0308-51662-0081', 'spec-0308-51662-0130', 'spec-0446-51899-0086', 'spec-0446-51899-0283', 'spec-0446-51899-0352',
            'spec-0351-51695-0217', 'spec-0309-51994-0288', 'spec-0420-51871-0156', 'spec-0420-51871-0474', 'spec-0423-51821-0402',
            'spec-0383-51818-0266', 'spec-0384-51821-0281', 'spec-0467-51901-0628', 'spec-0364-52000-0384', 'spec-0309-51994-0489',
            'spec-0310-51990-0136', 'spec-0443-51873-0542', 'spec-0456-51910-0195', 'spec-0456-51910-0306', 'spec-0456-51910-0546',
            'spec-0395-51783-0570', 'spec-0364-52000-0066', 'spec-0364-52000-0187', 'spec-0275-51910-0015', 'spec-0275-51910-0445',
            'spec-0501-52235-0361', 'spec-0367-51997-0561', 'spec-0312-51689-0508', 'spec-0331-52368-0006', 'spec-0345-51690-0572',
            'spec-0375-52140-0348', 'spec-0465-51910-0524', 'spec-0353-51703-0486', 'spec-0335-52000-0604', 'spec-0336-51999-0033',
            'spec-0505-52317-0228', 'spec-0507-52353-0087', 'spec-0412-52258-0441', 'spec-0412-52258-0518', 'spec-0472-51955-0235',
            'spec-0414-51869-0459', 'spec-0414-51869-0524', 'spec-0375-52140-0118', 'spec-0390-51900-0430', 'spec-0390-51900-0445',
            'spec-0530-52026-0184', 'spec-0394-51913-0402', 'spec-0394-51913-0469', 'spec-0394-51913-0472', 'spec-0341-51690-0256',
            'spec-0342-51691-0081', 'spec-0481-51908-0483', 'spec-0415-51810-0285', 'spec-0530-52026-0525', 'spec-0348-51671-0331',
            'spec-0487-51943-0583', 'spec-0439-51877-0086', 'spec-0450-51908-0377', 'spec-0546-52205-0419', 'spec-0550-51959-0092',
            'spec-0454-51908-0570', 'spec-0455-51909-0047', 'spec-0349-51699-0337', 'spec-0492-51955-0449', 'spec-0493-51957-0219', 'spec-0497-51989-0198', 'spec-0450-51908-0520', 'spec-0450-51908-0545', 'spec-0436-51883-0600', 'spec-0564-52224-0216', 'spec-0564-52224-0223', 'spec-0569-52264-0609', 'spec-0570-52266-0061', 'spec-0351-51695-0049', 'spec-0452-51911-0487', 'spec-0502-51957-0007', 'spec-0494-51915-0537', 'spec-0409-51871-0033', 'spec-0409-51871-0050', 'spec-0437-51876-0460', 'spec-0575-52319-0521', 'spec-0488-51914-0439', 'spec-0524-52027-0016', 'spec-0507-52353-0377', 'spec-0441-51868-0634', 'spec-0442-51882-0156', 'spec-0442-51882-0413', 'spec-0516-52017-0315', 'spec-0516-52017-0403', 'spec-0516-52017-0459', 'spec-0518-52282-0335', 'spec-0425-51898-0634', 'spec-0447-51877-0172', 'spec-0447-51877-0361', 'spec-0597-52059-0586', 'spec-0406-51817-0490', 'spec-0469-51913-0009', 'spec-0469-51913-0110', 'spec-0536-52024-0326', 'spec-0525-52295-0626', 'spec-0535-51999-0052', 'spec-0519-52283-0124', 'spec-0498-51984-0521', 'spec-0564-52224-0283', 'spec-0483-51924-0594', 'spec-0431-51877-0191', 'spec-0463-51908-0299', 'spec-0499-51988-0124', 'spec-0512-51992-0524', 'spec-0519-52283-0615', 'spec-0568-52254-0244', 'spec-0463-51908-0470', 'spec-0490-51929-0279', 'spec-0490-51929-0363', 'spec-0490-51929-0396', 'spec-0490-51929-0401', 'spec-0490-51929-0519', 'spec-0494-51915-0007', 'spec-0630-52050-0395', 'spec-0424-51893-0279', 'spec-0531-52028-0218', 'spec-0469-51913-0498', 'spec-0472-51955-0546', 'spec-0496-51988-0541', 'spec-0483-51924-0495', 'spec-0495-51988-0466', 'spec-0541-51959-0025', 'spec-0552-51992-0027', 'spec-0603-52056-0486', 'spec-0607-52368-0239', 'spec-0544-52201-0067', 'spec-0555-52266-0517', 'spec-0585-52027-0255', 'spec-0473-51929-0422', 'spec-0503-51999-0123', 'spec-0503-51999-0273', 'spec-0498-51984-0236', 'spec-0664-52174-0447', 'spec-0459-51924-0253', 'spec-0553-51999-0342', 'spec-0648-52559-0288', 'spec-0649-52201-0019', 'spec-0649-52201-0117', 'spec-0656-52148-0243', 'spec-0566-52238-0497', 'spec-0570-52266-0586', 'spec-0587-52026-0495', 'spec-0587-52026-0546', 'spec-0527-52342-0103', 'spec-0515-52051-0378', 'spec-0519-52283-0507', 'spec-0519-52283-0515', 'spec-0460-51924-0128', 'spec-0520-52288-0551', 'spec-0545-52202-0339', 'spec-0625-52145-0432', 'spec-0637-52174-0526', 'spec-0614-53437-0055', 'spec-0523-52026-0114', 'spec-0524-52027-0260', 'spec-0529-52025-0257', 'spec-0682-52525-0172', 'spec-0584-52049-0421', 'spec-0691-52199-0252', 'spec-0704-52205-0494', 'spec-0633-52079-0191', 'spec-0541-51959-0600', 'spec-0542-51993-0013', 'spec-0539-52017-0155', 'spec-0539-52017-0356', 'spec-0615-52347-0272', 'spec-0622-52054-0445', 'spec-0641-52176-0213', 'spec-0566-52238-0338', 'spec-0544-52201-0610', 'spec-0553-51999-0134', 'spec-0542-51993-0269', 'spec-0621-52055-0618', 'spec-0621-52055-0624', 'spec-0664-52174-0112', 'spec-0639-52146-0143', 'spec-0639-52146-0242', 'spec-0658-52146-0017', 'spec-0590-52057-0200', 'spec-0660-52177-0414', 'spec-0666-52149-0492', 'spec-0556-51991-0224', 'spec-0629-52051-0497', 'spec-0686-52519-0185', 'spec-0650-52143-0330', 'spec-0677-52606-0533', 'spec-0608-52081-0139', 'spec-0615-52347-0590', 'spec-0569-52264-0472', 'spec-0653-52145-0056', 'spec-0668-52162-0089', 'spec-0786-52319-0443', 'spec-0795-52378-0405', 'spec-0662-52147-0333', 'spec-0692-52201-0465', 'spec-0665-52168-0324', 'spec-0620-52375-0428', 'spec-0690-52261-0362', 'spec-0582-52045-0440', 'spec-0582-52045-0445', 'spec-0679-52177-0605', 'spec-0712-52199-0158', 'spec-0673-52162-0073', 'spec-0721-52228-0555', 'spec-0624-52377-0422', 'spec-0630-52050-0463', 'spec-0640-52200-0270', 'spec-0805-52586-0482', 'spec-0707-52177-0525', 'spec-0708-52175-0155', 'spec-0602-52072-0019', 'spec-0684-52523-0560', 'spec-0686-52519-0406', 'spec-0695-52202-0261', 'spec-0663-52145-0406', 'spec-0834-52316-0157', 'spec-0602-52072-0369', 'spec-0602-52072-0500', 'spec-0831-52294-0526', 'spec-0831-52294-0572', 'spec-0840-52374-0570', 'spec-0749-52226-0500', 'spec-0755-52235-0416', 'spec-0742-52263-0179', 'spec-0675-52590-0094', 'spec-0687-52518-0274', 'spec-0708-52175-0404', 'spec-0606-52365-0604', 'spec-0750-52235-0312', 'spec-0769-54530-0086', 'spec-0771-52370-0490', 'spec-0634-52164-0230', 'spec-0879-52365-0014', 'spec-0879-52365-0097', 'spec-0904-52381-0501', 'spec-0762-52232-0217', 'spec-0762-52232-0575', 'spec-0764-52238-0477', 'spec-0881-52368-0567', 'spec-0894-52615-0185', 'spec-0755-52235-0300', 'spec-0730-52466-0376', 'spec-0726-52226-0532', 'spec-0737-52518-0544', 'spec-0651-52141-0452', 'spec-0653-52145-0462', 'spec-0663-52145-0288', 'spec-0782-52320-0022', 'spec-0662-52147-0466', 'spec-0758-52253-0176', 'spec-0903-52400-0631', 'spec-0761-54524-0477', 'spec-0804-52286-0085', 'spec-0735-52519-0461', 'spec-0740-52263-0259', 'spec-0740-52263-0267', 'spec-0676-52178-0192', 'spec-0922-52426-0212', 'spec-0933-52642-0080', 'spec-0802-52289-0423', 'spec-0785-52339-0079', 'spec-0767-52252-0463', 'spec-0911-52426-0253', 'spec-0766-52247-0492', 'spec-0775-52295-0029', 'spec-0679-52177-0055', 'spec-0841-52375-0333', 'spec-0703-52209-0285', 'spec-0709-52205-0584', 'spec-0775-52295-0291', 'spec-0838-52378-0595', 'spec-0792-52353-0110', 'spec-0764-52238-0633', 'spec-0779-52342-0218', 'spec-0783-52325-0594', 'spec-0964-52646-0570', 'spec-0971-52644-0599', 'spec-0713-52178-0481', 'spec-0951-52398-0017', 'spec-0849-52439-0598', 'spec-0850-52338-0010', 'spec-0864-52320-0124', 'spec-0780-52370-0398', 'spec-0884-52374-0216', 'spec-0865-52323-0080', 'spec-0844-52378-0299', 'spec-0875-52354-0226', 'spec-0908-52373-0027', 'spec-0909-52379-0616', 'spec-1002-52646-0525', 'spec-0847-52426-0522', 'spec-0854-52373-0514', 'spec-0778-54525-0497', 'spec-0989-52468-0143', 'spec-0741-52261-0279', 'spec-0898-52606-0045', 'spec-0843-52378-0001', 'spec-0845-52381-0001', 'spec-0846-52407-0267', 'spec-1007-52706-0478', 'spec-0890-52583-0314', 'spec-0816-52379-0277', 'spec-0940-52670-0594', 'spec-0949-52427-0568', 'spec-0951-52398-0142', 'spec-0945-52652-0423', 'spec-0952-52409-0439', 'spec-0955-52409-0383', 'spec-0957-52398-0530', 'spec-0892-52378-0229', 'spec-0854-52373-0373', 'spec-0832-52312-0473', 'spec-1039-52707-0119', 'spec-0960-52425-0581', 'spec-0885-52379-0301', 'spec-0967-52636-0244', 'spec-0967-52636-0339', 'spec-0829-52296-0210', 'spec-0897-52605-0553', 'spec-1014-52707-0393', 'spec-1015-52734-0003', 'spec-0960-52425-0368', 'spec-0963-52643-0395', 'spec-0877-52353-0253', 'spec-0836-52376-0378', 'spec-0991-52707-0004', 'spec-0992-52644-0524', 'spec-1017-52706-0153', 'spec-0848-52669-0528', 'spec-1003-52641-0327', 'spec-1008-52707-0430', 'spec-0952-52409-0264', 'spec-0991-52707-0488', 'spec-0890-52583-0049', 'spec-0890-52583-0065', 'spec-0900-52637-0388', 'spec-0902-52409-0135', 'spec-0902-52409-0276', 'spec-0884-52374-0404', 'spec-0850-52338-0397', 'spec-0853-52374-0577', 'spec-0962-52620-0227', 'spec-1061-52641-0393', 'spec-0962-52620-0028', 'spec-1018-52672-0311', 'spec-0914-52721-0052', 'spec-0893-52589-0177', 'spec-1031-53172-0544', 'spec-0859-52317-0070', 'spec-0972-52435-0370', 'spec-0926-52413-0279', 'spec-0950-52378-0214', 'spec-0922-52426-0547', 'spec-0900-52637-0234', 'spec-0917-52400-0336', 'spec-0896-52592-0125', 'spec-0992-52644-0019', 'spec-1096-52974-0102', 'spec-1005-52703-0452', 'spec-0957-52398-0210', 'spec-0959-52411-0137', 'spec-0937-52707-0275', 'spec-0944-52614-0487', 'spec-1177-52824-0556', 'spec-1177-52824-0616', 'spec-1184-52641-0581', 'spec-1195-52724-0060', 'spec-1054-52516-0499', 'spec-0899-52620-0594', 'spec-0900-52637-0085', 'spec-0999-52636-0150', 'spec-1020-52721-0287', 'spec-0930-52618-0274', 'spec-1040-52722-0318', 'spec-1040-52722-0358', 'spec-1040-52722-0494', 'spec-1025-53239-0103', 'spec-1207-52672-0506', 'spec-1173-52790-0240', 'spec-1176-52791-0591', 'spec-1108-53227-0397', 'spec-1012-52649-0392', 'spec-0938-52708-0627', 'spec-1215-52725-0501', 'spec-1185-52642-0580', 'spec-1185-52642-0630', 'spec-1040-52722-0021', 'spec-1028-52884-0562', 'spec-1034-52813-0083', 'spec-1034-52813-0521', 'spec-1240-52734-0340', 'spec-0949-52427-0359', 'spec-1048-52736-0424', 'spec-1052-52466-0115', 'spec-0963-52643-0165', 'spec-0960-52425-0063', 'spec-1096-52974-0288', 'spec-1156-52641-0423', 'spec-1159-52669-0305', 'spec-1050-52721-0274', 'spec-1050-52721-0402', 'spec-0976-52413-0238', 'spec-0987-52523-0394', 'spec-1158-52668-0120', 'spec-1163-52669-0465', 'spec-1233-52734-0335', 'spec-1116-52932-0089', 'spec-1005-52703-0268', 'spec-1283-52762-0315', 'spec-0999-52636-0517', 'spec-1279-52736-0147', 'spec-1269-52937-0177', 'spec-1212-52703-0388', 'spec-1212-52703-0606', 'spec-1073-52649-0409', 'spec-1073-52649-0419', 'spec-1104-52912-0511', 'spec-1021-52460-0446', 'spec-1217-52672-0271', 'spec-1209-52674-0598', 'spec-1284-52736-0256', 'spec-1293-52765-0589', 'spec-1214-52731-0142', 'spec-1214-52731-0339', 'spec-1238-52761-0515', 'spec-1028-52884-0386', 'spec-1028-52884-0401', 'spec-1030-52914-0077', 'spec-1030-52914-0107', 'spec-1170-52756-0485', 'spec-1324-53088-0524', 'spec-1163-52669-0284', 'spec-1171-52753-0376', 'spec-1346-52822-0185', 'spec-1338-52765-0329', 'spec-1313-52790-0219', 'spec-1164-52674-0418', 'spec-1053-52468-0249', 'spec-1355-52823-0098', 'spec-1363-53053-0599', 'spec-1365-53062-0338', 'spec-1279-52736-0546', 'spec-1074-52937-0186', 'spec-1199-52703-0600', 'spec-1267-52932-0384', 'spec-1331-52766-0554', 'spec-1213-52972-0088', 'spec-1059-52618-0309', 'spec-1059-52618-0359', 'spec-1369-53089-0085', 'spec-1286-52725-0579', 'spec-1093-52591-0261', 'spec-1213-52972-0500', 'spec-1233-52734-0136', 'spec-1203-52669-0570', 'spec-1073-52649-0188', 'spec-1074-52937-0573', 'spec-1306-52996-0230', 'spec-1306-52996-0292', 'spec-1215-52725-0273', 'spec-1394-53108-0294', 'spec-1373-53063-0205', 'spec-1399-53172-0299', 'spec-1299-52972-0486', 'spec-1286-52725-0150', 'spec-1333-52782-0172', 'spec-1349-52797-0175', 'spec-1158-52668-0062', 'spec-1319-52791-0443', 'spec-1424-52912-0189', 'spec-1401-53144-0397', 'spec-1350-52786-0300', 'spec-1467-53115-0032', 'spec-1357-53034-0540', 'spec-1345-52814-0364', 'spec-1457-53116-0394', 'spec-1425-52913-0518', 'spec-1304-52993-0210', 'spec-1314-53050-0214', 'spec-1268-52933-0318', 'spec-1275-52996-0232', 'spec-1277-52765-0031', 'spec-1162-52668-0224', 'spec-1313-52790-0408', 'spec-1356-53033-0248', 'spec-1370-53090-0254', 'spec-1377-53050-0362', 'spec-1461-53062-0068', 'spec-1464-53091-0442', 'spec-1467-53115-0579', 'spec-1455-53089-0569', 'spec-1322-52791-0470', 'spec-1326-52764-0028', 'spec-1376-53089-0637', 'spec-1310-53033-0508', 'spec-1223-52781-0287', 'spec-1385-53108-0015', 'spec-1385-53108-0559', 'spec-1466-53083-0092', 'spec-1316-52790-0581', 'spec-1590-52974-0628', 'spec-1595-52999-0624', 'spec-1487-52964-0510', 'spec-1362-53050-0617', 'spec-1342-52793-0537', 'spec-1607-53083-0205', 'spec-1405-52826-0395', 'spec-1408-52822-0221', 'spec-1350-52786-0090', 'spec-1372-53062-0517', 'spec-1215-52725-0629', 'spec-1624-53386-0636', 'spec-1292-52736-0544', 'spec-1301-52976-0356', 'spec-1387-53118-0458', 'spec-1576-53496-0446', 'spec-1428-52998-0568', 'spec-1430-53002-0389', 'spec-1442-53050-0548', 'spec-1357-53034-0202', 'spec-1457-53116-0533', 'spec-1596-52998-0479', 'spec-1597-52999-0050', 'spec-1455-53089-0287', 'spec-1417-53141-0306', 'spec-1651-53442-0255', 'spec-1450-53120-0201', 'spec-1442-53050-0396', 'spec-1682-53173-0277', 'spec-1492-52932-0010', 'spec-1573-53226-0512', 'spec-1584-52943-0372', 'spec-1620-53137-0470', 'spec-1464-53091-0232', 'spec-1310-53033-0488', 'spec-1689-53177-0521', 'spec-1363-53053-0411', 'spec-1372-53062-0402', 'spec-1615-53166-0206', 'spec-1495-52944-0627', 'spec-1579-53473-0103', 'spec-1434-53053-0021', 'spec-1314-53050-0582', 'spec-1325-52762-0412', 'spec-1705-53848-0450', 'spec-1602-53117-0460', 'spec-1676-53147-0294', 'spec-1685-53463-0234', 'spec-1685-53463-0440', 'spec-1622-53385-0193', 'spec-1458-53119-0543', 'spec-1462-53112-0184', 'spec-1502-53741-0087', 'spec-1712-53531-0569', 'spec-1721-53857-0206', 'spec-1578-53496-0438', 'spec-1627-53473-0426', 'spec-1695-53473-0627', 'spec-1576-53496-0378', 'spec-1623-53089-0454', 'spec-1344-52792-0396', 'spec-1351-52790-0437', 'spec-1351-52790-0474', 'spec-1374-53083-0560', 'spec-1741-53052-0219', 'spec-1588-52965-0598', 'spec-1707-53885-0025', 'spec-1708-53503-0492', 'spec-1708-53503-0512', 'spec-1724-53859-0322', 'spec-1388-53119-0189', 'spec-1744-53055-0450', 'spec-1597-52999-0358', 'spec-1605-53062-0061', 'spec-1733-53047-0240', 'spec-1609-53142-0238', 'spec-1418-53142-0442', 'spec-1771-53498-0429', 'spec-1682-53173-0003', 'spec-1761-53376-0601', 'spec-1722-53852-0358', 'spec-1431-52992-0451', 'spec-1433-53035-0376', 'spec-1782-53299-0624', 'spec-1801-54156-0583', 'spec-1803-54152-0448', 'spec-1814-54555-0293', 'spec-1814-54555-0395', 'spec-1644-53144-0564', 'spec-1693-53446-0541', 'spec-1708-53503-0263', 'spec-1616-53169-0205', 'spec-1626-53472-0539', 'spec-1795-54507-0315', 'spec-1799-53556-0055', 'spec-1800-53884-0520', 'spec-1647-53531-0034', 'spec-1649-53149-0344', 'spec-1714-53521-0042', 'spec-1827-53531-0228', 'spec-1840-53472-0032', 'spec-1745-53061-0463', 'spec-1745-53061-0475', 'spec-1752-53379-0530', 'spec-1753-53383-0146', 'spec-1666-52991-0310', 'spec-1713-53827-0432', 'spec-1718-53850-0484', 'spec-1718-53850-0542', 'spec-1570-53149-0146', 'spec-1580-53145-0419', 'spec-1647-53531-0344', 'spec-1455-53089-0076', 'spec-1719-53876-0008', 'spec-1719-53876-0196', 'spec-1679-53149-0534', 'spec-1734-53034-0490', 'spec-1650-53174-0305', 'spec-1461-53062-0582', 'spec-1820-54208-0585', 'spec-1702-53144-0062', 'spec-1715-54212-0211', 'spec-1771-53498-0627', 'spec-1801-54156-0525', 'spec-1747-53075-0182', 'spec-1749-53357-0173', 'spec-1870-53383-0072', 'spec-1724-53859-0627', 'spec-1725-54266-0068', 'spec-1571-53174-0155', 'spec-1817-53851-0358', 'spec-1819-54540-0638', 'spec-1792-54270-0119', 'spec-1697-53142-0048', 'spec-1706-53442-0099', 'spec-1733-53047-0326', 'spec-1733-53047-0528', 'spec-1778-53883-0143', 'spec-1926-53317-0369', 'spec-1823-53886-0464', 'spec-1839-53471-0483', 'spec-1763-53463-0225', 'spec-1774-53759-0060', 'spec-1806-53559-0414', 'spec-1601-53115-0168', 'spec-1610-53144-0379', 'spec-1758-53084-0338', 'spec-1946-53432-0528', 'spec-1953-53358-0618', 'spec-1785-54439-0341', 'spec-1594-52992-0343', 'spec-1622-53385-0403', 'spec-1983-53442-0308', 'spec-1872-53386-0526', 'spec-1799-53556-0392', 'spec-1850-53786-0209', 'spec-1690-53475-0108', 'spec-1702-53144-0234', 'spec-1775-53847-0305', 'spec-1938-53379-0406', 'spec-1802-53885-0346', 'spec-1824-53491-0063', 'spec-1836-54567-0353', 'spec-1836-54567-0390', 'spec-1834-54562-0290', 'spec-1834-54562-0306', 'spec-2037-53446-0293', 'spec-1936-53330-0020', 'spec-1937-53388-0488', 'spec-1950-53436-0189', 'spec-1843-53816-0039', 'spec-1848-54180-0463', 'spec-1996-53436-0475', 'spec-2003-53442-0450', 'spec-1842-53501-0210', 'spec-1678-53433-0425', 'spec-1864-53313-0228', 'spec-1957-53415-0258', 'spec-1749-53357-0026', 'spec-1763-53463-0094', 'spec-2031-53848-0009', 'spec-2032-53815-0531', 'spec-1847-54176-0496', 'spec-1851-53524-0225', 'spec-1710-53504-0632', 'spec-1876-54464-0379', 'spec-1974-53430-0155', 'spec-1929-53349-0593', 'spec-1785-54439-0201', 'spec-1684-53239-0484', 'spec-1881-53261-0528', 'spec-1706-53442-0015', 'spec-1744-53055-0385', 'spec-1922-53315-0588', 'spec-1924-53330-0612', 'spec-1927-53321-0006', 'spec-2098-53460-0291', 'spec-1949-53433-0592', 'spec-1748-53112-0217', 'spec-1748-53112-0218', 'spec-2100-53713-0211', 'spec-1942-53415-0577', 'spec-1767-53436-0514', 'spec-1774-53759-0638', 'spec-1957-53415-0535', 'spec-2112-53534-0555', 'spec-1875-54453-0549', 'spec-1754-53385-0131', 'spec-1959-53440-0545', 'spec-2123-53793-0288', 'spec-2132-53493-0241', 'spec-1981-53463-0426', 'spec-1981-53463-0438', 'spec-1786-54450-0448', 'spec-2132-53493-0445', 'spec-2154-54539-0217', 'spec-1993-53762-0564', 'spec-1999-53503-0119', 'spec-2029-53819-0505', 'spec-1997-53442-0006', 'spec-1795-54507-0461', 'spec-2156-54525-0081', 'spec-2156-54525-0140', 'spec-2165-53917-0022', 'spec-2010-53495-0338', 'spec-2010-53495-0127', 'spec-1809-53792-0586', 'spec-1824-53491-0320', 'spec-2028-53818-0494', 'spec-1864-53313-0359', 'spec-2114-53848-0381', 'spec-2136-53494-0081', 'spec-1935-53387-0204', 'spec-2016-53799-0476', 'spec-2167-53889-0429', 'spec-2173-53874-0127', 'spec-2101-53858-0387', 'spec-1854-53566-0551', 'spec-2095-53474-0074', 'spec-2231-53816-0443', 'spec-2166-54232-0594', 'spec-2173-53874-0538', 'spec-2117-54115-0484', 'spec-2130-53881-0238', 'spec-1986-53475-0260', 'spec-1868-53318-0353', 'spec-2107-53786-0443', 'spec-1995-53415-0191', 'spec-2098-53460-0503', 'spec-2088-53493-0585', 'spec-2090-53463-0559', 'spec-2092-53460-0557', 'spec-2155-53820-0479', 'spec-2088-53493-0353', 'spec-2108-53473-0077', 'spec-1932-53350-0565', 'spec-2138-53757-0152', 'spec-2148-54526-0555', 'spec-2270-53714-0470', 'spec-2283-53729-0182', 'spec-2110-53467-0466', 'spec-2110-53467-0496', 'spec-2036-53446-0598', 'spec-2103-53467-0109', 'spec-2106-53714-0487', 'spec-2115-53535-0414', 'spec-2150-54510-0397', 'spec-2154-54539-0318', 'spec-2154-54539-0338', 'spec-2240-53823-0220', 'spec-2123-53793-0003', 'spec-2126-53794-0436', 'spec-2196-53534-0331', 'spec-2153-54212-0357', 'spec-1927-53321-0560', 'spec-1935-53387-0141', 'spec-2137-54206-0012', 'spec-1974-53430-0425', 'spec-2143-54184-0161', 'spec-2171-53557-0183', 'spec-2171-53557-0324', 'spec-2147-53491-0514', 'spec-2154-54539-0639', 'spec-1994-53845-0417', 'spec-2008-53473-0505', 'spec-2167-53889-0210', 'spec-2169-53556-0041', 'spec-2268-53682-0518', 'spec-2152-53874-0110', 'spec-2152-53874-0157', 'spec-1982-53436-0532', 'spec-1987-53765-0573', 'spec-2163-53823-0546', 'spec-2018-53800-0096', 'spec-2018-53800-0202', 'spec-2218-53816-0522', 'spec-2226-53819-0157', 'spec-2230-53799-0014', 'spec-2101-53858-0304', 'spec-2214-53794-0230', 'spec-1995-53415-0045', 'spec-2241-54169-0451', 'spec-2246-53767-0295', 'spec-2364-53737-0324', 'spec-2289-53708-0204', 'spec-2187-54270-0155', 'spec-2234-53823-0301', 'spec-2009-53904-0640', 'spec-2292-53713-0389', 'spec-2135-53827-0450', 'spec-2245-54208-0594', 'spec-2025-53431-0571', 'spec-2353-53794-0526', 'spec-2144-53770-0359', 'spec-2144-53770-0550', 'spec-2148-54526-0212', 'spec-2034-53466-0618', 'spec-2237-53828-0259', 'spec-2237-53828-0382', 'spec-2093-53818-0320', 'spec-2224-53815-0337', 'spec-2347-53757-0501', 'spec-2350-53765-0146', 'spec-2227-53820-0445', 'spec-2109-53468-0261', 'spec-2267-53713-0082', 'spec-2275-53709-0316', 'spec-2275-53709-0349', 'spec-2280-53680-0176', 'spec-2125-53795-0471', 'spec-2429-53799-0069', 'spec-2371-53762-0356', 'spec-2356-53786-0172', 'spec-2202-53566-0001', 'spec-2205-53793-0297', 'spec-2097-53491-0355', 'spec-2419-54139-0611', 'spec-2427-53815-0117', 'spec-2433-53820-0275', 'spec-2290-53727-0527', 'spec-2293-53730-0005', 'spec-2359-53826-0222', 'spec-2212-53789-0156', 'spec-2109-53468-0607', 'spec-2356-53786-0466', 'spec-2229-53823-0031', 'spec-2119-53792-0135', 'spec-2122-54178-0393', 'spec-2504-54179-0371', 'spec-2372-53768-0508', 'spec-2483-53852-0254', 'spec-2519-54570-0212', 'spec-2479-54174-0510', 'spec-2482-54175-0389', 'spec-2365-53739-0497', 'spec-2426-53795-0015', 'spec-2167-53889-0598', 'spec-2498-54169-0212', 'spec-2500-54178-0063', 'spec-2217-53794-0352', 'spec-2364-53737-0435', 'spec-2364-53737-0441', 'spec-2364-53737-0546', 'spec-2513-54141-0309', 'spec-2273-53709-0464', 'spec-2407-53771-0068', 'spec-2549-54523-0248', 'spec-2232-53827-0062', 'spec-2571-54055-0543', 'spec-2522-54570-0387', 'spec-2208-53880-0306', 'spec-2425-54139-0203', 'spec-2429-53799-0383', 'spec-2579-54068-0388', 'spec-2527-54569-0571', 'spec-2529-54585-0291', 'spec-2529-54585-0414', 'spec-2483-53852-0399', 'spec-2483-53852-0572', 'spec-2430-53815-0117', 'spec-2613-54481-0507', 'spec-2493-54115-0164', 'spec-2234-53823-0182', 'spec-2587-54138-0565', 'spec-2291-53714-0567', 'spec-2480-53851-0136', 'spec-2548-54152-0200', 'spec-2522-54570-0226', 'spec-2344-53740-0444', 'spec-2478-54097-0585', 'spec-2489-53857-0209', 'spec-2496-54178-0241', 'spec-2505-53856-0213', 'spec-2660-54504-0419', 'spec-2577-54086-0315', 'spec-2499-54176-0080', 'spec-2495-54175-0057', 'spec-2271-53726-0443', 'spec-2278-53711-0411', 'spec-2510-53877-0560', 'spec-2516-54241-0245', 'spec-2505-53856-0072', 'spec-2515-54180-0507', 'spec-2527-54569-0630', 'spec-2586-54169-0585', 'spec-2524-54568-0146', 'spec-2286-53700-0472', 'spec-2642-54232-0571', 'spec-2322-53727-0365', 'spec-2606-54154-0474', 'spec-2608-54474-0124', 'spec-2611-54477-0634', 'spec-2520-54584-0205', 'spec-2582-54139-0510', 'spec-2428-53801-0393', 'spec-2528-54571-0513', 'spec-2361-53762-0304', 'spec-2567-54179-0190', 'spec-2770-54510-0583', 'spec-2588-54174-0369', 'spec-2574-54084-0323', 'spec-2583-54095-0384', 'spec-2585-54097-0334', 'spec-2592-54178-0330', 'spec-2601-54144-0342', 'spec-2776-54554-0372', 'spec-2780-54557-0245', 'spec-2618-54506-0623', 'spec-2481-54086-0081', 'spec-2376-53770-0030', 'spec-2745-54231-0004', 'spec-2587-54138-0509', 'spec-2398-53768-0241', 'spec-2592-54178-0449', 'spec-2574-54084-0620', 'spec-2576-54086-0384', 'spec-2421-54153-0494', 'spec-2779-54540-0261', 'spec-2742-54233-0048', 'spec-2604-54484-0439', 'spec-2606-54154-0383', 'spec-2612-54480-0516', 'spec-2762-54533-0008', 'spec-2661-54505-0036', 'spec-2646-54479-0622', 'spec-2864-54467-0595', 'spec-2769-54527-0135', 'spec-2785-54537-0284', 'spec-2645-54477-0042', 'spec-2645-54477-0070', 'spec-2654-54231-0527', 'spec-2580-54092-0541', 'spec-2643-54208-0236', 'spec-2489-53857-0154', 'spec-2480-53851-0364', 'spec-2488-54149-0396', 'spec-2488-54149-0535', 'spec-2653-54230-0380', 'spec-2748-54234-0552', 'spec-2729-54419-0548', 'spec-2499-54176-0604', 'spec-2501-54084-0470', 'spec-2661-54505-0382', 'spec-2663-54234-0324', 'spec-2881-54502-0550', 'spec-2742-54233-0385', 'spec-2742-54233-0565', 'spec-2508-53875-0615', 'spec-2510-53877-0330', 'spec-2744-54272-0182', 'spec-2926-54625-0136', 'spec-2748-54234-0343', 'spec-2760-54506-0362', 'spec-2756-54508-0185', 'spec-2782-54592-0293', 'spec-2772-54529-0566', 'spec-2970-54589-0280', 'spec-2777-54554-0633', 'spec-2763-54507-0169', 'spec-2788-54553-0205', 'spec-2768-54265-0449', 'spec-2770-54510-0330', 'spec-2773-54533-0362', 'spec-2794-54537-0201', 'spec-2759-54534-0475', 'spec-2647-54495-0242', 'spec-2603-54479-0322', 'spec-2767-54243-0015', 'spec-2865-54503-0304', 'spec-2954-54561-0037', 'spec-2918-54554-0391', 'spec-2927-54621-0442', 'spec-2877-54523-0344', 'spec-2747-54233-0293', 'spec-2750-54242-0045', 'spec-2882-54498-0236', 'spec-2788-54553-0451', 'spec-2926-54625-0538', 'spec-2947-54533-0059', 'spec-2974-54592-0526', 'spec-2880-54509-0095', 'spec-2911-54631-0344', 'spec-2960-54561-0379', 'spec-3239-54888-0207', 'spec-2757-54509-0068', 'spec-2764-54535-0172', 'spec-2969-54586-0054', 'spec-3242-54889-0544', 'spec-2974-54592-0360', 'spec-2881-54502-0455', 'spec-2956-54525-0457', 'spec-2956-54525-0606', 'spec-3783-55246-0123', 'spec-3317-54908-0565', 'spec-3587-55182-0148', 'spec-3767-55214-0459', 'spec-3691-55274-0169', 'spec-3925-55338-0501', 'spec-3615-56544-0376', 'spec-3848-55647-0395', 'spec-3952-55330-0899', 'spec-3961-55654-0019', 'spec-3825-55533-0627', 'spec-3801-55509-0539', 'spec-3828-55539-0803', 'spec-3973-55323-0647', 'spec-3830-55574-0381', 'spec-3840-55574-0619', 'spec-3663-55176-0387', 'spec-3932-55337-0561', 'spec-3857-55272-0747', 'spec-3864-55649-0861', 'spec-3942-55338-0939', 'spec-3857-55272-0595', 'spec-4181-55685-0071', 'spec-3968-55590-0761', 'spec-4020-55332-0027', 'spec-3926-55327-0774', 'spec-4006-55328-0874', 'spec-4296-55499-0913', 'spec-4321-55504-0583', 'spec-4093-55475-0275', 'spec-4189-55679-0519', 'spec-4317-55480-0057', 'spec-4371-55830-0779', 'spec-4454-55536-0235', 'spec-4297-55806-0323', 'spec-4490-55629-0573', 'spec-4478-55600-0939', 'spec-4276-55505-0279', 'spec-4351-55538-0807', 'spec-4379-55881-0577', 'spec-4424-55532-0766', 'spec-4485-55836-0407', 'spec-4232-55447-0200', 'spec-4570-55623-0245', 'spec-4532-55559-0651', 'spec-4548-55565-0823', 'spec-4320-55894-0831', 'spec-4445-55869-0993', 'spec-4446-55589-0041', 'spec-4698-55623-0355', 'spec-4531-55563-0115', 'spec-4557-55588-0819', 'spec-4724-55742-0563', 'spec-4602-55644-0075', 'spec-4537-55806-0899', 'spec-4483-55587-0301', 'spec-4706-55705-0119', 'spec-4538-55860-0881', 'spec-4658-55592-0859', 'spec-4655-55620-0387', 'spec-4627-55626-0187', 'spec-4724-55742-0101', 'spec-4775-55708-0842', 'spec-4733-55649-0430', 'spec-4770-55928-0857', 'spec-4862-55685-0045', 'spec-4895-55708-0827', 'spec-4977-56047-0378', 'spec-4764-55646-0337', 'spec-4760-55656-0232', 'spec-4763-55869-0183', 'spec-4862-55685-0295', 'spec-4770-55928-0407', 'spec-4691-55651-0745', 'spec-4801-55653-0952', 'spec-4903-55927-0800', 'spec-4966-55712-0793', 'spec-5165-56063-0345', 'spec-4783-55652-0053', 'spec-5177-56245-0697', 'spec-5200-56091-0681', 'spec-5045-56181-0997', 'spec-5184-56352-0557', 'spec-5184-56352-0343', 'spec-5128-55912-0299', 'spec-5301-55987-0625', 'spec-5391-56000-0177', 'spec-5316-55955-0133', 'spec-5398-56011-0062', 'spec-5363-55956-0305', 'spec-5371-55976-0568', 'spec-5398-56011-0229', 'spec-5466-56033-0441', 'spec-5649-55912-0445', 'spec-5488-56013-0339', 'spec-5177-56245-0742', 'spec-5362-56017-0341', 'spec-5444-56038-0595', 'spec-5450-55986-0047', 'spec-5410-56016-0303', 'spec-5466-56033-0198', 'spec-5427-56001-0289', 'spec-5436-56015-0463', 'spec-5321-55945-0361', 'spec-5495-55896-0529', 'spec-5777-56280-0893', 'spec-5866-56035-0797', 'spec-5808-56325-0383', 'spec-5344-55924-0017', 'spec-5741-55980-0039', 'spec-5767-56245-0126', 'spec-5865-56067-0822', 'spec-5785-56269-0387', 'spec-5963-56191-0105', 'spec-5852-56034-0153', 'spec-5438-56002-0426', 'spec-5888-56041-0705', 'spec-5796-56274-0105', 'spec-6140-56189-0595', 'spec-5743-56011-0352', 'spec-6005-56090-0873', 'spec-5880-56042-0587', 'spec-6050-56089-0885', 'spec-6145-56266-0163', 'spec-6385-56356-0119', 'spec-6281-56295-0829', 'spec-6293-56561-0573', 'spec-6167-56189-0231', 'spec-6278-56266-0215', 'spec-6114-56209-0437', 'spec-6125-56273-0741', 'spec-6402-56334-0033', 'spec-6431-56311-0161', 'spec-6308-56215-0641', 'spec-6488-56364-0823', 'spec-6176-56264-0143', 'spec-6380-56340-0182', 'spec-6409-56306-0367', 'spec-6439-56358-0552', 'spec-6461-56329-0529', 'spec-6475-56337-0385', 'spec-6519-56566-0861', 'spec-6622-56365-0320', 'spec-6458-56274-0655', 'spec-6495-56339-0478', 'spec-6664-56383-0845', 'spec-6688-56412-0259', 'spec-6730-56425-0898', 'spec-6620-56368-0337', 'spec-6813-56419-0498', 'spec-6456-56339-0275', 'spec-6671-56388-0473', 'spec-6706-56385-0307', 'spec-6727-56369-0267', 'spec-6976-56448-0597', 'spec-6821-56396-0363', 'spec-6716-56401-0835', 'spec-7049-56570-0467', 'spec-6649-56364-0285', 'spec-6728-56426-0277', 'spec-6739-56393-0901', 'spec-7132-56565-0369', 'spec-7128-56567-0023', 'spec-7261-56603-0267', 'spec-7132-56565-0949', 'spec-7236-56605-0764', 'spec-6790-56430-0085', 'spec-7295-57067-0731', 'spec-7315-56685-0165', 'spec-7153-56904-0984', 'spec-7564-56804-0220', 'spec-7604-56947-0944', 'spec-8057-57190-0321', 'spec-8068-57185-0168', 'spec-8199-57428-0845', 'spec-8208-57430-0244', 'spec-7870-57016-0041', 'spec-7705-57332-0036', 'spec-3851-55302-0325', 'spec-4075-55352-0777', 'spec-4871-55928-0299', 'spec-6290-56238-0843']
    synthetic_runs = ['synthetic1', 'synthetic2', 'synthetic3', 'synthetic4', 'synthetic5', 'synthetic6', 'synthetic7', 'synthetic8']
    # Set which galaxy to run
    #rungal = 'AOS2015'
    #rungal = 'test'
    rungal = 'spec-0266-51630-0407'
    #rungal = 'synthetic'

    # First, remove 'all_output' file containing old output, if it exists
    if os.path.exists('all_output'):
        os.remove('all_output')
    # Otherwise, create 'all_output' to store MCMC output
    outfile = open('all_output', 'w')
    outfile.write(
        'Object y+ y+_p y+_m dens dens_p dens_m aHe aHe_p aHe_m tauHe tauHe_p tauHe_m temp temp_p temp_m cHb cHb_p cHb_m aH aH_p aH_m xi xi_p xi_m\n')
    outfile.close()

    if rungal == 'AOS2015':
        # Run MCMC on all galaxies
        galfail = []
        for gal in AOS2015:
            try:
                print ('Working on', gal)
                MCMCgal(gal)
            except IOError:
                print('ERROR :: The following galaxy data could not be found: {0:s}'.format(gal))
                galfail += [gal]
            except TypeError:
                print('ERROR :: The following galaxy is not known: {0:s}'.format(gal))
                galfail += [gal]
            except ValueError:
                print('ERROR :: The following galaxy failed: {0:s}'.format(gal))
                galfail += [gal]
        print('The following galaxies failed:\n' + '\n'.join(galfail))

    elif rungal == 'ours':
        # First, remove the file containing old output
        if os.path.exists('all_output'):
            os.remove('all_output')
        outfile = open('all_output', 'w')
        outfile.write('Object y+ y+_p y+_m dens dens_p dens_m aHe aHe_p aHe_m tauHe tauHe_p tauHe_m temp temp_p temp_m cHb cHb_p cHb_m aH aH_p aH_m xi xi_p xi_m\n')
        outfile.close()
        # Run MCMC on all galaxies
        galfail = []
        for gal in ours:
            try:
                print ('Working on', gal)
                MCMCgal(gal)
            except IOError:
                print('ERROR :: The following galaxy data could not be found: {0:s}'.format(gal))
                galfail += [gal]
            except TypeError:
                print('ERROR :: The following galaxy is not known: {0:s}'.format(gal))
                galfail += [gal]
            except ValueError:
                print('ERROR :: The following galaxy failed: {0:s}'.format(gal))
                galfail += [gal]
        print('The following galaxies failed:\n' + '\n'.join(galfail))

    elif rungal == 'SDSS':
        galfail = []
        for gal in SDSS:
            try:
                print ('Working on', gal)
                MCMCgal(gal)
            except IOError:
                print('ERROR :: The following galaxy data could not be found: {0:s}'.format(gal))
                galfail += [gal]
            except TypeError:
                print('ERROR :: The following galaxy is not known: {0:s}'.format(gal))
                galfail += [gal]
            except ValueError:
                print('ERROR :: The following galaxy failed: {0:s}'.format(gal))
                galfail += [gal]
        print('The following galaxies failed:\n' + '\n'.join(galfail))

    elif rungal == 'synthetic':
        synfail = []
        for syn in synthetic_runs:
            try:
                MCMCgal(syn)
            except IOError:
                print('ERROR :: The following galaxy data could not be found: {0:s}'.format(syn))
                synfail += [syn]
            except TypeError:
                print('ERROR :: The following galaxy is not known: {0:s}'.format(syn))
                synfail += [syn]
            except ValueError:
                print('ERROR :: The following galaxy failed: {0:s}'.format(syn))
                synfail += [syn]
        print('The following galaxies failed:\n' + '\n'.join(synfail))

    else:
        if rungal in AOS2015 or rungal in ours or rungal in SDSS or rungal == 'test':
            MCMCgal(rungal)
        else:
            print('Invalid Galaxy name. Select one of the following:\n' + '\n'.join(names))
