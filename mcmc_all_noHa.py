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

        galdict = galaxy.load_ours_noHaHb(self.galaxyname) # Our LRIS sample assuming no Ha/Hb

        self._full_tbl = galdict['full_tbl']
        self._T_OIII = galdict['T_OIII']

        # Read in measured data (wavelength, flux ratios, and EWs)
        if 10833.306 in self._full_tbl['Wavelength']:
            self._flux_ratios = self._full_tbl[:-1]  # Ignore the entry for P-gamma for MCMC'
        else:
            self._flux_ratios = self._full_tbl

        print (self._flux_ratios)
        # Names of wavelenghts of interest for MCMC
        self._hydrogen_lines = np.array([10941.082, 6564.612, 4862.721, 4341.684, 4102.891, 3890.166])  # Pa-g, Ha, Hb, Hg, Hd, H8
        self._helium_lines = np.array([10833.306, 7067.198, 6679.994, 5877.299, 5017.079, 4472.755, 4027.328, 3890.151])

        self._allowed_lines = np.sort(np.concatenate((self._hydrogen_lines, self._helium_lines)))[1:]  # [1:] to remove the duplicate ~3890 wavelength

        self._emis_lines = self._flux_ratios['Wavelength']

        # Measured data from spectra
        self._EWs_meas = np.array(self._flux_ratios['EW'])

        self._y = np.array(self._flux_ratios['Flux Ratio'])  # F(lambda) / F(H-beta)
        try:
#            self._y_error = np.array(self._flux_ratios['Flux Ratio Errors'])
            self._y_error =  np.sqrt( np.array(self._flux_ratios['Flux Ratio Errors'])**2. + ( 0.09 * np.array(self._flux_ratios['Flux Ratio']) )**2. )

        except:
            self._y_error = np.sqrt(np.array(self._flux_ratios['Flux Ratio'] * 0.002) ** 2. + np.array(self._flux_ratios['Flux Ratio'] * 0.02) ** 2.)

        self._x = np.zeros(self._y.size)

        # Range of values for 8 parameters: y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, xi/n_HI
        self._min_y_plus, self._max_y_plus = 0.05, 0.1  # fraction of singly ionized He; y+ = He+/H+
        self._min_temp, self._max_temp = 10000, 22000  # electron temperature (K)
        self._min_log_dens, self._max_log_dens = 0, 3  # log10(electron density) (cm^-3)
        self._min_c_Hb, self._max_c_Hb = 0, 0.5  # reddening
        self._min_a_H, self._max_a_H = 0, 10  # underlying stellar H absorption (Angstroms)
        self._min_a_He, self._max_a_He = 0, 4  # underlying stellar HeI absorption (Angstroms)
        self._min_tau_He, self._max_tau_He = 0, 5  # optical depth; range of values from Izotov & Thuan 2010
        self._min_log_xi, self._max_log_xi = -6, -0.0969  # equals to xi=0-0.8; ratio of neutral hydrogen to singly ionized hydrogen densities; xi = n(HI)/n(HII)
        self.mcmc_steps()

    def __call__(self, theta, x, y, yerr):
        return self.lnprob(theta, x, y, yerr)

    # Set up MCMC
    def get_model(self, theta):
        y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi = theta

        # Reparameterize density and neutral hydrogen input
        dens = 10 ** log_dens
        xi = 10 ** log_xi

        # Take into account error on EW(Hb) by perturbing EW(Hb) by its measured EW error
        try:
            EW_Hb = np.random.normal(
                self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0],
                self._flux_ratios['EW Errors'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0])
        except:
            EW_Hb = np.random.normal(
                self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0],
                0.1 * self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 4862.721)[0]][0])

        # Continuum level ratio; Eq. 2.4 of AOS2012
        h = self._y * EW_Hb / self._EWs_meas  # relative to H-beta; i.e., h(lambda) / h(Hbeta)

        # Model flux
        model_flux = np.zeros(self._y.size)  # emission lines we want to model

        # Some values, calculated at Hbeta, for later use in model flux
        collisional_to_recomb_Hbeta = mfr.hydrogen_collision_to_recomb(xi, self._hydrogen_lines[2], temp, method='A2002')
        f_lambda_at_Hbeta = mfr.f_lambda_avg_interp(self._hydrogen_lines[2])
        # AHbeta_Av = mfr.reddening_coefficient(hydrogen_lines[2]) # CCM 1989 reddening curve

        done_3889 = False
        for w in range(len(self._emis_lines)):
            # Determine if working with hydrogen or helium line; within 3.5 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
            nearest_wave = self._allowed_lines[np.where(np.abs(self._allowed_lines - self._emis_lines[w]) < 3.5)[0]][0]
            # The above line is redundant for my input waves, but allows for cases where emis_lines[w] is some other array, say waves_of_interest[w],
            # and not exactly at the wavelengths given in the allowed_lines array (which is concatenated from arrays hydrogen_lines and helium_lines)

            # Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166) and P-gamma
            if nearest_wave in self._hydrogen_lines and nearest_wave != 3890.166 and nearest_wave != 10941.082:  # and nearest_wave != 4862.721:
                line_species = 'hydrogen'

                # Separate -- if HI line is Ha, and we should generate its flux ratio w.r.t itself, i.e. flux = 1.0, since we don't have red/blue side info
                if nearest_wave == 6564.612:
                    flux = 1.0

                # Otherwise, do F(HI)/F(Hb) as usual
                else:
                    emissivity_ratio = mfr.hydrogen_emissivity_S2018(self._emis_lines[w], temp, dens)
                    a_H_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_H, ion=line_species)
                    collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, self._emis_lines[w], temp, method='A2002')
                    reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.

                    #            flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
                    #                ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    #                10**-(reddening_function * c_Hb)
                    # Reparameterization of flux to use the continuum level
                    flux = (emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) *
                            10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_H_at_wave / EW_Hb) * (h[w]))

            # Any HeI line besides the blended HeI+H8 line (HeI at 3890.151)
            elif nearest_wave in self._helium_lines and nearest_wave != 3890.151 and nearest_wave != 10833.306:

                # Separate -- if HeI line is on 'red' side, need a HeI/Ha ratio
                if nearest_wave >= 5877.299:
                    # First, theoretical F(Halpha)/F(Hbeta) ratio (e.g., Eq. 2.3 of AOS11)
                    line_species = 'hydrogen'

                    emissivity_ratio = mfr.hydrogen_emissivity_S2018(6564.612, temp, dens)
                    a_H_at_wave = mfr.stellar_absorption(6564.612, a_H, ion=line_species)
                    collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, 6564.612, temp, method='A2002')
                    reddening_function = (mfr.f_lambda_avg_interp(6564.612) / f_lambda_at_Hbeta) - 1.

                    # EW_Ha = self._full_tbl['EW'][np.where(self._full_tbl['Wavelength'] == 6564.612)[0][0]]
                    # Take into account error on EW(Ha) by perturbing EW(Ha) by its measured EW error
                    try:
                        EW_Ha = np.random.normal(
                            self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 6564.612)[0]][0],
                            self._flux_ratios['EW Errors'][np.where(self._flux_ratios['Wavelength'] == 6564.612)[0]][0])
                    except:
                        EW_Ha = np.random.normal(
                            self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 6564.612)[0]][0],
                            0.1 * self._flux_ratios['EW'][np.where(self._flux_ratios['Wavelength'] == 6564.612)[0]][0])

                    Ha_to_Hb_flux = emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) * \
                                    10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb)) / ((EW_Ha + a_H_at_wave) / (EW_Ha))

                    # Now, F(HeI)/F(Hb) ratio
                    line_species = 'helium'

                    emissivity_ratio = mfr.helium_emissivity_PFSD2012(self._emis_lines[w], temp, dens)
                    a_He_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_He, ion=line_species)
                    optical_depth_at_wave = mfr.optical_depth_function(self._emis_lines[w], temp, dens, tau_He)
                    reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.

                    HeI_to_Hb_flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                                      10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_He_at_wave / EW_Hb) * (h[w] * Ha_to_Hb_flux))

                    # Want to get theoretical F(HeI)/F(Ha) to match that in the input table of flux_ratios -- can do this by ( F(HeI)/F(Hbeta) ) / ( F(Halpha)/F(Hbeta) )!
                    flux = HeI_to_Hb_flux / Ha_to_Hb_flux

                # Otherwise, do F(HeI)/F(Hb) as usual
                else:
                    line_species = 'helium'

                    emissivity_ratio = mfr.helium_emissivity_PFSD2012(self._emis_lines[w], temp, dens)
                    a_He_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_He, ion=line_species)
                    optical_depth_at_wave = mfr.optical_depth_function(self._emis_lines[w], temp, dens, tau_He)
                    reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.

                    flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                            10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_He_at_wave / EW_Hb) * (h[w]))

            # The blended HeI+H8 line
            elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
                if done_3889: continue
                done_3889 = True
                reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.
                # reddening_function = ( mfr.reddening_coefficient(self._emis_lines[w]) / AHbeta_Av ) - 1. # CCM 1989 reddening curve

                # HeI 3890.151 contribution:
                line_species = 'helium'

                emissivity_ratio = mfr.helium_emissivity_PFSD2012(self._emis_lines[w], temp, dens)
                a_He_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_He, ion=line_species)
                optical_depth_at_wave = mfr.optical_depth_function(self._emis_lines[w], temp, dens, tau_He)

                flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                        10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / EW_Hb)) - (
                                   (a_He_at_wave / EW_Hb) * (h[w]))

                # H8 contribution:
                line_species = 'hydrogen'

                emissivity_ratio = mfr.hydrogen_emissivity_S2018(self._emis_lines[w], temp, dens)
                a_H_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_H, ion=line_species)
                collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, self._emis_lines[w], temp, method='A2002')

                flux += (emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) *
                         10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - (
                                    (a_H_at_wave / EW_Hb) * (h[w]))

            # Infrared HeI10830 line
            elif nearest_wave == 10833.306:
                # Theoretical F(Pg)/F(Hb) ratio, aka the 'model-dependent scaling ratio' (e.g., Eq. 2.3 of AOS11)
                line_species = 'hydrogen'

                emissivity_ratio = mfr.hydrogen_emissivity_S2018(10941.082, temp, dens)  # hard-coded Pg wavelength; could also be hydrogen_lines[0]
                a_H_at_wave = mfr.stellar_absorption(10941.082, a_H, ion=line_species)
                collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, self._hydrogen_lines[0], temp, method='A2002')
                reddening_function = (mfr.f_lambda_avg_interp(10941.082) / f_lambda_at_Hbeta) - 1.

                #EW_Pg = self._full_tbl[np.where(self._full_tbl['Wavelength'] == 10941.082)[0][0]]['EW']
                # Take into account error on EW(Pg) by perturbing EW(Pg) by its measured EW error
                try:
                    EW_Pg = np.random.normal(
                        self._full_tbl['EW'][np.where(self._full_tbl['Wavelength'] == 10941.082)[0]][0],
                        self._full_tbl['EW Errors'][np.where(self._full_tbl['Wavelength'] == 10941.082)[0]][0])
                except:
                    EW_Pg = np.random.normal(
                        self._full_tbl['EW'][np.where(self._full_tbl['Wavelength'] == 10941.082)[0]][0],
                        0.1 * self._full_tbl['EW'][np.where(self._full_tbl['Wavelength'] == 10941.082)[0]][0])

                Pg_to_Hb_flux = emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) * \
                                ((EW_Hb + a_H) / (EW_Hb)) / ((EW_Pg + a_H_at_wave) / (EW_Pg)) * 10 ** -(reddening_function * c_Hb)

                # Theoretical F(HeI10830)/F(Hbeta) ratio
                line_species = 'helium'

                emissivity_ratio = mfr.helium_emissivity_PFSD2012(self._emis_lines[w], temp, dens)
                a_He_at_wave = mfr.stellar_absorption(self._emis_lines[w], a_He, ion=line_species)
                optical_depth_at_wave = mfr.optical_depth_function(self._emis_lines[w], temp, dens, tau_He)
                reddening_function = (mfr.f_lambda_avg_interp(self._emis_lines[w]) / f_lambda_at_Hbeta) - 1.
                # reddening_function = ( mfr.reddening_coefficient(self._emis_lines[w]) / AHbeta_Av ) - 1. # CCM 1989 reddening curve

                # The way h is defined above and given the format of the input fluxes gives ( F(HeI10830)/F(Pg) ) * ( EW(Hb) / EW(HeI10830) ) here; must be multiplied by the calculated
                # theoretical F(Pg)/F(Hb) ratio from above to get the HeI10830 to Hbeta continuum level ratio, which is the definition of h, from Eq. 2.4 of AOS2012
                HeI10830_to_Hb_flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                                       10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - (
                                                  (a_He_at_wave / EW_Hb) * (h[w] * Pg_to_Hb_flux))

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

        print('Running MCMC on', self.galaxyname, '...')
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
        prenams = ['y+', 'temperature', '$log(n_{e})$', 'c(H\\beta)', '$a_{H}$', '$a_{He}$', '$\\tau_{He}',
                   '$log(\\xi)$']
        input_vals = np.array([0.08, 18000, 2, 0.1, 1.0, 1.0, 1.0, -2])  # Input parameters for fake spectra

        print ('Best parameter values:')
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
        print ('\n Input parameter values:')
        print (input_vals)

        dens = 10.0 ** (samples[:, 2])
        v = np.percentile(dens, [16, 50, 84])
        dens_mcmc = (v[1], v[2] - v[1], v[1] - v[0],)
        # Save some output
        allpars = np.hstack(
            (y_plus_mcmc, dens_mcmc, a_He_mcmc, tau_He_mcmc, temp_mcmc, c_Hb_mcmc, a_H_mcmc, log_xi_mcmc))
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
    ours = ['J0000p3052A', 'J0000p3052B', 'J0018p2345', 'J0118p3512', 'J0220p2044A',
             'J0452m0541', 'J0743p4807', 'J0757p4750', 'J0943p3326', 'J1044p6306',
             'J1204p5259', 'J1214p1245', 'J1322p5425', 'J1425p4441', 'J1655p6337',
             'J1705p3527', 'J1732p4452', 'J1757p6454', 'J2030m1343', 'J2213p1722', 'J2319p1616',
            'J2230m0531', 'J2339p3230', 'KJ2', 'KJ29', 'KJ5', 'KJ5B', 'KJ97']
    # Removed ['J0140p2951', 'J0214m0835', 'J0220p2044B', 'J1414m0208'] for too many missing lines; Leo P is in mcmc_all.py
    # J1214p1245, J1425p4441, J1757p6454, J2339p3230 have EW(Hb) < 50
    # Set which galaxy to run
    rungal = 'test'
    #rungal = 'ours'

    if rungal == 'ours':
        # First, remove the file containing old output
        if os.path.exists('all_output'):
            os.remove('all_output')
        outfile = open('all_output', 'w')
        outfile.write(
            'Object y+ y+_p y+_m dens dens_p dens_m aHe aHe_p aHe_m tauHe tauHe_p tauHe_m temp temp_p temp_m cHb cHb_p cHb_m aH aH_p aH_m xi xi_p xi_m\n')
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

    else:
        if rungal in ours or rungal == 'test':
            MCMCgal(rungal)
        else:
            print('Invalid Galaxy name. Select one of the following:\n' + '\n'.join(names))
