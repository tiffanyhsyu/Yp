###########
# Updates #
###########
# 2019-06-25: Created to start and test MCMC without Ha/Hb flux ratios (the 'arm' setup)

import emcee
import time
import os
import numpy as np
import model_flux_ratio as mfr
import pdb

from astropy.table import Table

# Read in measured data (wavelength, flux ratios, and EWs)
full_tbl = Table.read(os.getcwd() + '/test_no_Ha', format='ascii', delimiter=' ')
# NIR
flux_ratios = full_tbl[:-1]  # Ignore the last entry, assumed to be for P-gamma, for MCMC
# Optical
# flux_ratios = full_tbl

# Names of wavelengths of interest for MCMC
y_names = ['HeI+H83890', 'HeI4027', 'Hd', 'Hg', 'HeI4472', 'Hb', 'HeI5017', 'HeI5877', 'Ha', 'HeI6679', 'HeI7067', 'HeI10830']

# Hydrogen and Helium lines of possible interest for MCMC
hydrogen_lines = np.array([10941.082, 6564.612, 4862.721, 4341.684, 4102.891, 3890.166]) # Pg, Ha, Hb, Hg, Hd, H8
helium_lines = np.array([10833.306, 7067.198, 6679.994, 5877.299, 5017.079, 4472.755, 4027.328, 3890.151])
allowed_lines = np.sort(np.concatenate((hydrogen_lines, helium_lines)))[1:] # Remove the duplicate H8, HeI3890 line

# Emission lines actually measured and to be used in the MCMC
emis_lines = np.array(flux_ratios['Wavelength'])

print (flux_ratios)

# Measured data from spectra
EWs_meas = np.array(flux_ratios['EW'])

y = np.array(flux_ratios['Flux Ratio'])  # F(lambda) / F(H-beta)
y_error = np.array(flux_ratios['Flux Ratio'] * 0.1)  # For test_output flux, where we have no EW errors
# y_error = np.array(flux_ratios['Flux Ratio Errors'])
x = np.zeros(y.size)

# Range of values for 8 parameters: y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, xi/n_HI
min_y_plus, max_y_plus = 0.05, 0.1  # fraction of singly ionized He; y+ = He+/H+
min_temp, max_temp = 10000, 22000  # restricted from 5000, 25000 to match finemesh HeI emissivities # electron temperature (K)
min_log_dens, max_log_dens = 0, 4  # log10(electron density) (cm^-3)
min_c_Hb, max_c_Hb = 0, 0.5  # reddening
min_a_H, max_a_H = 0, 10  # underlying stellar H absorption (Angstroms)
min_a_He, max_a_He = 0, 4  # underlying stellar HeI absorption (Angstroms)
min_tau_He, max_tau_He = 0, 5  # optical depth; range of values from Izotov & Thuan 2010
min_log_xi, max_log_xi = -6, -0.0969  # equals to xi=0-0.8; ratio of neutral hydrogen to singly ionized hydrogen densities; xi = n(HI)/n(HII)


# Set up MCMC
def get_model(theta):
    #    y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, n_HI = theta
    y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi = theta

    # Reparameterize density and neutral hydrogen input
    dens = 10 ** log_dens
    xi = 10 ** log_xi

    # Take into account error on EW(Hb) by perturbing EW(Hb) by some EW error
    EW_Hb = np.random.normal(flux_ratios['EW'][np.where(flux_ratios['Wavelength'] == 4862.721)[0]][0],
                             0.1 * flux_ratios['EW'][np.where(flux_ratios['Wavelength'] == 4862.721)[0]][0]) # flux_ratios['EW Errors'][np.where(flux_ratios['Wavelength'] == 4862.721)[0]][0])

    # Continuum level ratio; Eq. 2.4 of AOS2012
    h = y * EW_Hb / EWs_meas  # relative to H-beta; i.e., h(lambda) / h(H-beta)

    # Model flux
    model_flux = np.zeros(y.size)  # emission line fluxes we want to model

    # Some values, calculated at Hbeta, for later use in model flux
    collisional_to_recomb_Hbeta = mfr.hydrogen_collision_to_recomb(xi, 4862.721, temp, method='A2002')
    f_lambda_at_Hbeta = mfr.f_lambda_avg_interp(4862.721)

    #### Should emis_lines here be flux_ratios['Wavelengths']??
    for w in range(len(emis_lines)):
        # Determine if working with hydrogen or helium line; within 3.5 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
        nearest_wave = allowed_lines[np.where(np.abs(allowed_lines - emis_lines[w]) < 3.5)[0]][0]

        # Any HI line besides the blended HeI+H8 line (H8 at 3890.166) and Pg
        if nearest_wave in hydrogen_lines and nearest_wave != 3890.166 and nearest_wave != 10941.082:
            line_species = 'hydrogen'

            # Separate -- if HI line is Ha, and we should generate its flux ratio w.r.t itself, i.e. flux = 1.0, since we don't have red/blue side info
            if nearest_wave == 6564.612:
                flux = 1.0

            # Otherwise, do F(HI)/F(Hb) as usual
            else:
                emissivity_ratio = mfr.hydrogen_emissivity_S2018(emis_lines[w], temp, dens)
                a_H_at_wave = mfr.stellar_absorption(emis_lines[w], a_H, ion=line_species)
                collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, emis_lines[w], temp, method='A2002')
                reddening_function = (mfr.f_lambda_avg_interp(emis_lines[w]) / f_lambda_at_Hbeta) - 1.

                #            flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
                #                ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                #                10**-(reddening_function * c_Hb)
                # Reparameterization of flux to use the continuum level
                flux = (emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) *
                        10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_H_at_wave / EW_Hb) * (h[w]))

        # Any HeI line besides the blended HeI+H8 line (HeI at 3890.151) and NIR HeI10830 line
        elif nearest_wave in helium_lines and nearest_wave != 3890.151 and nearest_wave != 10833.306:

            # Separate -- if HeI line is on 'red' side, need a HeI/Ha ratio
            if nearest_wave >= 5877.299:
                # First, theoretical F(Halpha)/F(Hbeta) ratio
                line_species = 'hydrogen'

                emissivity_ratio = mfr.hydrogen_emissivity_S2018(6564.612, temp, dens)
                a_H_at_wave = mfr.stellar_absorption(6564.612, a_H, ion=line_species)
                collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, 6564.612, temp, method='A2002')
                reddening_function = (mfr.f_lambda_avg_interp(6564.612) / f_lambda_at_Hbeta) - 1.

                EW_Ha = full_tbl['EW'][np.where(full_tbl['Wavelength'] == 6564.612)[0][0]]
                Ha_to_Hb_flux = emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) * \
                                10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb)) / ((EW_Ha + a_H_at_wave) / (EW_Ha))

                # Now, F(HeI)/F(Hb) ratio
                line_species = 'helium'

                emissivity_ratio = mfr.helium_emissivity_PFSD2012(emis_lines[w], temp, dens)
                a_He_at_wave = mfr.stellar_absorption(emis_lines[w], a_He, ion=line_species)
                optical_depth_at_wave = mfr.optical_depth_function(emis_lines[w], temp, dens, tau_He)
                reddening_function = (mfr.f_lambda_avg_interp(emis_lines[w]) / f_lambda_at_Hbeta) - 1.

                HeI_to_Hb_flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                                   10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_He_at_wave / EW_Hb) * (h[w] * Ha_to_Hb_flux))

                # Want to get theoretical F(HeI)/F(Ha) to match that in the input table of flux_ratios -- can do this by ( F(HeI)/F(Hbeta) ) / ( F(Halpha)/F(Hbeta) )!
                flux = HeI_to_Hb_flux / Ha_to_Hb_flux

            # Otherwise, do F(HeI)/F(Hb) as usual
            else:
                line_species = 'helium'

                emissivity_ratio = mfr.helium_emissivity_PFSD2012(emis_lines[w], temp, dens)
                a_He_at_wave = mfr.stellar_absorption(emis_lines[w], a_He, ion=line_species)
                optical_depth_at_wave = mfr.optical_depth_function(emis_lines[w], temp, dens, tau_He)
                reddening_function = (mfr.f_lambda_avg_interp(emis_lines[w]) / f_lambda_at_Hbeta) - 1.

                flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) * \
                        10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_He_at_wave / EW_Hb) * (h[w]))

        # The blended HeI+H8 line
        elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
            reddening_function = (mfr.f_lambda_avg_interp(emis_lines[w]) / f_lambda_at_Hbeta) - 1.

            # HeI 3890.151 contribution:
            line_species = 'helium'

            emissivity_ratio = mfr.helium_emissivity_PFSD2012(emis_lines[w], temp, dens)
            a_He_at_wave = mfr.stellar_absorption(emis_lines[w], a_He, ion=line_species)
            optical_depth_at_wave = mfr.optical_depth_function(emis_lines[w], temp, dens, tau_He)

            flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) * \
                    10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_He_at_wave / EW_Hb) * (h[w]))

            # H8 contribution:
            line_species = 'hydrogen'

            emissivity_ratio = mfr.hydrogen_emissivity_S2018(emis_lines[w], temp, dens)
            a_H_at_wave = mfr.stellar_absorption(emis_lines[w], a_H, ion=line_species)
            collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, emis_lines[w], temp, method='A2002')

            flux += (emissivity_ratio * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) *
                     10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_H_at_wave / EW_Hb) * (h[w]))

        # Infrared HeI10830 line
        elif nearest_wave == 10833.306:
            # First, theoretical F(Pg)/F(Hb) ratio, aka 'model-dependent scaling' factor
            line_species = 'hydrogen'

            emissivity_ratio = mfr.hydrogen_emissivity_S2018(10941.082, temp,
                                                             dens)  # hard-coded Pg wavelength; could also be hydrogen_lines[0]
            a_H_at_wave = mfr.stellar_absorption(10941.082, a_H, ion=line_species)
            collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, 10941.082, temp, method='A2002')
            reddening_function = (mfr.f_lambda_avg_interp(10941.082) / f_lambda_at_Hbeta) - 1.  # hard-coded Pg wavelength; could also be hydrogen_lines[0]

            # This ratio is all theoretical, hence no reformulation of the flux equation to use the continuum level here
            EW_Pg = full_tbl['EW'][np.where(full_tbl['Wavelength'] == 10941.082)[0][0]]
            Pg_to_Hb_flux = emissivity_ratio * ((EW_Hb + a_H) / (EW_Hb)) / ((EW_Pg + a_H_at_wave) / (EW_Pg)) * \
                            ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) * 10 ** -(reddening_function * c_Hb)

            # Theoretical F(HeI10830)/F(Hbeta) ratio
            line_species = 'helium'

            emissivity_ratio = mfr.helium_emissivity_PFSD2012(emis_lines[w], temp, dens)
            a_He_at_wave = mfr.stellar_absorption(emis_lines[w], a_He, ion=line_species)
            optical_depth_at_wave = mfr.optical_depth_function(emis_lines[w], temp, dens, tau_He)
            reddening_function = (mfr.f_lambda_avg_interp(emis_lines[w]) / f_lambda_at_Hbeta) - 1.

            # The way h is defined above and given the format of the input fluxes gives ( F(HeI10830)/F(Pg) ) * ( EW(Hb) / EW(HeI10830) ) here; must be multiplied by the calculated
            # theoretical F(Pg)/F(Hb) ratio from above to get the HeI10830 to Hbeta continuum level ratio, which is the definition of h, from Eq. 2.4 of AOS2012
            HeI10830_to_Hb_flux = (y_plus * emissivity_ratio * optical_depth_at_wave * (1 / (1 + collisional_to_recomb_Hbeta)) *
                                   10 ** -(reddening_function * c_Hb) * ((EW_Hb + a_H) / (EW_Hb))) - ((a_He_at_wave / EW_Hb) * (h[w] * Pg_to_Hb_flux))

            # Want to get theoretical F(HeI10830)/F(Pg) to match that in the input table of flux_ratios -- can do this by ( F(HeI10830)/F(Hbeta) ) / ( F(Hbeta)/F(Pg) )!
            flux = HeI10830_to_Hb_flux / Pg_to_Hb_flux

        else:
            print ('Check your input wavelength -- not a recognized hydrogen or helium line for MCMC analysis')
            pdb.set_trace()

        model_flux[w] = flux

    return model_flux


# Define the probability function as likelihood * prior.
def lnprior(theta):
    #    y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, n_HI = theta
    y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi = theta

    if min_y_plus <= y_plus <= max_y_plus and \
            min_temp <= temp <= max_temp and \
            min_log_dens <= log_dens <= max_log_dens and \
            min_c_Hb <= c_Hb <= max_c_Hb and \
            min_a_H <= a_H <= max_a_H and \
            min_a_He <= a_He <= max_a_He and \
            min_tau_He <= tau_He <= max_tau_He and \
            min_log_xi <= log_xi <= max_log_xi:
        return -(temp - 16500) ** 2 / (0.2 * 16500) ** 2  # 0.0

    return -np.inf


def lnlike(theta, x, y, yerr):
    model = get_model(theta)
    inv_sigma2 = 1.0 / yerr ** 2

    return -0.5 * (np.sum((y - model) ** 2 * inv_sigma2 - np.log(inv_sigma2)))


def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)

    if not np.isfinite(lp):
        return -np.inf

    return lp + lnlike(theta, x, y, yerr)


# Set up sampler
ndim, nwalkers = 8, 500

pos = [np.array([np.random.uniform(min_y_plus, max_y_plus),
                 np.random.uniform(min_temp, max_temp),
                 np.random.uniform(min_log_dens, max_log_dens),
                 np.random.uniform(min_c_Hb, max_c_Hb),
                 np.random.uniform(min_a_H, max_a_H),
                 np.random.uniform(min_a_He, max_a_He),
                 np.random.uniform(min_tau_He, max_tau_He),
                 np.random.uniform(min_log_xi, max_log_xi)]) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, y_error), threads=ndim)

print('Running MCMC...')
nmbr = 1000
a = time.time()
for i, result in enumerate(sampler.run_mcmc(pos, nmbr, rstate0=np.random.get_state())):
    if True:  # (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i) / nmbr))
print('Done!')
print((time.time() - a) / 60.0, 'mins')

print('Saving samples')
np.save('{0:s}_{1:d}walkers_{2:d}steps'.format('test_no_Ha', nwalkers, nmbr), sampler.chain)

print('Making plots')
burnin = int(0.8 * nmbr)

samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
# Names of 8 parameters and input 'true' parameter values
prenams = ['y+', 'temperature', '$log(n_{e})$', 'c(H\\beta)', '$a_{H}$', '$a_{He}$', '$\\tau_{He}', '$log(\\xi)$']
# Input parameters for fake spectra
input_vals = np.array([0.08, 18000, 2, 0.1, 1.0, 1.0, 1.0, -4])

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
print ('log(xi)', log_xi_mcmc)
print ('\n Input parameter values:')
print (input_vals)

'''
fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, 12))
for i in range(ndim):
    axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
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