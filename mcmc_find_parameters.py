import corner
import emcee
import time
import numpy as np
import matplotlib.pyplot as plt
import model_flux_ratio as mfr

from astropy.table import Table
from matplotlib.ticker import MaxNLocator

# Read in generated/fake flux ratios
flux_ratios = Table.read('/Users/thsyu/Dropbox/BCDs/primordial_helium/test_output_flux', format='ascii', delimiter=' ')

# Names of wavelenghts of interest
# y_names = ['H12', 'H11', 'H10', 'H9', 'H8+NeIII', 'HeI3890', 'Heps+HeI3970', 'Hd', 'Hg', 'HeI4472', \
#      'Hb', 'HeI4923', 'HeI5027', 'HeI5877', 'Ha', 'HeI6679', 'HeI7067', 'HeI7283']
y_names = ['HeI3890', 'HeI4027', 'Hd', 'Hg', 'HeI4472', 'Hb', 'HeI4923', 'HeI5027', 'HeI5877', 'Ha', 'HeI6679',
           'HeI7067', 'HeI7283']

# Balmer and Helium lines
# balmer_lines = np.array([6564.612, 4862.721, 4341.684, 4102.891, 3971.195, 3869.81, 3836.472, 3798.976, 3771.701, 3751.217])
balmer_lines = np.array([6564.612, 4862.721, 4341.684, 4102.891])
helium_lines = np.array([7283.356, 7067.198, 6679.994, 5877.299, 5017.079, 4923.304, 4472.755, 4027.328, 3890.151])
emis_lines = np.sort(np.concatenate((balmer_lines, helium_lines)))

# Measured data from spectra
EWs = np.array(flux_ratios['EW'])
EW_Hb = flux_ratios['EW'][np.where(flux_ratios['Wavelength'] == 4862.721)[0]]
y = np.array(flux_ratios['Flux Ratio'])
y_error = np.array(flux_ratios['Flux Ratio'] * 0.02)
x = np.zeros(y.size)

# Range of values for 8 parameters: y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, n_HI
min_y_plus, max_y_plus = 0.05, 0.1  # fraction of singly ionized He; y+ = He+/H+
min_temp, max_temp = 5000, 25000  # electron temperature (K)
min_dens, max_dens = 0, 14  # log10(electron density) (cm^-3)
min_c_Hb, max_c_Hb = 0, 0.5  # reddening
min_a_He, max_a_He = 0, 4  # underlying stellar HeI absorption (Angstroms)
min_a_H, max_a_H = 0, 10  # underlying stellar H absorption (Angstroms)
min_tau_He, max_tau_He = 0, 5  # optical depth; range of values from Izotov & Thuan 2010
min_n_HI, max_n_HI = 1e-4, 1e-1  # neutral hydrogen density (cm^-3)


# Set up MCMC
def get_model(theta):
    y_plus, temp, log_dens, c_Hb, a_He, a_H, tau_He, n_HI = theta

    model_flux = []
    dens = 10 ** log_dens
    eta = n_HI / dens

    for w in range(len(emis_lines)):
        # Determine if working with hydrogen or helium line
        nearest_wave = emis_lines[np.where(np.abs(emis_lines - emis_lines[w]) < 3)[0]][0]

        if nearest_wave in balmer_lines:
            line_species = 'hydrogen'

            emissivity_ratio = mfr.hydrogen_emissivity(emis_lines[w], temp, dens)
            a_H_at_wave = mfr.stellar_absorption(emis_lines[w], a_H, ion=line_species)
            collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(eta, emis_lines[w], temp)
            collisional_to_recomb_Hbeta = mfr.hydrogen_collision_to_recomb(eta, balmer_lines[1], temp)
            reddening_function = (mfr.f_lambda_avg_interp(emis_lines[w]) / mfr.f_lambda_avg_interp(4862.721)) - 1.

            flux = emissivity_ratio * (((EW_Hb + a_H) / (EW_Hb)) / ((EWs[w] + a_H_at_wave) / (EWs[w]))) * \
                   ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) * \
                   10 ** -(reddening_function * c_Hb)

        elif nearest_wave in helium_lines:
            line_species = 'helium'

            emissivity_ratio = mfr.helium_emissivity(emis_lines[w], temp, dens)
            a_He_at_wave = mfr.stellar_absorption(emis_lines[w], a_He, ion=line_species)
            optical_depth_at_wave = mfr.optical_depth_function(emis_lines[w], temp, dens, tau_He)
            collisional_to_recomb_ratio = mfr.helium_collision_to_recomb(emis_lines[w], temp, dens)
            collisional_to_recomb_Hbeta = mfr.hydrogen_collision_to_recomb(eta, balmer_lines[1], temp)
            reddening_function = (mfr.f_lambda_avg_interp(emis_lines[w]) / mfr.f_lambda_avg_interp(4862.721)) - 1.
            flux = y_plus * emissivity_ratio * (((EW_Hb + a_H) / (EW_Hb)) / ((EWs[w] + a_He_at_wave) / (EWs[w]))) * \
                   optical_depth_at_wave * ((1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta)) * \
                   10 ** -(reddening_function * c_Hb)

        model_flux.append(flux)

    return model_flux


# Define the probability function as likelihood * prior.
def lnprior(theta):
    y_plus, temp, log_dens, c_Hb, a_He, a_H, tau_He, n_HI = theta

    if min_y_plus <= y_plus <= max_y_plus and \
            min_temp <= temp <= max_temp and \
            min_dens <= dens <= max_dens and \
            min_c_Hb <= c_Hb <= max_c_Hb and \
            min_a_He <= a_He <= max_a_He and \
            min_a_H <= a_H <= max_a_H and \
            min_tau_He <= tau_He <= max_tau_He and \
            min_n_HI <= n_HI <= max_n_HI:
        return 0.0

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
ndim, nwalkers = 8, 100

pos = [np.array([np.random.uniform(min_y_plus, max_y_plus),
                 np.random.uniform(min_temp, max_temp),
                 np.random.uniform(min_dens, max_dens),
                 np.random.uniform(min_c_Hb, max_c_Hb),
                 np.random.uniform(min_a_He, max_a_He),
                 np.random.uniform(min_a_H, max_a_H),
                 np.random.uniform(min_tau_He, max_tau_He),
                 np.random.uniform(min_n_HI, max_n_HI)]) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, y_error), threads=ndim)

print('Running MCMC...')
nmbr = 100
a = time.time()
print('Started at', a)
for i, result in enumerate(sampler.run_mcmc(pos, nmbr, rstate0=np.random.get_state())):
    if True:  # (i+1) % 100 == 0:
        print("{0:5.1%}".format(float(i) / nmbr))
print('Done!')
print((time.time() - a) / 60.0, 'mins')

print('Saving samples')
np.save("{0:s}_samples{1:d}.npy".format('test_MCMC', nmbr), sampler.chain)

print('Making plots')
burnin = 0.1 * nmbr

samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
# Names of 8 parameters and input 'true' parameter values
prenams = ['y+', 'temperature', '$n_{e}$', 'c(H\\beta)', '$a_{He}$', '$a_{H}$', '$\\tau_{He}', '$n_{HI}$']
input_vals = np.array([0.08, 18000, 100, 0.1, 1.0, 1.0, 1.0, 1e-2])

print ('Best parameter values:')
y_plus_mcmc, temp_mcmc, dens_mcmc, c_Hb_mcmc, a_H_mcmc, a_He_mcmc, tau_He_mcmc, n_HI_mcmc = map(
    lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
print ('y+', y_plus_mcmc)
print ('T', temp_mcmc)
print ('n_e', dens_mcmc)
print ('c(Hb)', c_Hb_mcmc)
print ('a_H', a_H_mcmc)
print ('a_He', a_He_mcmc)
print ('tau_He', tau_He_mcmc)
print ('n_HI', n_HI_mcmc)
print (input_vals)

fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, 12))
for i in range(ndim):
    axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
    axes[i].yaxis.set_major_locator(MaxNLocator(5))
    axes[i].axvline(10, color='red')
    axes[i].set_ylabel(prenams[i])
axes[7].set_xlabel('Steps')
fig.tight_layout(h_pad=0.0)
fig.savefig("{0:s}_samples{1:d}.pdf".format('test_MCMC_time_evol', nmbr), overwrite=True)

fig = corner.corner(samples, labels=prenams, truths=input_vals)
fig.savefig('{0:s}_samples{1:d}.pdf'.format('test_MCMC_params', nmbr), overwrite=True)