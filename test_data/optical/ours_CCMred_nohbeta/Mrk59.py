import corner
import emcee
import time
import os
import numpy as np
import matplotlib.pyplot as plt
import model_flux_ratio as mfr

from astropy.table import Table
from matplotlib.ticker import MaxNLocator

# Read in measured data (wavelength, flux ratios, and EWs)
flux_ratios = Table.read(os.getcwd()+'/test_data/Mrk59', format='ascii', delimiter=' ')

# Names of wavelengths of interest for MCMC
y_names = ['HeI+H83890', 'HeI4027', 'Hd', 'Hg', 'HeI4472', 'Hb', 'HeI5877', 'Ha', 'HeI6679', 'HeI7067']#, 'HeI10830']

# Balmer and Helium lines of interest for MCMC
balmer_lines = np.array([6564.612, 4862.721, 4341.684, 4102.891, 3890.166]) # Ha, Hb, Hg, Hd, H8
helium_lines = np.array([7067.198, 6679.994, 5877.299, 4472.755, 4027.328, 3890.151])

# Wavelengths we care about for MCMC
emis_lines = np.sort(np.concatenate((balmer_lines, helium_lines)))[1:] # [1:] to remove the duplicate ~3890 wavelength

# Measured data from spectra
EW_meas = np.array(flux_ratios['EW'])

y = np.array(flux_ratios['Flux Ratio']) # F(lambda) / F(H-beta)
y_error = np.array(flux_ratios['Flux Ratio Errors'])
x = np.zeros(y.size)

# Range of values for 8 parameters: y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, xi/n_HI
min_y_plus, max_y_plus = 0.05, 0.1  # fraction of singly ionized He; y+ = He+/H+
min_temp, max_temp = 5000, 25000  # electron temperature (K)
min_log_dens, max_log_dens = 0, 3  # log10(electron density) (cm^-3)
min_c_Hb, max_c_Hb = 0, 0.5  # reddening
min_a_H, max_a_H = 0, 10  # underlying stellar H absorption (Angstroms)
min_a_He, max_a_He = 0, 4  # underlying stellar HeI absorption (Angstroms)
min_tau_He, max_tau_He = 0, 5  # optical depth; range of values from Izotov & Thuan 2010
min_log_xi, max_log_xi = -6, -0.0969 # equals to xi=0-0.8; ratio of neutral hydrogen to singly ionized hydrogen densities; xi = n(HI)/n(HII)

# Set up MCMC
def get_model(theta):
	y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi = theta
	
	# Get density and neutral hydrogen out of log space
	dens = 10 ** log_dens
	xi = 10 ** log_xi

	# Take into account error on EW(Hb) by perturbing EW(Hb) by its measured EW error
	EW_Hb = np.random.normal(flux_ratios['EW'][np.where(flux_ratios['Wavelength'] == 4862.721)[0]][0], flux_ratios['EW Errors'][np.where(flux_ratios['Wavelength'] == 4862.721)[0]][0])

	# Continuum level ratio
	h = y * EW_Hb / EW_meas # relative to H-beta; i.e., h(lambda) / h(H-beta)
	
	# Model flux
	model_flux = np.zeros(10) # 10 emission line fluxes we want to model

	# Some values, calculated at Hbeta, for later use
	collisional_to_recomb_Hbeta = mfr.hydrogen_collision_to_recomb(xi, balmer_lines[1], temp)
	#AHbeta_AV = mfr.reddening_coefficient(balmer_lines[1])

	for w in range(len(emis_lines)):
		# Determine if working with hydrogen or helium line; within 3 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
		nearest_wave = emis_lines[np.where(np.abs(emis_lines - emis_lines[w]) < 3)[0]][0]
		# The above line is redundant for my input waves, but allows for cases where emis_lines[w] is some other array, say waves_of_interest[w], 
		# and not exactly at the wavelengths given in the emis_lines array (which is concatenated from arrays balmer_lines and helium_lines)
	
		# Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166)
		if nearest_wave in balmer_lines and nearest_wave != 3890.166:
			line_species = 'hydrogen'
			
			emissivity_ratio = mfr.hydrogen_emissivity_HS1987(emis_lines[w], temp, dens)
			a_H_at_wave = mfr.stellar_absorption(emis_lines[w], a_H, ion=line_species)
			collisional_to_recomb_ratio = mfr.hydrogen_collision_to_recomb(xi, emis_lines[w], temp)
			reddening_function = (mfr.reddening_coefficient(emis_lines[w]) ) - 1.
			
#			flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
#				( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
#				10**-(reddening_function * c_Hb)
			# Reparameterization of flux to use the continuum level
			flux = ( emissivity_ratio * ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
				10**-(reddening_function * c_Hb) * ( (EW_Hb + a_H)/(EW_Hb) ) ) - ( (a_H_at_wave / EW_Hb) * (h[w]) )
                  
		# Any HeI line besides the blended HeI+H8 line (HeI at 3890.151)
		elif nearest_wave in helium_lines and nearest_wave != 3890.151:
			line_species = 'helium'
            
			emissivity_ratio = mfr.helium_emissivity(emis_lines[w], temp, dens)
			a_He_at_wave = mfr.stellar_absorption(emis_lines[w], a_He, ion=line_species)
			optical_depth_at_wave = mfr.optical_depth_function(emis_lines[w], temp, dens, tau_He)
			collisional_to_recomb_ratio = 0. #mfr.helium_collision_to_recomb(emis_lines[w], temp, dens)
			reddening_function = (mfr.reddening_coefficient(emis_lines[w]) ) - 1.
			
			# Reparameterization of flux to use the continuum level
			flux = ( y_plus * emissivity_ratio *  optical_depth_at_wave * ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
				10**-(reddening_function * c_Hb) * ( (EW_Hb + a_H)/(EW_Hb) ) ) - ( (a_He_at_wave/EW_Hb) * (h[w]) )

		# The blended HeI+H8 line
		elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
			reddening_function = (mfr.reddening_coefficient(emis_lines[w]) ) - 1.
			
			# HeI 3890.151 contribution:
			line_species = 'helium'
			
			emissivity_ratio = mfr.helium_emissivity(emis_lines[w], temp, dens)
			a_He_at_wave = mfr.stellar_absorption(emis_lines[w], a_He, ion=line_species)            
			optical_depth_at_wave = mfr.optical_depth_function(emis_lines[w], temp, dens, tau_He)            
			collisional_to_recomb_ratio = 0. #mfr.helium_collision_to_recomb(emis_lines[w], temp, dens)            

			flux = ( y_plus * emissivity_ratio *  optical_depth_at_wave * ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                                10**-(reddening_function * c_Hb) * ( (EW_Hb + a_H)/(EW_Hb) ) ) - ( (a_He_at_wave/EW_Hb) * (h[w]) )

			# H8 contribution:
			line_species = 'hydrogen'

			emissivity_ratio = mfr.hydrogen_emissivity_HS1987(emis_lines[w], temp, dens)
			a_H_at_wave = mfr.stellar_absorption(emis_lines[w], a_H, ion=line_species)            
			collisional_to_recomb_factor = np.exp(( -13.6 * ((1/5**2)-(1/8**2)) ) / (8.6173303e-5 * temp)) # scale factor for going from C/R(Hg) to C/R(H8)
			collisional_to_recomb_ratio = collisional_to_recomb_factor * mfr.hydrogen_collision_to_recomb(xi, 4341.684, temp) # Calculate C/R(Hg) and multiply by above scale factor

			flux += ( emissivity_ratio * ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                                10**-(reddening_function * c_Hb) * ( (EW_Hb + a_H)/(EW_Hb) ) ) - ( (a_H_at_wave / EW_Hb) * (h[w]) )

		model_flux[w] = flux

	return model_flux

# Define the probability function as likelihood * prior.
def lnprior(theta):
	y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi = theta

	if min_y_plus <= y_plus <= max_y_plus and \
		min_temp <= temp <= max_temp and \
		min_log_dens <= log_dens <= max_log_dens and \
		min_c_Hb <= c_Hb <= max_c_Hb and \
		min_a_H <= a_H <= max_a_H and \
		min_a_He <= a_He <= max_a_He and \
		min_tau_He <= tau_He <= max_tau_He and \
		min_log_xi <= log_xi <= max_log_xi:
		return 0.0 #( -(temp - 18000)**2 / (0.2*18000)**2 ) # Test

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
np.save('{0:s}_{1:d}walkers_{2:d}steps'.format('Mrk59', nwalkers, nmbr), sampler.chain)

print('Making plots')
burnin = int(0.8 * nmbr)

samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
# Names of 8 parameters and input 'true' parameter values
prenams = ['y+', 'temperature', '$log(n_{e})$', 'c(H\\beta)', '$a_{H}$', '$a_{He}$', '$\\tau_{He}', '$log(\\xi)$'] # '$n_{HI}$']

print ('Best parameter values:')
y_plus_mcmc, temp_mcmc, log_dens_mcmc, c_Hb_mcmc, a_H_mcmc, a_He_mcmc, tau_He_mcmc, log_xi_mcmc = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))

print ('y+', y_plus_mcmc)
print ('T', temp_mcmc)
print ('log(n_e)', log_dens_mcmc)
print ('c(Hb)', c_Hb_mcmc)
print ('a_H', a_H_mcmc)
print ('a_He', a_He_mcmc)
print ('tau_He', tau_He_mcmc)
print ('log(xi)', log_xi_mcmc)
