#########
# yMCMC #
#########
# Code to model the flux ratio of an emission line, given 8 parameters:
## y+, helium abundance
## T, temperature
## n_e, density
## c(Hb), reddening
## a_H, hydrogen stellar absorption
## a_He, helium stellar absorption
## tau_He, optical depth
## xi, ratio of neutral to ionized hydrogen densities
# Flux Ratio = Emissivity Ratio * EW+absorption Ratio * Optical Depth * Collisional to Recombination Ratio * Reddening

###########
# Updates #
###########
# 2019-06-04: adding hydrogen emissivity to return either an emissivity as a ratio to Hb (default) or just emissivity
#             changed C/R using R+2015 and A+2002&O1983 to use emissivities as a substitute for recombination rate (removing need for HS1987 recombination rates)
# 2019-03-28: fixed bug in underlying stellar absorption (grabbing wrong value for H8 and HeI5017)
#             edits to generate_emission_line_ratio
#              - combined generate_nir_emission_line_ratio and generate_optical_emission_line_ratio
#              - added option for using S2018 emissivities (default) or HS1987 emissivities
#             cleaned up formatting
# 2019-03-07: removed Porter's HeI emissivities on older coarse grid with larger parameter range T=5000-25000, log(n_e)=0-14
#               - changed temperature priors to be 10000-22000 to match fine mesh
#               - hoping to speed up code by eliminating reading in these 'old' emissivities..!
# 2019-03-05: added Porter's HeI emissivities on a finer parametric mesh at T=10000-22000, log(n_e)=0-4
# 2019-03-01: added options for methods of calculating Hbeta emissivity to use as ratio for a HeI emissivity in helium_emissivity_PFSD2012
#               - defaults to S2018 Hbeta emissivity (as of 03-28)
#               - able to reproduce input parameters to within 1-sigma!
# 2019-02-25: edited model flux generator to produce F(HeI10830)/F(Pg) instead F(Hbeta) to mirror our measured data
#               - separated this into generate_nir_emission_line_ratio, which uses HS1987 emissivities for all but Pg, which uses S2018
#               - able to reproduce input parameters to within 1-sigma!
# 2019-02-12: added F(HeI10830)/F(Hbeta) to the model flux
#               - note: still using HS1987 emissivities to match AOS2015, which means needing (1+C/R(Hbeta))
#               - returns F(HeI10830)/F(Hbeta), not F(HeI10830)/F(Pg)
#               - able to reproduce input parameters to within 1-sigma!
# 2019-02-11: created to begin NIR development!

#############
# Test Flux #
#############
# mfr.generate_emission_line_ratio('test_output_flux', [3890.166, 4027.328, 4102.891, 4341.684, 4472.755, 4862.721, 5017.079, 5877.299, 6564.612, 6679.994, 7067.198, 10833.306], [10, 10, 75, 100, 10, 250, 5, 10, 350, 10, 5, 200], 250, 0.08, 18000, 2, 0.1, 1.0, 1.0, 1.0, -4, EW_Pg=50.)

# Imports
import os
import pdb
import numpy as np
import scipy.interpolate as interp
from functools import reduce
from astropy.table import Table


##########
# Tables #
##########
# Load in tables we'll need
path = os.getcwd()

hydrogen_emis = Table.read(path+'/tables/hydrogen_emissivity_S2018', format='ascii', delimiter='\t')
hydrogencoeff = Table.read(path+'/tables/hydrogen_emissivity_HS1987', format='ascii', delimiter='\t')
helium_emis = Table.read(path+'/tables/helium_emissivity', format='ascii', delimiter='\t')
helium_finemesh_emis = Table.read(path+'/tables/helium_emissivity_finemesh', format='ascii', delimiter='\t')
heliumcoeff = Table.read(path+'/tables/helium_emissivity_coeff', format='ascii', delimiter='\t')
helium_optical_depth = Table.read(path+'/tables/helium_optical_depth', format='ascii', delimiter='\t')
hydrogen_CR_coeff = Table.read(path+'/tables/hydrogen_CR_coeff', format='ascii', delimiter='\t')
helium_CR_coeff = Table.read(path+'/tables/helium_CR_coeff', format='ascii', delimiter='\t')

# Vacuum wavelengths of Balmer lines Ha, Hb, Hg, Hd, H8 for MCMC
hydrogen_lines = np.array([10941.082, 6564.612, 4862.721, 4341.684, 4102.891, 3890.166, 18756.096, 12821.578, 40522.79]) #last 3 lines are Pa, Pb, Br-a (4-3, 5-3, 5-4)

# Vacuum wavelengths of Helium lines for MCMC
helium_lines = np.array([10833.306, 7067.198, 6679.994, 5877.299, 5017.079, 4472.755, 4027.328, 3890.151])

##############
# Emissivity #
##############

# Hydrogen
# --------
# Interpolated Storey 2018 hydrogen emissivities
ha_RBS = np.zeros((21,6)) # H-alpha, 3-->2
hb_RBS = np.zeros((21,6)) # H-beta, 4-->2
hg_RBS = np.zeros((21,6)) # H-gamma, 5-->2
hd_RBS = np.zeros((21,6)) # H-delta, 6-->2
h8_RBS = np.zeros((21,6)) # H8, 8-->2

pa_RBS = np.zeros((21,6)) # Paschen-alpha, 4-->3
pb_RBS = np.zeros((21,6)) # Paschen-beta, 5-->3
pg_RBS = np.zeros((21,6)) # Paschen-gamma, 6-->3

bra_RBS = np.zeros((21,6)) # Brackett-alpha, 5-->4

# Find corresponding data in hydrogen emissivity table for each emission line
for t in range(len(np.arange(5000, 26000, 1000))):
    ha_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 3)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 2)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    hb_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 4)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 2)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    hg_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 5)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 2)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    hd_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 6)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 2)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    h8_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 8)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 2)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    pa_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 4)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 3)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    pb_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 5)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 3)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    pg_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 6)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 3)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    bra_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 5)[0], \
                                                          np.where(hydrogen_emis['Nl'] == 4)[0], \
                                                          np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]

# Linear interpolation on the hydrogen emissivities
S2018_ha_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), ha_RBS, kx=1, ky=1)
S2018_hb_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), hb_RBS, kx=1, ky=1)
S2018_hg_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), hg_RBS, kx=1, ky=1)
S2018_hd_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), hd_RBS, kx=1, ky=1)
S2018_h8_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), h8_RBS, kx=1, ky=1)
S2018_pa_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), pa_RBS, kx=1, ky=1)
S2018_pb_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), pb_RBS, kx=1, ky=1)
S2018_pg_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), pg_RBS, kx=1, ky=1)
S2018_bra_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), bra_RBS, kx=1, ky=1)
# Cubic interpolation on the hydrogen emissivities
S2018_ha_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), ha_RBS, kx=3, ky=3)
S2018_hb_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), hb_RBS, kx=3, ky=3)
S2018_hg_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), hg_RBS, kx=3, ky=3)
S2018_hd_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), hd_RBS, kx=3, ky=3)
S2018_h8_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), h8_RBS, kx=3, ky=3)
S2018_pa_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), pa_RBS, kx=3, ky=3)
S2018_pb_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), pb_RBS, kx=3, ky=3)
S2018_pg_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), pg_RBS, kx=3, ky=3)
S2018_bra_cubic = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0,6), bra_RBS, kx=3, ky=3)

def hydrogen_emissivity_S2018(wave, temp, dens, deg='linear', ratio=True):
    '''
    Calculate the emissivity of a hydrogen line,
    defaulted to be relative to H-beta, using
    P. Storey's 2018 hydrogen emissivities.
    Option to return just a hydrogen line
    emissivity.

    These are interpolated using a
    RectBivariateSpline() -- linear by default,
    option of cubic

    Parameters
    ----------
    wave : float
        Wavelength of the hydrogen line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (in cm^-3)
    deg : str
        Degree of RBS interpolation; default is linear
    ratio : True/False
        Return the emissivity as a ratio to H-beta?

    Returns
    -------
    emissivity : float
        The E(lambda)/E(Hbeta) ratio (default);
        optional return of just E(lambda)
    '''
    # Reformat the density
    logdens = np.log10(dens)

    # Match to Balmer line of interest
    idx = np.where(np.abs(wave - hydrogen_lines) < 3.5)[0][0]

    if deg == 'linear':
        # H-beta emissivity, for calculating the ratio of emissivities
        Hbeta_emis = S2018_hb_lin(temp, logdens).flatten()
        # Hydrogen emissivity
        if idx == 0: # P-gamma
            Xt = S2018_pg_lin(temp, logdens).flatten()
        elif idx == 1: # H-alpha:
            Xt = S2018_ha_lin(temp, logdens).flatten()
        elif idx == 2: # H-beta:
            Xt = Hbeta_emis
        elif idx == 3: # H-gamma
            Xt = S2018_hg_lin(temp, logdens).flatten()
        elif idx == 4: # H-delta
            Xt = S2018_hd_lin(temp, logdens).flatten()
        elif idx == 5: # H8
            Xt = S2018_h8_lin(temp, logdens).flatten()
        elif idx == 6: # Pa
            Xt = S2018_pa_lin(temp, logdens).flatten()
        elif idx == 7: # Pb
            Xt = S2018_pb_lin(temp, logdens).flatten()
        elif idx == 8: # Br-a
            Xt = S2018_bra_lin(temp, logdens).flatten()
        else:
            print('Not ready for this hydrogen line!')
            pdb.set_trace()
    elif deg == 'cubic':
        # H-beta emissivity, for calculating the ratio of emissivities
        Hbeta_emis = S2018_hb_cubic(temp, logdens).flatten()
        # Hydrogen emissivity
        if idx == 0: # P-gamma
            Xt = S2018_pg_cubic(temp, logdens).flatten()
        elif idx == 1: # H-alpha:
            Xt = S2018_ha_cubic(temp, logdens).flatten()
        elif idx == 2: # H-beta:
            Xt = Hbeta_emis
        elif idx == 3: # H-gamma
            Xt = S2018_hg_cubic(temp, logdens).flatten()
        elif idx == 4: # H-delta
            Xt = S2018_hd_cubic(temp, logdens).flatten()
        elif idx == 5: # H8
            Xt = S2018_h8_cubic(temp, logdens).flatten()
        elif idx == 6: # Pa
            Xt = S2018_pa_cubic(temp, logdens).flatten()
        elif idx == 7: # Pb
            Xt = S2018_pb_cubic(temp, logdens).flatten()
        elif idx == 8: # Br-a
            Xt = S2018_bra_cubic(temp, logdens).flatten()
        else:
            print('Not ready for this hydrogen line!')
            pdb.set_trace()
    else:
        print ('Not ready for this degree of interpolation!')
        pdb.set_trace()

    if ratio is True:
        Xt = Xt[0] / Hbeta_emis[0]

    return Xt


def hydrogen_emissivity_HS1987(wave, temp, dens):
    '''
    Calculate the emissivity of a Balmer line
    relative to H(beta) using Hummer and Storey's
    1987 hydrogen emissivities, reparameterized
    by Aver et al. in 2010

    These do not include the collisional to
    recombination correction

    Parameters
    ----------
    wave : float
        Wavelength of the hydrogen line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)

    Returns
    -------
    emissivity : float
        The E(lambda)/E(H(beta)) ratio
    '''
    # Redefine the temperature
    T4 = temp / 10000.

    # Match Balmer line of interest to relevant rows in Table 3 of AOS 2010
    idx = np.where(np.abs(wave - hydrogen_lines) < 3)[0][0]

    if idx == 1:
        line = str('Ha')
    elif idx == 2:
        line = str('Hb')
    elif idx == 3:
        line = str('Hg')
    elif idx == 4:
        line = str('Hd')
    elif idx == 5:
        line = str('H8')
    else:
        print('Not ready for this hydrogen line')
        pdb.set_trace()

    if line == 'Hb':
        Xt = 1.
    # From AOS 2010, a little above Eq. 4.1; fitted to data from Hummer & Storey
    elif line == 'H8':
        Xt = 0.104 * (T4**0.046)
    else:
        cij = np.array(hydrogencoeff[line]).reshape((3, 3))
        Xt = 0.

        for i in range(cij.shape[0]):
            for j in range(cij.shape[1]):
                # Balmer emissivity; from Equation in Section 3.1 Hydrogen emission of AOS 2010 Citation (3)
                Xt += cij[i][j] * (np.log10(T4) ** i) * (np.log10(dens) ** j)

    return Xt


# Helium
# ------
# Interpolated PFSD 2012/2013 erratum HeI emissivities; T=5000-25000, log(n_e)=1-14
#### These are *not* interpolated on a log(n_e) scale!
#### Currently not being used -- using fine mesh HeI emissivities below instead
# Linear interpolation on the helium emissivities
#HeI3889_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['3889A'].reshape((14, 21)), kx=1, ky=1)
#HeI4026_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4026A'].reshape((14, 21)), kx=1, ky=1)
#HeI4471_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4471A'].reshape((14, 21)), kx=1, ky=1)
#HeI5016_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5016A'].reshape((14, 21)), kx=1, ky=1)
#HeI5876_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5876A'].reshape((14, 21)), kx=1, ky=1)
#HeI6678_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['6678A'].reshape((14, 21)), kx=1, ky=1)
#HeI7065_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['7065A'].reshape((14, 21)), kx=1, ky=1)
#HeI10833_lin = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['10830A'].reshape((14, 21)), kx=1, ky=1)
# Cubic interpolation on the helium emissivities
#HeI3889_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['3889A'].reshape((14, 21)), kx=3, ky=3)
#HeI4026_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4026A'].reshape((14, 21)), kx=3, ky=3)
#HeI4471_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4471A'].reshape((14, 21)), kx=3, ky=3)
#HeI5016_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5016A'].reshape((14, 21)), kx=3, ky=3)
#HeI5876_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5876A'].reshape((14, 21)), kx=3, ky=3)
#HeI6678_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['6678A'].reshape((14, 21)), kx=3, ky=3)
#HeI7065_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['7065A'].reshape((14, 21)), kx=3, ky=3)
#HeI10833_cubic = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['10830A'].reshape((14, 21)), kx=3, ky=3)
# ------
# Interpolated Porter HeI emissivities on finer mesh, from AOPS 2013; T=10000-22000, log(n_e)=0-4
dens_finemesh = np.array(np.unique(helium_finemesh_emis['log n_e']))
temp_finemesh = np.arange(10000, 22250, step=250)
# Linear interpolation on the helium emissivities
HeI3889_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['3889A'].reshape((31, 49)), kx=1, ky=1)
HeI4026_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['4026A'].reshape((31, 49)), kx=1, ky=1)
HeI4471_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['4471A'].reshape((31, 49)), kx=1, ky=1)
HeI5016_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['5016A'].reshape((31, 49)), kx=1, ky=1)
HeI5876_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['5876A'].reshape((31, 49)), kx=1, ky=1)
HeI6678_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['6678A'].reshape((31, 49)), kx=1, ky=1)
HeI7065_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['7065A'].reshape((31, 49)), kx=1, ky=1)
HeI10833_finemesh_lin = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['10830A'].reshape((31, 49)), kx=1, ky=1)
# Cubic interpolation on the helium emissivities
HeI3889_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['3889A'].reshape((31, 49)), kx=3, ky=3)
HeI4026_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['4026A'].reshape((31, 49)), kx=3, ky=3)
HeI4471_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['4471A'].reshape((31, 49)), kx=3, ky=3)
HeI5016_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['5016A'].reshape((31, 49)), kx=3, ky=3)
HeI5876_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['5876A'].reshape((31, 49)), kx=3, ky=3)
HeI6678_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['6678A'].reshape((31, 49)), kx=3, ky=3)
HeI7065_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['7065A'].reshape((31, 49)), kx=3, ky=3)
HeI10833_finemesh_cubic = interp.RectBivariateSpline(dens_finemesh, temp_finemesh, helium_finemesh_emis['10830A'].reshape((31, 49)), kx=3, ky=3)

def helium_emissivity_PFSD2012(wave, temp, dens, deg='linear', ratio='porter'):
    '''
    Calculate the emissivity of a HeI line
    using Porter's 2013 fine mesh emissivities

    These include the collisional to recombination
    correction and are interpolated using a 
    RectBivariateSpline() -- linear by default,
    option of cubic

    Parameters
    ----------
    wavelength : float
        Wavelength of the HeI line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3), not in log!
    deg : str
        Degree of RBS interpolation; default is linear
    ratio : str
        Method of determining Hbeta emissivity for
        emissivitity ratio. Defaults to Porter's 
        reparamerization given in Eq. 3.1 of AOS2010.
        Option for Storey's 2018 hydrogen emissivities,
        taking on same degree of interpolation as HeI.

    Returns
    -------
    emissivity : float
        The E(HeI)/E(H-beta) ratio (default).
        Units of emissivity are in ergs*cm^3*s^-1
    '''
    # Find which helium line we are working on
    HeI_lines = np.array([3889, 4026, 4471, 5016, 5876, 6678, 7065, 10833])
    try:
        HeI_line = str(HeI_lines[np.where(np.abs(HeI_lines - wave) < 3.5)[0]][0])
    except:
        print ('Not ready for this helium line!')
        pdb.set_trace()

    logdens = np.log10(dens)

    # Temperature and density fall within finer parametric mesh
    if deg == 'linear':
        if HeI_line == '3889':
            HeI_emis = 10 ** HeI3889_finemesh_lin(logdens, temp)[0][0]
        elif HeI_line == '4026':
            HeI_emis = 10 ** HeI4026_finemesh_lin(logdens, temp)[0][0]
        elif HeI_line == '4471':
            HeI_emis = 10 ** HeI4471_finemesh_lin(logdens, temp)[0][0]
        elif HeI_line == '5016':
            HeI_emis = 10 ** HeI5016_finemesh_lin(logdens, temp)[0][0]
        elif HeI_line == '5876':
            HeI_emis = 10 ** HeI5876_finemesh_lin(logdens, temp)[0][0]
        elif HeI_line == '6678':
            HeI_emis = 10 ** HeI6678_finemesh_lin(logdens, temp)[0][0]
        elif HeI_line == '7065':
            HeI_emis = 10 ** HeI7065_finemesh_lin(logdens, temp)[0][0]
        elif HeI_line == '10833':
            HeI_emis = 10 ** HeI10833_finemesh_lin(logdens, temp)[0][0]
    elif deg == 'cubic':
        if HeI_line == '3889':
            HeI_emis = 10 ** HeI3889_finemesh_cubic(logdens, temp)[0][0]
        elif HeI_line == '4026':
            HeI_emis = 10 ** HeI4026_finemesh_cubic(logdens, temp)[0][0]
        elif HeI_line == '4471':
            HeI_emis = 10 ** HeI4471_finemesh_cubic(logdens, temp)[0][0]
        elif HeI_line == '5016':
            HeI_emis = 10 ** HeI5016_finemesh_cubic(logdens, temp)[0][0]
        elif HeI_line == '5876':
            HeI_emis = 10 ** HeI5876_finemesh_cubic(logdens, temp)[0][0]
        elif HeI_line == '6678':
            HeI_emis = 10 ** HeI6678_finemesh_cubic(logdens, temp)[0][0]
        elif HeI_line == '7065':
            HeI_emis = 10 ** HeI7065_finemesh_cubic(logdens, temp)[0][0]
        elif HeI_line == '10833':
            HeI_emis = 10 ** HeI10833_finemesh_cubic(logdens, temp)[0][0]
    else:
        print ('Not ready for this degree of interpolation!')
        pdb.set_trace()

    # Calculate Hbeta emissivity for a ratio
    if ratio == 'porter':
        # H beta emissivity; from Equation 3.1 of Citation (3) AOS 2010
        Hbeta_emis = (-2.6584e5 - (1420.9 * (np.log(temp) ** 2.)) + (35546 * np.log(temp)) + (6.5669e5 / np.log(temp))) \
                     * (1 / temp) * 1e-25
    elif ratio == 'storey':
        if deg == 'linear':
            Hbeta_emis = S2018_hb_lin(temp, np.log10(dens))[0][0]
        elif deg == 'cubic':
            Hbeta_emis = S2018_hb_cubic(temp, np.log10(dens))[0][0]
    else:
        print ('Not prepared for this method of H-beta emissivity')

    return HeI_emis / Hbeta_emis

def helium_emissivity_coarse_PFSD2012(wave, temp, dens, deg='linear', ratio='porter'):
    '''
    Calculate the emissivity of a HeI line
    using Porter et al.'s 2013 emissivities
    (which are edits to emissivities in Porter,
    Ferland, Storey, & Detisch 2012 MNRAS 425, L28)

    These include the collisional to recombination
    correction and are interpolated using a 
    RectBivariateSpline() -- linear by default,
    option of cubic

    Parameters
    ----------
    wavelength : float
        Wavelength of the HeI line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3), not in log!
    deg : str
        Degree of RBS interpolation; default is linear
    ratio : str
        Method of determining Hbeta emissivity for
        emissivitity ratio. Defaults to Porter's 
        reparamerization given in Eq. 3.1 of AOS2010.
        Option for Storey's 2018 hydrogen emissivities,
        taking on same degree of interpolation as HeI.

    Returns
    -------
    emissivity : float
        The E(HeI)/E(H-beta) ratio (default).
        Units of emissivity are in ergs*cm^3*s^-1
    '''
    # Find which helium line we are working on
    HeI_lines = np.array([3889, 4026, 4471, 5016, 5876, 6678, 7065, 10833])
    try:
        HeI_line = str(HeI_lines[np.where(np.abs(HeI_lines - wave) < 3.5)[0]][0])
    except:
        print ('Not ready for this helium line!')
        pdb.set_trace()

    logdens = np.log10(dens)

    # Temperature or density fall outside the fine mesh grid; interpolate amongst online PFSD2012 emissivities
    if temp < 10000 or temp > 22000 or logdens > 4:
        if deg == 'linear':
            if HeI_line == '3889':
                HeI_emis = 10 ** HeI3889_lin(dens, temp)[0][0]
            elif HeI_line == '4026':
                HeI_emis = 10 ** HeI4026_lin(dens, temp)[0][0]
            elif HeI_line == '4471':
                HeI_emis = 10 ** HeI4471_lin(dens, temp)[0][0]
            elif HeI_line == '5016':
                HeI_emis = 10 ** HeI5016_lin(dens, temp)[0][0]
            elif HeI_line == '5876':
                HeI_emis = 10 ** HeI5876_lin(dens, temp)[0][0]
            elif HeI_line == '6678':
                HeI_emis = 10 ** HeI6678_lin(dens, temp)[0][0]
            elif HeI_line == '7065':
                HeI_emis = 10 ** HeI7065_lin(dens, temp)[0][0]
            elif HeI_line == '10833':
                HeI_emis = 10 ** HeI10833_lin(dens, temp)[0][0]
        elif deg == 'cubic':
            if HeI_line == '3889':
                HeI_emis = 10 ** HeI3889_cubic(dens, temp)[0][0]
            elif HeI_line == '4026':
                HeI_emis = 10 ** HeI4026_cubic(dens, temp)[0][0]
            elif HeI_line == '4471':
                HeI_emis = 10 ** HeI4471_cubic(dens, temp)[0][0]
            elif HeI_line == '5016':
                HeI_emis = 10 ** HeI5016_cubic(dens, temp)[0][0]
            elif HeI_line == '5876':
                HeI_emis = 10 ** HeI5876_cubic(dens, temp)[0][0]
            elif HeI_line == '6678':
                HeI_emis = 10 ** HeI6678_cubic(dens, temp)[0][0]
            elif HeI_line == '7065':
                HeI_emis = 10 ** HeI7065_cubic(dens, temp)[0][0]
            elif HeI_line == '10833':
                HeI_emis = 10 ** HeI10833_cubic(dens, temp)[0][0]
        else:
            print ('Not ready for this degree of interpolation!')
            pdb.set_trace()
    # Temperature and density fall within finer parametric mesh
    else:
        print ('Using emissivities on fine mesh')
        if deg == 'linear':
            if HeI_line == '3889':
                HeI_emis = 10 ** HeI3889_finemesh_lin(logdens, temp)[0][0]
            elif HeI_line == '4026':
                HeI_emis = 10 ** HeI4026_finemesh_lin(logdens, temp)[0][0]
            elif HeI_line == '4471':
                HeI_emis = 10 ** HeI4471_finemesh_lin(logdens, temp)[0][0]
            elif HeI_line == '5016':
                HeI_emis = 10 ** HeI5016_finemesh_lin(logdens, temp)[0][0]
            elif HeI_line == '5876':
                HeI_emis = 10 ** HeI5876_finemesh_lin(logdens, temp)[0][0]
            elif HeI_line == '6678':
                HeI_emis = 10 ** HeI6678_finemesh_lin(logdens, temp)[0][0]
            elif HeI_line == '7065':
                HeI_emis = 10 ** HeI7065_finemesh_lin(logdens, temp)[0][0]
            elif HeI_line == '10833':
                HeI_emis = 10 ** HeI10833_finemesh_lin(logdens, temp)[0][0]
        elif deg == 'cubic':
            if HeI_line == '3889':
                HeI_emis = 10 ** HeI3889_finemesh_cubic(logdens, temp)[0][0]
            elif HeI_line == '4026':
                HeI_emis = 10 ** HeI4026_finemesh_cubic(logdens, temp)[0][0]
            elif HeI_line == '4471':
                HeI_emis = 10 ** HeI4471_finemesh_cubic(logdens, temp)[0][0]
            elif HeI_line == '5016':
                HeI_emis = 10 ** HeI5016_finemesh_cubic(logdens, temp)[0][0]
            elif HeI_line == '5876':
                HeI_emis = 10 ** HeI5876_finemesh_cubic(logdens, temp)[0][0]
            elif HeI_line == '6678':
                HeI_emis = 10 ** HeI6678_finemesh_cubic(logdens, temp)[0][0]
            elif HeI_line == '7065':
                HeI_emis = 10 ** HeI7065_finemesh_cubic(logdens, temp)[0][0]
            elif HeI_line == '10833':
                HeI_emis = 10 ** HeI10833_finemesh_cubic(logdens, temp)[0][0]
        else:
            print ('Not ready for this degree of interpolation!')
            pdb.set_trace()

    # Calculate Hbeta emissivity for a ratio
    if ratio == 'porter':
        # H beta emissivity; from Equation 3.1 of Citation (3) AOS 2010
        Hbeta_emis = (-2.6584e5 - (1420.9 * (np.log(temp) ** 2.)) + (35546 * np.log(temp)) + (6.5669e5 / np.log(temp))) \
                     * (1 / temp) * 1e-25
    elif ratio == 'storey':
        if deg == 'linear':
            Hbeta_emis = S2018_hb_lin(temp, np.log10(dens))[0][0]
        elif deg == 'cubic':
            Hbeta_emis = S2018_hb_cubic(temp, np.log10(dens))[0][0]
    else:
        print ('Not prepared for this method of H-beta emissivity')

    return HeI_emis / Hbeta_emis


def helium_emissivity_PFM2007(wave, temp):
    '''
    Calculate the emissivity of a HeI line
    using Porter et al.'s 2007 work

    Defaults to returning the emissivity as
    a ratio relative to H-beta and does not
    include the collisional to recombination
    correction
    
    Parameters
    ----------
    wavelength : float
        Wavelength of the HeI line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
        
    Returns
    -------
    emissivity : float
        The E(HeI)/E(H-beta) ratio (default).
        Units of emissivity are in ergs*cm^3*s^-1
    '''
    # Find the row in Porter's 2007 emissivities corresponding to HeI wavelength of interest
    idx = np.abs(heliumcoeff['Wavelength'] - wave).argmin()

    a = heliumcoeff['a'][idx]
    b = heliumcoeff['b'][idx]
    c = heliumcoeff['c'][idx]
    d = heliumcoeff['d'][idx]
    
    # HeI emissivity; from Equation A1 of Citation (2) PFM 2007
    HeI_emis = ( a + ( b*(np.log(temp)**2.) ) + ( c*np.log(temp) ) + ( d / np.log(temp) ) ) * (1./temp) * 1e-25
    
    # H beta emissivity; from Equation 3.1 of Citation (3) AOS 2010
    Hbeta_emis = ( -2.6584e5 - ( 1420.9*(np.log(temp)**2.) ) + ( 35546*np.log(temp) ) + ( 6.5669e5 / np.log(temp) ) ) \
                                                                * (1./temp) * 1e-25
    
    return HeI_emis / Hbeta_emis

def helium_emissivity_BSS2007(wave, temp, dens):
    '''
    Calculate the emissivity of a HeI line using
    Porter 2005/2007 data, but re-parameterized
    to the form of BSS as given in AOS 2010
    
    Returns the emissivity as a ratio relative
    to H-beta and includes the collisional
    correction for helium
        
    Parameters
    ----------
    wavelength : float
        Wavelength of the HeI line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3), not in log!
        
    Returns
    -------
    1/phi : float
        E(HeI)/E(H-beta) * (1+C/R(HeI))
        Units of emissivity are in ergs*cm^3*s^-1
    '''
    # Redefine the temperature (checked with new HeI emis to confirm that temp needs to be in T4)
    T4 = temp / 1e4
    
    # Find which Helium line we are working with
    HeI_lines = np.array([3889, 4026, 4471, 5016, 5876, 6678, 7065])
    HeI_line = str(HeI_lines[np.where(np.abs(HeI_lines - wave) < 3)[0]][0])

    if HeI_line == '3889':
        phi = 0.8799 * ( T4**(-0.128 - 0.00041*dens) )
    elif HeI_line == '4026':
        phi = 4.233 * ( T4**(0.085 - 0.00012*dens) )
    elif HeI_line == '4471':
        phi = 2.021 * ( T4**(0.121 - 0.00020*dens) )
    elif HeI_line == '5876':
        phi = 0.754 * ( T4**(0.212 - 0.00051*dens) )
    elif HeI_line == '6678':
        phi = 2.639 * ( T4**(0.244 - 0.00054*dens) )
    elif HeI_line == '7065':
        phi = 5.903 * (T4**-0.519) / ( 1.462 - ((0.127 - 0.00076*dens + 0.000000255*dens**2)*T4) )

    return 1. / phi


###############################
# EWs + Underlying Absorption #
###############################
# Linear fit to normalizations to Hb from Equation 5.1 of AOS 2010, excluding downward trend of Ha
fit_balmer_factor = np.polyfit(np.array([4862.721, 4341.684, 4102.891]), np.array([1.00, 0.959, 0.896]), deg=1) # Fit only to Hb, Hg, Hd; Ha causes a turnaround
a_HI_balmer_fit = (fit_balmer_factor[0] * hydrogen_lines) + fit_balmer_factor[1]

# Linear fit to normalizations to HeI4472 from Equation 5.2 of AOS 2010
fit_he_factor = np.polyfit(np.array([7067.198, 6679.994, 5877.299, 4472.755, 4027.328, 3890.151]), np.array([0.4, 0.525, 0.874, 1.0, 1.347, 1.4]), deg=1)
a_He_fit = (fit_he_factor[0] * helium_lines) + fit_he_factor[1]


def stellar_absorption(wave, a_default, ion=None):
    '''
    Calculate the amount of underlying hydrogen
    or helium stellar absorption at the
    wavelength of interest

    Normalizations for optical stellar absorption given in Equations 5.1, 5.2 of AOS2010
    Normalizations for optical+NIR stellar absorption given in Equations 4.2, 4.3 of AOS2015

    Parameters
    ----------
    wave : float
        Wavelength of the Balmer line (in Angstroms)
    a_default : float
        The amount of stellar absorption (in Angstroms)
        at the 'default' wavelength. Default
        wavelength for hydrogen is at Hb and
        at HeI4472 for helium
    ion : str
        Ion of interest, hydrogen or helium
        Accepted formats:
        hydrogen, Hydrogen, H
        helium, Helium, He

    Returns
    -------
    a_at_wave : float
        Amount of stellar absorption at
        the input wavelength
    '''
    # Underlying HI stellar absorption
    if ion in ['hydrogen', 'Hydrogen', 'H']:
        # Match to closest Balmer line
        H_idx = np.where(np.abs(hydrogen_lines - wave) < 3.5)[0][0]
        if H_idx == 0: # P-gamma
            a_at_wave = 0.4 * a_default # From b/t Eq. 4.1, 4.2 of AOS2015
        elif H_idx == 1: # H-alpha
            a_at_wave = 0.942 * a_default
        elif H_idx == 2: # H-beta
            a_at_wave = a_default
        elif H_idx == 3: # H-gamma
            a_at_wave = 0.959 * a_default
        elif H_idx == 4: # H-delta
            a_at_wave = 0.896 * a_default
        elif H_idx == 5: # H8
            a_at_wave = a_HI_balmer_fit[H_idx]*a_default
        else:
            print ('Cannot identify this hydrogen line for a stellar absorption at this wavelength')
            pdb.set_trace()

    # Underlying HeI stellar absorption
    elif ion in ['helium', 'Helium', 'He']:
        # Match to closest Helium line
        He_idx = np.where(np.abs(helium_lines - wave) < 3.5)[0][0]

        # Multiply underlying stellar absorption by normalization
        if He_idx == 0: # HeI10833
            a_at_wave = 0.800 * a_default # From end of Section 2 in AOS2015
        elif He_idx == 1:  # HeI7067
            a_at_wave = 0.400 * a_default
        elif He_idx == 2:  # HeI6679
            a_at_wave = 0.525 * a_default
        elif He_idx == 3:  # HeI5877
            a_at_wave = 0.874 * a_default
        elif He_idx == 4:  # HeI5017
            a_at_wave = a_He_fit[He_idx] * a_default
        elif He_idx == 5:  # HeI4472
            a_at_wave = a_default
        elif He_idx == 6:  # HeI4027
            a_at_wave = 1.347 * a_default
        elif He_idx == 7:  # HeI3890
            a_at_wave = 1.400 * a_default
        else:
            print ('Cannot identify this helium line for a stellar absorption at this wavelength')
            pdb.set_trace()
    else:
        print ('Please supply ion of interest, either hydrogen or helium')

    return a_at_wave


##########################
# Optical Depth Function #
##########################
def optical_depth_function(wave, temp, dens, tau):
    '''
    Determine the wavelength-dependent
    expression for calculating the
    helium optical depth

    Parameters
    ----------
    wave : float
        Wavelength of the Balmer line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3)
    tau : float
        Optical depth at HeI3889

    Returns
    -------
    f_tau : float
        Optical depth function
    '''
    # Redefine the temperature
    T4 = temp / 10000.

    # Match wavelength to relevant column in optical depth table
    idx = np.where(np.abs(helium_optical_depth['Wave'] - wave) < 3.5)[0]

    if len(idx) == 0:
        # HeI5017 is not on the helium_optical_depth table, but its optical depth is 1
        if np.abs(5017.079 - wave) < 3.5:
            f_tau = 1
        else:
            print ('No expression for optical depth at this wavelength')
    # HeI7067 requires a different functional form, given in BSS2002 Section 3.2, paragraph 4 (before Eq. 5)
    elif idx == 9: # 9th index corresponds to the row in Table helium_optical_depth that HeI7067 is in
        f_tau = 1 + ((tau/2) * (0.359 + ((-3.46e-2 - (1.84e-4 * dens) + (3.039e-7 * dens ** 2)) * T4)))
    # HeI10833 optical depth formulation is introducedin AOS2015, Eq. 2.2; I added it into the helium_optical_depth table for consistency in idx retrieval
    elif idx == 11: # 11th index corresponds to the row in Table helium_optical_depth that HeI10833 is in
        f_tau = 1 + ((tau/2) * (0.0149 + ((4.45e-3 - (6.34e-5 * dens) + (9.20e-8 * dens ** 2)) * T4)))
    else:
        a = helium_optical_depth['a'][idx][0]
        b0 = helium_optical_depth['b0'][idx][0]
        b1 = helium_optical_depth['b1'][idx][0]
        b2 = helium_optical_depth['b2'][idx][0]
        f_tau = 1 + ((tau/2) * (a + ((b0 + (b1 * dens) + (b2 * dens ** 2)) * T4)))

    return f_tau


################################
# Collisional to Recombination #
################################

# Hydrogen
# --------
#grid_temp_A2002 = np.array([5802.26125726, 11604.52251451, 34813.56754354, 58022.61257257, 116045.22514514, 174067.83771772, 232090.45029029, 290113.06286286])
# The following are collision strengths for various 1s-->n'l' transitions from Anderson et al. 2002, at T = [5802.26125726, 11604.52251451, 34813.56754354]
#upsilon_3s1s = np.array([0.0651, 0.0696, 0.0776])
#upsilon_3p1s = np.array([0.112, 0.126, 0.186])
#upsilon_3d1s = np.array([0.0621, 0.0658, 0.0782])
#upsilon_4s1s = np.array([0.0223, 0.0255, 0.0319])
#upsilon_4p1s = np.array([0.0402, 0.0479, 0.074])
#upsilon_4d1s = np.array([0.03, 0.0319, 0.0404])
#upsilon_4f1s = np.array([0.0123, 0.0114, 0.0105])
#upsilon_5s1s = np.array([0.0145, 0.0172, 0.0192])
#upsilon_5p1s = np.array([0.0269, 0.0315, 0.0404])
#upsilon_5d1s = np.array([0.0208, 0.0222, 0.0247])
#upsilon_5f1s = np.array([0.00919, 0.00914, 0.00952])
#upsilon_5g1s = np.array([0.00466, 0.00403, 0.00285])
# These are the resulting coefficients for the 2nd degree polyfits to the above collision strengths
upsilon_3s1s_coeff = np.array([0.00233701, -0.00334599,  0.04458994]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_3s1s, 2)
upsilon_3p1s_coeff = np.array([0.10184036, -0.75072221,  1.49488154]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_3p1s, 2)
upsilon_3d1s_coeff = np.array([0.01760335, -0.12551182,  0.28513044]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_3d1s, 2)
upsilon_4s1s_coeff = np.array([0.00357721, -0.01737304,  0.03701513]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_4s1s, 2)
upsilon_4p1s_coeff = np.array([0.03742747, -0.26741176,  0.51648268]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_4p1s, 2)
upsilon_4d1s_coeff = np.array([0.01478313, -0.10941403,  0.23239221]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_4d1s, 2)
upsilon_4f1s_coeff = np.array([0.001418  , -0.0140902 ,  0.04524426]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_4f1s, 2)
upsilon_5s1s_coeff = np.array([-0.00613942,  0.05702998, -0.11317503]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_5s1s, 2)
upsilon_5p1s_coeff = np.array([0.00433421, -0.01864829,  0.03569204]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_5p1s, 2)
upsilon_5d1s_coeff = np.array([0.000757  , -0.00127525,  0.01487691]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_5d1s, 2)
upsilon_5f1s_coeff = np.array([0.00123696, -0.00984927,  0.02873762]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_5f1s, 2)
upsilon_5g1s_coeff = np.array([-0.00048879,  0.00173353,  0.00505922]) #np.polyfit(np.log10(grid_temp_A2002[0:3]), upsilon_5g1s, 2)

def hydrogen_collision_to_recomb(xi, wave, temp, method='AOS2010'):
    '''
    Calculate the factor that corrects the
    measured hydrogen flux for emission due
    to collisional excitation of neutral
    hydrogen

    Assumes that at these densities
    and temperatures, all neutral hydrogen is
    excited from the ground state

    Uses collision strengths from Anderson
    et al. 2002 and branching ratios from
    Omidvar 1983. Recombination rates are
    from Hummer & Storey 1987.

    Parameters
    ----------
    xi : float
        n(HI)/n(HII); ratio of neutral hydrogen
        to ionized hydrogen densities
    wave : float
        Wavelength of the Balmer line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    method : string (optional)
        Method of calculating the C/R ratio;
        options include
        AOS2010:
            - collision strength: Anderson et al. 2002
            - branching ratio: Omidvar 1983
            - recombination rate: Hummer & Storey 1987
        R2015:
            - collision strength: Chianti database
            - recombination rate: Seaton 1959a
        A2002:
            - collision strength: Anderson et al. 2002
            - branching ratio: Omidvar 1983
            - recombination rate: calculated from Peter Storey's 'S18' hydrogen emissivities
    Returns
    -------
    hydrogen_CR : float
        Relative amount of collisional to
        recombination emission for a given
        hydrogen line:
        C/R(wavelength) = eta*K_eff/alpha_eff
    '''
    # Boltzmann's constant
    kB = 8.61733e-5
    # Redefine the temperature
    T4 = temp / 10000.

    # Identify hydrogen line of interest
    idx = np.where(np.abs(hydrogen_lines - wave) < 3.5)[0][0]

    if method == 'AOS2010':
        if idx == 0:
            line = str('Pg')
        elif idx == 1:
            line = str('Ha')
        elif idx == 2:
            line = str('Hb')
        elif idx == 3:
            line = str('Hg')
        elif idx == 4:
            line = str('Hd')
        elif idx == 5:
            line = ('H8')
        else:
            print ('No C/R information for ', wave, 'from AOS2010')
            pdb.set_trace()
        #    print ('Hydrogen C/R for', line)

        rows = np.where(line == hydrogen_CR_coeff['Line'])[0]

        # Calculate the total K_eff/alpha_eff for relevant energy levels -- collisional sum includes an infinite
        # number of levels, but probabilities fall off quickly. This sum excludes terms contributing < 1%
        Keff_alphaeff = 0.
        for i in range(1, 9):  # 1-9 here is to grab the 'Term1', 'Term2', etc. column names
            a, b, c = hydrogen_CR_coeff['Term ' + str(i)][rows]
            Keff_alphaeff += (a * np.exp(b / T4) * (T4 ** c))

        # Amount of collisional to recombination emission; from Equation 6.1 of AOS 2010
        hydrogen_CR = Keff_alphaeff * xi * 1e4

    elif method == 'A2002':
        hc = 1.986445824171758e-18  # Planck constant * speed of light [ergs * m]

        #print ('Using A+2002 collision strengths, O1983 branching ratios, HS1987 recombination rates')
        # Temperatures and reported recombination rates from Hummer & Storey 1987
        #grid_temp = np.array([5000., 7500., 10000., 12500., 15000., 20000., 30000.])
        #recomb_42 = np.array([5.380e-14, 3.863e-14, 3.022e-14, 2.482e-14, 2.105e-14, 1.610e-14, 1.087e-14])  # recombination rate for 4-->2 transition

        if idx == 0: # Pgamma
            # This would be the scaling from recomb_42 to Pgamma
            # scale = np.array([9.87e-2, 9.39e-2, 9.04e-2, 8.77e-2, 8.56e-2, 8.23e-2, 7.79e-2])
            # But instead we use Pbeta recombination scaling, branching ratio, collision strengths, then apply an Energy difference scaling factor to Pgamma
            #scale = np.array([1.84e-1, 1.72e-1, 1.63e-1, 1.57e-1, 1.52e-1, 1.45e-1, 1.36e-1])
            i = np.array([5, 5, 5, 5, 5])
            # Pbeta: 1-->5s, 1-->5p, 1-->5d, 1-->5f 1-->5g
            # 1-->11, 1-->12, 1-->13, 1-->14, 1-->15
            Ediff_factor = np.exp((-13.6 * ((1 / 5 ** 2) - (1 / 6 ** 2))) / (kB * temp)) # Energy difference b/t Pb and Pg
            upsilon = np.array([np.array([np.polyval(upsilon_5s1s_coeff, np.log10(temp)), np.polyval(upsilon_5p1s_coeff, np.log10(temp)), \
                                          np.polyval(upsilon_5d1s_coeff, np.log10(temp)), np.polyval(upsilon_5f1s_coeff, np.log10(temp)), \
                                          np.polyval(upsilon_5g1s_coeff, np.log10(temp))])])
            branching_ratio = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
            emiss_transitions = np.array([hydrogen_emissivity_S2018(hydrogen_lines[7], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[8], temp, 1, ratio=False)])
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([( hc / (hydrogen_lines[7]*1e-10)),
                                       ( hc / (hydrogen_lines[8]*1e-10))])[:, None]
        elif idx == 1:  # Halpha
            # Recombination rate scaling factors from Hummer & Storey 1987
            #scale = np.array([3.04, 2.93, 2.86, 2.82, 2.79, 2.75, 2.70])
            i = np.array([3, 3, 3, 4, 4, 4, 4])
            # Collision strengths from Table 1 of Anderson et al. 2002, interpolated to be at our temperature
            # Halpha: 1-->3s, 1-->3p, 1-->3d, 1-->4s, 1-->4p, 1-->4d, 1-->4f
            # In A2000's definition of i and j, these correspond to: 1-->4, 1-->5, 1-->6, 1-->7, 1-->8, 1-->9, 1-->10
            Ediff_factor = 1.
            upsilon = np.array([np.polyval(upsilon_3s1s_coeff, np.log10(temp)), np.polyval(upsilon_3p1s_coeff, np.log10(temp)),
                                np.polyval(upsilon_3d1s_coeff, np.log10(temp)),
                                np.polyval(upsilon_4s1s_coeff, np.log10(temp)), np.polyval(upsilon_4p1s_coeff, np.log10(temp)),
                                np.polyval(upsilon_4d1s_coeff, np.log10(temp)), np.polyval(upsilon_4f1s_coeff, np.log10(temp))])
            # Branching Ratios from Table 2 of Omidvar 1983
            branching_ratio = np.array([1.0, 1.0, 1.0, 4.16e-1, 4.2e-2, 2.54e-1, 1.0])  # last 4 are from 40,41,42,43 --> 3
            # For Halpha, we care about the 3-->2 recombinations
            emiss_transitions = hydrogen_emissivity_S2018(hydrogen_lines[1], temp, 1, ratio=False)
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([hc / (hydrogen_lines[1]*1e-10)]) # 1e-10 converts Angstrom to meters
        elif idx == 2:  # Hbeta
            #scale = np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00])
            i = np.array([4, 4, 4, 4, 5, 5, 5, 5, 5])
            # Hbeta: 1-->4s, 1-->4p, 1-->4d, 1-->4f, 1-->5s, 1-->5p, 1-->5d, 1-->5f 1-->5g
            # 1-->7, 1-->8, 1-->9, 1-->10, 1-->11, 1-->12, 1-->13, 1-->14, 1-->15
            Ediff_factor = 1.
            upsilon = np.array([np.polyval(upsilon_4s1s_coeff, np.log10(temp)), np.polyval(upsilon_4p1s_coeff, np.log10(temp)),
                      np.polyval(upsilon_4d1s_coeff, np.log10(temp)), np.polyval(upsilon_4f1s_coeff, np.log10(temp)),
                      np.polyval(upsilon_5s1s_coeff, np.log10(temp)), np.polyval(upsilon_5p1s_coeff, np.log10(temp)),
                      np.polyval(upsilon_5d1s_coeff, np.log10(temp)), np.polyval(upsilon_5f1s_coeff, np.log10(temp)),
                      np.polyval(upsilon_5g1s_coeff, np.log10(temp))])
            branching_ratio = np.array([1.0, 1.0, 1.0, 1.0, 2.27e-1, 2.20e-2, 1.07e-1, 3.63e-1, 1.0])  # last 5 are from 50,51,52,53,54 --> 4
            # For Hbeta, we care about the 4--3, 4-->2 recombinations
            emiss_transitions = np.array([hydrogen_emissivity_S2018(hydrogen_lines[2], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[6], temp, 1, ratio=False)])
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([( hc / (hydrogen_lines[2]*1e-10)),
                                       ( hc / (hydrogen_lines[6]*1e-10))])[:, None]
        elif idx == 3:  # Hgamma
            #scale = np.array([4.58e-1, 4.65e-1, 4.68e-01, 4.71e-1, 4.73e-1, 4.75e-1, 4.78e-1])
            i = np.array([5, 5, 5, 5, 5])
            # Hgamma: 1-->5s, 1-->5p, 1-->5d, 1-->5f 1-->5g
            # 1-->11, 1-->12, 1-->13, 1-->14, 1-->15
            Ediff_factor = 1.
            upsilon = np.array([np.array([np.polyval(upsilon_5s1s_coeff, np.log10(temp)), np.polyval(upsilon_5p1s_coeff, np.log10(temp)),
                                          np.polyval(upsilon_5d1s_coeff, np.log10(temp)), np.polyval(upsilon_5f1s_coeff, np.log10(temp)),
                                          np.polyval(upsilon_5g1s_coeff, np.log10(temp))])])
            branching_ratio = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
            # For Hgamma, we care about the 5--4, 5-->3, 5-->2 recombinations
            emiss_transitions = np.array([hydrogen_emissivity_S2018(hydrogen_lines[7], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[8], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[3], temp, 1, ratio=False)])
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([( hc / (hydrogen_lines[7]*1e-10)),
                                       ( hc / (hydrogen_lines[8]*1e-10)),
                                       ( hc / (hydrogen_lines[3]*1e-10))])[:, None]
        elif idx == 4:  # Hdelta
            #scale = np.array([2.51e-1, 2.56e-1, 2.59e-1, 2.61e-1, 2.62e-1, 2.64e-1, 2.66e-1]) # This scales recomb_42 to Hdelta, but since we are applying an Ediff_factor from Hg to Hd, we don't need this
            #scale = np.array([4.58e-1, 4.65e-1, 4.68e-01, 4.71e-1, 4.73e-1, 4.75e-1, 4.78e-1])
            i = np.array([5, 5, 5, 5, 5])
            Ediff_factor = np.exp((-13.6 * ((1 / 5 ** 2) - (1 / 6 ** 2))) / (kB * temp))  # Energy difference b/t Hg and Hd
            upsilon = np.array([np.array([np.polyval(upsilon_5s1s_coeff, np.log10(temp)), np.polyval(upsilon_5p1s_coeff, np.log10(temp)),
                                          np.polyval(upsilon_5d1s_coeff, np.log10(temp)), np.polyval(upsilon_5f1s_coeff, np.log10(temp)),
                                          np.polyval(upsilon_5g1s_coeff, np.log10(temp))])])
            branching_ratio = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
            emiss_transitions = np.array([hydrogen_emissivity_S2018(hydrogen_lines[7], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[8], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[3], temp, 1, ratio=False)])
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([( hc / (hydrogen_lines[7]*1e-10)),
                                       ( hc / (hydrogen_lines[8]*1e-10)),
                                       ( hc / (hydrogen_lines[3]*1e-10))])[:, None]
        elif idx == 5:  # H8
            #scale = np.array([1.02e-1, 1.04e-1, 1.05e-1, 1.06e-1, 1.06e-1, 1.07e-1, 1.08e-1]) # This scales recomb_42 to H8, but since we are applying an Ediff_factor from Hg to Hd, we don't need this
            #scale = np.array([4.58e-1, 4.65e-1, 4.68e-01, 4.71e-1, 4.73e-1, 4.75e-1, 4.78e-1])
            i = np.array([5, 5, 5, 5, 5])
            Ediff_factor = np.exp((-13.6 * ((1 / 5 ** 2) - (1 / 8 ** 2))) / (kB * temp))  # Energy difference b/t Hg and H8
            upsilon = np.array([np.array([np.polyval(upsilon_5s1s_coeff, np.log10(temp)), np.polyval(upsilon_5p1s_coeff, np.log10(temp)),
                                          np.polyval(upsilon_5d1s_coeff, np.log10(temp)), np.polyval(upsilon_5f1s_coeff, np.log10(temp)),
                                          np.polyval(upsilon_5g1s_coeff, np.log10(temp))])])
            branching_ratio = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
            emiss_transitions = np.array([hydrogen_emissivity_S2018(hydrogen_lines[7], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[8], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[3], temp, 1, ratio=False)])
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([( hc / (hydrogen_lines[7]*1e-10)),
                                       ( hc / (hydrogen_lines[8]*1e-10)),
                                       ( hc / (hydrogen_lines[3]*1e-10))])[:, None]
        else:
            print ('Not ready for this hydrogen transition')
            pdb.set_trace()

        #### ****************
        #### May want to do a different type of interpolation in the future; right now at specific grid_temps, we are not recovering the exact recombination rate
        #### ****************
        #grid_alpha = recomb_42 * scale
        #coeff = np.polyfit(np.log10(grid_temp), np.log10(grid_alpha), 2)
        #alpha = 10 ** np.polyval(coeff, np.log10(temp))

        if len(emiss_transitions.shape) == 1:
            alpha = emiss_transitions / photon_energy
        elif len(emiss_transitions.shape) == 2:
            alpha = np.sum(emiss_transitions / photon_energy, axis=0)

        K = 4.004e-8 * np.sqrt(1 / (kB * temp)) * np.exp(-13.6 * (1 - (1 / i ** 2)) / (kB * temp)) * upsilon
        numerator = K * branching_ratio

        hydrogen_CR = xi * Ediff_factor * np.sum(numerator) / alpha[0]


    elif method == 'R2015':
        hc = 1.986445824171758e-18 # Planck constant * speed of light [ergs * m]

        #print ('Using R+2015 formulations for collision strengths and recombination rates')
        if idx == 0:  # Pgamma
            # R+2015 only has fits for Balmer lines; taking coeffs from Hgamma and scaling it to Pgamma with Ediff_factor
            i = 5
            acoeffs = np.array([0.0773, 0.0678, -0.0945, 0.0796, -0.0177, 0.0013])
            bcoeffs = np.array([-13.6820, -0.8629, -0.1957, -0.0375, 0.0199])
            Ediff_factor = np.exp( ( -13.6 * (-19/150) ) / ( kB * temp ) ) # Energy difference b/t Hg and Pg
        elif idx == 1:  # Halpha
            i = 3
            # Fits to collision strength from Table 1 of R+2015
            acoeffs = np.array([0.2500, 0.2461, 0.3297, 0.3892, -0.0928, 0.0071])
            # Fits for the recombination coefficients from Table 2 of R+2015
            #bcoeffs = np.array([-13.3377, -0.7161, -0.1435, -0.0386, 0.0077])
            Ediff_factor = 1.0
            # For Halpha, we care about the 3-->2 recombinations
            emiss_transitions = hydrogen_emissivity_S2018(hydrogen_lines[1], temp, 1, ratio=False)
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([hc / (hydrogen_lines[1]*1e-10)]) # 1e-10 converts Angstrom to meters

        elif idx == 2:  # Hbeta
            i = 4
            acoeffs = np.array([0.1125, 0.1370, -0.1152, 0.1209, -0.0276, 0.0020])
            #bcoeffs = np.array([-13.5225, -0.7928, -0.1749, -0.0412, 0.0154])
            Ediff_factor = 1.0
            # For Hbeta, we care about the 4--3, 4-->2 recombinations
            emiss_transitions = np.array([hydrogen_emissivity_S2018(hydrogen_lines[2], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[6], temp, 1, ratio=False)])
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([( hc / (hydrogen_lines[2]*1e-10)),
                                       ( hc / (hydrogen_lines[6]*1e-10))])[:, None]
        elif idx == 3:  # Hgamma
            i = 5
            acoeffs = np.array([0.0773, 0.0678, -0.0945, 0.0796, -0.0177, 0.0013])
            #bcoeffs = np.array([-13.6820, -0.8629, -0.1957, -0.0375, 0.0199])
            Ediff_factor = 1.0
            # For Hgamma, we care about the 5--4, 5-->3, 5-->2 recombinations
            emiss_transitions = np.array([hydrogen_emissivity_S2018(hydrogen_lines[7], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[8], temp, 1, ratio=False),
                                          hydrogen_emissivity_S2018(hydrogen_lines[3], temp, 1, ratio=False)])
            photon_energy = np.ones(emiss_transitions.shape)
            photon_energy *= np.array([( hc / (hydrogen_lines[7]*1e-10)),
                                       ( hc / (hydrogen_lines[8]*1e-10)),
                                       ( hc / (hydrogen_lines[3]*1e-10))])[:, None]
        elif idx == 4: # Hdelta
            i = 6
            acoeffs = np.array([0.0773, 0.0678, -0.0945, 0.0796, -0.0177, 0.0013])
            bcoeffs = np.array([-13.6820, -0.8629, -0.1957, -0.0375, 0.0199])
            Ediff_factor = np.exp((-13.6 * ((1 / 5 ** 2) - (1 / 6 ** 2))) / (kB * temp))  # Energy difference b/t Hg and Hd
        elif idx == 5: # H8
            i = 6
            acoeffs = np.array([0.0773, 0.0678, -0.0945, 0.0796, -0.0177, 0.0013])
            bcoeffs = np.array([-13.6820, -0.8629, -0.1957, -0.0375, 0.0199])
            Ediff_factor = np.exp((-13.6 * ((1 / 5 ** 2) - (1 / 8 ** 2))) / (kB * temp))  # Energy difference b/t Hg and H8
        else:
            print("Line not ready yet")
            pdb.set_trace()

        # Energy difference from ground state to collisionall excited state
        ediff = -13.6 * (1 - (1 / i ** 2))
        # Reparameterization of the temperature, defined after Eq. A11 of R+2015
        tval = np.log10(temp / 1.0E4)

        # Omega, the collision strength, calculated using Eq. A11 of R+2015
        omega1k = np.zeros(tval.size)

        for aa in range(acoeffs.size):
            omega1k += acoeffs[aa] * tval ** aa
        # Note 4.004E-8/np.sqrt(kB) = 0.5 * 8.629E-6  (i.e. Erik's Eq. 6.3 = Raga's Eq. A12, for the constant)
        q1k = 0.5 * 8.629E-6 * omega1k[0] * np.exp(ediff / (kB * temp)) / np.sqrt(temp)

        # Alpha, the recombination coefficient, calculated using Eq. A13 of R+2015

        if idx in (1,2,3):
            if len(emiss_transitions.shape) == 1:
                alphak = emiss_transitions / photon_energy
            elif len(emiss_transitions.shape) == 2:
                alphak = np.sum(emiss_transitions / photon_energy, axis=0)
        else:
            alphak = np.zeros(tval.size)

            for aa in range(bcoeffs.size):
                alphak += bcoeffs[aa] * tval ** aa
            alphak = 10.0 ** (alphak[0])

        # C/R should be eta * q/alpha
        hydrogen_CR = Ediff_factor * xi * q1k / alphak

    return hydrogen_CR

# Helium
# ------
def helium_collision_to_recomb(wave, temp, dens):
    '''
    Calculate the factor that corrects the
    measured helium flux for emission due
    to collisional excitation of neutral
    helium

    Assumes that at these densities
    and temperatures, all neutral helium is
    excited from the ground state

    Parameters
    ----------
    wave : float
        Wavelength of the HeI line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3)

    Returns
    -------
    helium_CR : float
        Relative amount of collisional to
        recombination emission for a given
        helium line
    '''
    # Redefine the temperature
    T4 = temp / 10000.

    # Find which Upper Level this wavelength comes from
    wave_idx = np.where(np.abs(heliumcoeff['Wavelength'] - wave) < 3.5)[0]
    upper_level = heliumcoeff['Upper Level'][wave_idx]
    #    print ('Upper Level of HeI', wave, 'is', upper_level)

    # Match coefficients for C/R calculation
    idx = np.where(helium_CR_coeff['n 2S+1L (Upper)'] == upper_level)[0]

    # Calculate the summation term in Equation A2 of PFM 2007 Citation (2)
    sum_upper_lev = np.sum(helium_CR_coeff['ai'][idx] * (T4 ** helium_CR_coeff['bi'][idx]) * np.exp(helium_CR_coeff['ci'][idx] / T4))

    # Amount of collisional to recombination emission; from Equation A2 of PFM 2007
    helium_CR = ((1 + (3552. * (T4 ** -0.55) / dens))**-1) * sum_upper_lev

    return helium_CR


#############
# Reddening #
#############
# Seaton 1979 extinction curve + interpolation over it
#f_lambda_avg_old = Table.read(path+'/tables/average_extinction_curve', format='ascii', delimiter=' ')
#f_lambda_avg_interp_old = interp.interp1d(f_lambda_avg_old['wavelength'], f_lambda_avg_old['X(x)'])

# Fine interpolation
f_lambda_avg = Table.read(path+'/tables/extinction_table', format='ascii', delimiter=' ')
f_lambda_avg_interp = interp.interp1d(f_lambda_avg['wavelength'], f_lambda_avg['X(x)'])

# Cardelli, Clayton, & Mathis 1989 A(lambda)/A_v calculation
# http://adsabs.harvard.edu/doi/10.1086/167900
def reddening_coefficient(wave):
    '''
    Calculate the reddening coefficient
    at a given wavelength. As a reminder,
    the optical depth at a wavelength,
    tau(lambda), can be reparameterized as
    tau(lambda) = c * f(lambda), such that
    c depends on the star/stellar population
    but f(lambda) is constant for the same
    star/stellar population (Eq. 7.3 of 
    Osterbrock).

    A(lambda) is the albedo of dust particle
    at a given wavelength (from Osterbrock's
    glossary of physical symbols).

    A(lambda)/Av relates to f(lambda) by:
    
    Parameters
    ----------
    wave : float
        Wavelength of interest (in Angstroms)
    '''

    # Change our wavelengths into microns, which is required for the following parameterizations, and define x
    x = 1 / (wave*0.0001)
    
    # Rv dependent extinction law; infrared regime from Eq. 2a, 2b of CCM1989
    if x >= 0.3 and x < 1.1:
        a = 0.574 * x**1.61
        b = -0.527 * x**1.61
    # Optical+NIR regime from Eq. 3a, 3b of CCM1989
    elif x >= 1.1 and x <= 3.3:
        y = x - 1.82
        a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    else:
        print ('No parameterization of the reddening law in this wavelength regime')
        pdb.set_trace()
        
    Rv = 3.1
    Alambda_Av = a + (b/Rv)
    
    return Alambda_Av


################################
# Generate table of model flux #
################################
# General
# -------
def generate_emission_line_ratio(filename, waves, EWs, EW_Hb, y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi, EW_Pg=0., hydrogen_method='S2018'):
    '''
    Generate the predicted flux ratio
    F(lambda)/F(Hbeta)

    Parameters
    ----------
    waves : array or float
        Wavelength of line(s) of interest (in Angstroms)
    EWs : array or float
        Equivalent width of line(s) of interest (in Angstroms)
    EW_Hb : float
        Equivalent width of Hbeta (in Angstroms)
    EW_Pg : float
        Equivalent width of Paschen-gamma (in Angstroms)
    y_plus : float
        He+/H+ ratio; abundance of
        singly ionized helium
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3)
    c_Hb : float
        Amount of reddening, in both
        the system and in the line
        of sight
    a_H : float
        Underlying hydrogen stellar
        absorption (in Angstroms) at Hbeta
    a_He : float
        Underlying helium stellar
        absorption (in Angstroms) at HeI4472
    tau_He : float
        Optical depth at HeI3889
#    n_HI : float
#        Density of neutral hydrogen (cm^-3)
    xi : float
        Ratio of neutral to singly ionized hydrogen density

    Returns
    -------
    flux_ratio : float
        The model emission line flux ratio, relative to Hbeta
    '''
    emis_lines = np.sort(np.concatenate((hydrogen_lines, helium_lines)))[1:] # 1: to remove the duplicate HeI+H8 3890 line

    wavelength = []
    species = []
    flux_ratio = []
    
    dens = 10**log_dens
    xi = 10**log_xi

    collisional_to_recomb_Hbeta = hydrogen_collision_to_recomb(xi, hydrogen_lines[2], temp, method='A2002')
    f_lambda_at_Hbeta = f_lambda_avg_interp(hydrogen_lines[2])

    for w in range(len(waves)):
        # Determine if working with hydrogen or helium line; within 3 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
        nearest_wave = emis_lines[np.where(np.abs(emis_lines - waves[w]) < 3.5)[0]][0]

        print ('Working on ', nearest_wave)
        
        # Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166) and Pg
        if nearest_wave in hydrogen_lines and nearest_wave != 3890.166 and nearest_wave != 10941.082:
            line_species = 'hydrogen'

            if hydrogen_method == 'S2018':
                emissivity_ratio = hydrogen_emissivity_S2018(waves[w], temp, dens)
                collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp, method='A2002')
            elif hydrogen_method == 'HS1987':
                emissivity_ratio = hydrogen_emissivity_HS1987(waves[w], temp, dens)
                collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp, method='A2002') # HS1987 does not include C/R

            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            reddening_function = ( f_lambda_avg_interp(waves[w]) / f_lambda_at_Hbeta ) - 1.

            flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
        # Any HeI line besides the blended HeI+H8 line (HeI at 3890.151) and NIR HeI10833 line
        elif nearest_wave in helium_lines and nearest_wave != 3890.151 and nearest_wave != 10833.306:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity_PFSD2012(waves[w], temp, dens, ratio='storey')
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            reddening_function = ( f_lambda_avg_interp(waves[w]) / f_lambda_at_Hbeta ) - 1.

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_He_at_wave)/(EWs[w]) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
        
        # The blended HeI+H8 line
        elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
            # For purposes of generating fake fluxes, assume fractional contribution of HeI vs H8 to the blended line is 50/50
            frac_of_he = 0.5
            EW_HeI = frac_of_he * EWs[w]
            EW_H8 = (1-frac_of_he) * EWs[w]

            reddening_function = ( f_lambda_avg_interp(waves[w]) / f_lambda_at_Hbeta ) - 1.

            # HeI 3890.151 contribution:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity_PFSD2012(waves[w], temp, dens)
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_HeI + a_He_at_wave)/(EW_HeI) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
            # H8 contribution:
            line_species = 'hydrogen'

            if hydrogen_method == 'S2018':
                emissivity_ratio = hydrogen_emissivity_S2018(waves[w], temp, dens)
                collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp, method='A2002')
            elif hydrogen_method == 'HS1987':
                emissivity_ratio = hydrogen_emissivity_HS1987(waves[w], temp, dens)
                collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp, method='A2002')
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)

            flux += emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_H8 + a_H_at_wave)/(EW_H8) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
        
        # Infrared HeI10830 line
        elif nearest_wave == 10833.306:
            if EW_Pg == 0.:
                print ('Need to specify the EW of P-gamma!')
                pdb.set_trace()

            else:
                # First, theoretical F(Pg)/F(Hb) ratio, aka 'model-dependent scaling' factor
                line_species = 'hydrogen'

                emissivity_ratio = hydrogen_emissivity_S2018(10941.082, temp, dens) # hard-coded Pg wavelength; could also be hydrogen_lines[0]
                a_H_at_wave = stellar_absorption(10941.082, a_H, ion=line_species)
                reddening_function = ( f_lambda_avg_interp(10941.082) / f_lambda_at_Hbeta ) - 1. # hard-coded Pg wavelength; could also be hydrogen_lines[0]
                collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, 10940.082, temp, method='A2002')

                Pg_to_Hb_flux = emissivity_ratio * ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_Pg + a_H_at_wave)/(EW_Pg) ) * \
                                ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                                10**-(reddening_function * c_Hb)

                # Now, theoretical F(HeI10830)/F(Hbeta) ratio
                line_species = 'helium'

                emissivity_ratio = helium_emissivity_PFSD2012(waves[w], temp, dens)
                a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)
                optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)
                reddening_function = ( f_lambda_avg_interp(waves[w]) / f_lambda_at_Hbeta ) - 1.

                HeI10830_to_Hb_flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_He_at_wave)/(EWs[w]) ) ) * \
                        optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                        10**-(reddening_function * c_Hb)

                # This now gives F(HeI10830)/F(Pg), format as measured
                flux = HeI10830_to_Hb_flux / Pg_to_Hb_flux

        else:
            print ('Check your input wavelength -- not a recognized hydrogen or helium line for MCMC analysis')
            pdb.set_trace()

        wavelength.append(waves[w])
        species.append(line_species)
        flux_ratio.append(flux)
        
    flux_table = Table([wavelength, species, np.array(flux_ratio), EWs], names=('Wavelength', 'Species', 'Flux Ratio', 'EW'))
    flux_table.write(path+'/'+filename, format='ascii', overwrite=True)
    
    return