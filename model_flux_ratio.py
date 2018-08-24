import numpy as np
import os
import scipy.interpolate as interp
from astropy.table import Table

# Flux Ratio = Emissivity Ratio * EW+absorption Ratio * Optical Depth * Collisional to Recombination Ratio * Reddening

# Load in tables we'll need
path = os.getcwd()
hydrogen_emis_coeff = Table.read(path+'/tables/AOS_2010_hydrogen_emissivities', format='ascii', delimiter='\t')
helium_emis = Table.read(path+'/tables/porter_2013_helium_emissivities.dat', format='ascii', delimiter='\t')
helium_emis_coeff = Table.read(path+'/tables/PFM_2007_emissivities', format='ascii', delimiter='\t')
helium_optical_depth = Table.read(path+'/tables/benjamin_1999_optical_depths', format='ascii', delimiter='\t')
hydrogen_CR_coeff = Table.read(path+'/tables/aver_2010_hydrogen_CR_coeff', format='ascii', delimiter='\t')
helium_CR_coeff = Table.read(path+'/tables/PFM_2007_helium_CR_coeff', format='ascii', delimiter='\t')

# Vacuum wavelengths of Balmer lines Ha, Hb, Hg, Hd #, Heps, H8, H9, H10, H11, H12
balmer_lines = np.array([6564.612, 4862.721, 4341.684, 4102.891])  # , 3971.195, 3869.81, 3836.472, 3798.976, 3771.701, 3751.217])

# Vacuum wavelengths of Helium lines
helium_lines = np.array([7283.356, 7067.198, 6679.994, 5877.299, 5017.079, 4923.304, 4472.755, 4027.328, 3890.151])


###########
# Emissivity
###########
def hydrogen_emissivity(wave, temp, dens):
    '''
    Calculate the emissivity of a Balmer line
    relative to H(beta).

    Parameters
    ----------
    wave : float
        Wavelength of the Balmer line (in Angstroms)
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
    idx = np.where(np.abs(wave - balmer_lines) < 3)[0]

    if idx == 0:
        line = str('Ha')
    elif idx == 1:
        line = str('Hb')
    elif idx == 2:
        line = str('Hg')
    elif idx == 3:
        line = str('Hd')
    # print ('Emissivity for', line)

    if line == 'Hb':
        Xt = 1.
    else:
        cij = np.array(hydrogen_emis_coeff[line]).reshape((3, 3))
        Xt = 0.

        for i in range(cij.shape[0]):
            for j in range(cij.shape[0]):
                # Balmer emissivity; from Equation in Section 3.1 Hydrogen emission of AOS 2010 Citaion (3)
                Xt += cij[i][j] * (np.log10(T4) ** i) * (np.log10(dens) ** j)

    # print ('Emissivity for', line, 'is', Xt)

    return Xt


# Interpolated He emissivities
HeI_emis_3889 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['3889A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_4026 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4026A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_4471 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4471A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_4922 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4922A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_5016 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5016A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_5876 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5876A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_6678 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['6678A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_7065 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['7065A'].reshape((14, 21)), kx=1, ky=2)
HeI_emis_7281 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['7281A'].reshape((14, 21)), kx=1, ky=2)


def helium_emissivity(wave, temp, dens, ratio=True):
    '''
    Calculate the emissivity of a HeI line.
    Defaults to returning the emissivity as
    a ratio relative to H(beta).

    Parameters
    ----------
    wavelength : float
        Wavelength of the HeI line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3)
    ratio : True/False (optional)
        If False, returns the HeI emissivity value

    Returns
    -------
    emissivity : float
        The E(HeI)/E(H(beta)) ratio (default).
        Units of emissivity are in ergs*cm^3*s^-1
    '''
    # Helium lines in format of column names
    HeI_lines = np.array([3889, 4026, 4471, 4922, 5016, 5876, 6678, 7065, 7281])

    # Find column in Porter's 2013 emissivities corresponding to HeI wavelength of interest
    HeI_line = str(HeI_lines[np.where(np.abs(HeI_lines - wave) < 3)[0]][0])

    if HeI_line == '3889':
        HeI_emis = 10 ** HeI_emis_3889(dens, temp)[0][0]
    elif HeI_line == '4026':
        HeI_emis = 10 ** HeI_emis_4026(dens, temp)[0][0]
    elif HeI_line == '4471':
        HeI_emis = 10 ** HeI_emis_4471(dens, temp)[0][0]
    elif HeI_line == '4922':
        HeI_emis = 10 ** HeI_emis_4922(dens, temp)[0][0]
    elif HeI_line == '5016':
        HeI_emis = 10 ** HeI_emis_5016(dens, temp)[0][0]
    elif HeI_line == '5876':
        HeI_emis = 10 ** HeI_emis_5876(dens, temp)[0][0]
    elif HeI_line == '6678':
        HeI_emis = 10 ** HeI_emis_6678(dens, temp)[0][0]
    elif HeI_line == '7065':
        HeI_emis = 10 ** HeI_emis_7065(dens, temp)[0][0]
    elif HeI_line == '7281':
        HeI_emis = 10 ** HeI_emis_7281(dens, temp)[0][0]

    # H beta emissivity; from Equation 3.1 of Citation (3) AOS 2010
    Hbeta_emis = (-2.6584e5 - (1420.9 * (np.log(temp) ** 2.)) + (35546 * np.log(temp)) + (6.5669e5 / np.log(temp))) \
                 * (1. / temp) * 1e-25

    return HeI_emis / Hbeta_emis


#############################
# EWs + Underlying Absorption
#############################
# Linear fit to normalizations to Hb from Equation 5.1 of AOS 2010, excluding downward trend of Ha
fit_balmer_factor = np.polyfit(balmer_lines[1:4], np.array([1.00, 0.959, 0.896]), deg=1)
a_HI_balmer_fit = (fit_balmer_factor[0] * balmer_lines) + fit_balmer_factor[1]

# Linear fit to normalizations to HeI4472 from Equation 5.2 of AOS 2010
fit_he_factor = np.polyfit(np.array([7067.198, 6679.994, 5877.299, 4472.755, 4027.328, 3890.151]), np.array([0.4, 0.525, 0.874, 1.0, 1.347, 1.4]), deg=1)
a_He_fit = (fit_he_factor[0] * helium_lines) + fit_he_factor[1]


def stellar_absorption(wave, a_default, ion=None):
    '''
    Calculate the amount of underlying hydrogen
    or helium stellar absorption at the
    wavelength of interest

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
        Ion of interest, hydrogen or helium.
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
        H_idx = np.where(np.abs(balmer_lines - wave) < 3)[0]

        if H_idx == 0:
            a_at_wave = 0.942 * a_default
        elif H_idx == 1:
            a_at_wave = a_default
        elif H_idx == 2:
            a_at_wave = 0.959 * a_default
        elif H_idx == 3:
            a_at_wave = 0.896 * a_default

    # Underlying HeI stellar absorption
    elif ion in ['helium', 'Helium', 'He']:
        # Match to closest Helium line
        He_idx = np.abs(wave - helium_lines).argmin()

        # Multiply underlying stellar absorption by normalization
        if He_idx == 0:  # HeI7283
            a_at_wave = a_He_fit[0] * a_default

        elif He_idx == 1:  # HeI7067
            a_at_wave = 0.400 * a_default

        elif He_idx == 2:  # HeI6679
            a_at_wave = 0.525 * a_default

        elif He_idx == 3:  # HeI5877
            a_at_wave = 0.874 * a_default

        elif He_idx == 4:  # HeI5017
            a_at_wave = a_He_fit[4] * a_default

        elif He_idx == 5:  # HeI4923
            a_at_wave = a_He_fit[5] * a_default

        elif He_idx == 6:  # HeI4472
            a_at_wave = a_default

        elif He_idx == 7:  # HeI4027
            a_at_wave = 1.347 * a_default

        elif He_idx == 8:  # HeI3890
            a_at_wave = 1.400 * a_default
    else:
        print ('Please supply ion of interest, either hydrogen or helium')

    return a_at_wave


#######################
# Optical Depth Function
#######################
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
    idx = np.where(np.abs(helium_optical_depth['Wave'] - wave) < 3)[0]

    if len(idx) == 0:
        if np.abs(5017.079 - wave) < 3:
            f_tau = 1
        else:
            print ('No expression for optical depth at this wavelength')
    #    else:
    #        print ('Optical depth function for HeI', helium_optical_depth['Wave'][idx][0])

    # HeI7065 requires a different functional form; given in BSS 2002 Section 3.2, paragraph 4 (before Eq 5)
    elif idx == 9:
        f_tau = 1 + ((tau / 2) * (0.359 + ((-3.46e-2 - (1.84e-4 * dens) + 3.039e-7 * dens ** 2) * T4)))

    else:
        a = helium_optical_depth['a'][idx][0]
        b0 = helium_optical_depth['b0'][idx][0]
        b1 = helium_optical_depth['b1'][idx][0]
        b2 = helium_optical_depth['b2'][idx][0]

        f_tau = 1 + ((tau / 2) * (a + ((b0 + (b1 * dens) + (b2 * dens ** 2)) * T4)))

    return f_tau


#############################
# Collisional to Recombination
#############################
def hydrogen_collision_to_recomb(eta, wave, temp):
    '''
    Calculate the factor that corrects the
    measured hydrogen flux for emission due
    to collisional excitation of neutral
    hydrogen

    Assumes that at these densities
    and temperatures, all neutral hydrogen is
    excited from the ground state

    Parameters
    ----------
    eta : float
        n(HI)/n(HII); ratio of neutral hydrogen
        to ionized hydrogen densities
    wave : float
        Wavelength of the Balmer line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)

    Returns
    -------
    hydrogen_CR : float
        Relative amount of collisional to
        recombination emission for a given
        Balmer line
        C/R(wavelength) = eta*K_eff/alpha_eff
    '''
    # Redefine the temperature
    T4 = temp / 10000.

    # Match Balmer line of interest to relevant rows in Table 3 of AOS 2010
    idx = np.where(np.abs(balmer_lines - wave) < 3)[0]

    if idx == 0:
        line = str('Ha')
    elif idx == 1:
        line = str('Hb')
    elif idx == 2:
        line = str('Hg')
    elif idx == 3:
        line = str('Hd')
    #    print ('Hydrogen C/R for', line)

    rows = np.where(line == hydrogen_CR_coeff['Line'])[0]

    # Calculate the total K_eff/alpha_eff for relevant energy levels -- collisional sum includes an infinite
    # number of levels, but probabilities fall off quickly. This sum excludes terms contributing < 1%
    Keff_alphaeff = 0.
    for i in range(1, 9):
        a, b, c = hydrogen_CR_coeff['Term ' + str(i)][rows]
        Keff_alphaeff += (a * np.exp(-b / T4) * (T4 ** c))
    #    print (Keff_alphaeff)

    # Amount of collisional to recombination emission; from Equation 6.1 of AOS 2010
    hydrogen_CR = eta * Keff_alphaeff

    return hydrogen_CR


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
        Wavelength of the Balmer line (in Angstroms)
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
    wave_idx = np.where(np.abs(helium_emis_coeff['Wavelength'] - wave) < 3)[0]
    upper_level = helium_emis_coeff['Upper Level'][wave_idx]
    #    print ('Upper Level of HeI', wave, 'is', upper_level)

    # Match coefficients for C/R calculation
    idx = np.where(helium_CR_coeff['n 2S+1L (Upper)'] == upper_level)[0]

    # Calculate the summation term in Equation A2 of PFM 2007 Citation (2)
    sum_upper_lev = np.sum(
        helium_CR_coeff['ai'][idx] * (T4 ** helium_CR_coeff['bi'][idx]) * np.exp(helium_CR_coeff['ci'][idx] / T4))

    # Amount of collisional to recombination emission; from Equation A2 of PFM 2007
    helium_CR = (1 + (3552. * (T4 ** -0.55) / dens)) ** -1 * sum_upper_lev

    return helium_CR


###########
# Reddening
###########
# Seaton 1979 extinction curve + interpolation over it
f_lambda_avg = Table.read(/tables/average_extinction_curve', format='ascii', delimiter=' ')
f_lambda_avg_interp = interp.interp1d(f_lambda_avg['wavelength'], f_lambda_avg['X(x)'])