import numpy as np
import os
import pdb
import scipy.interpolate as interp
from astropy.table import Table

# Flux Ratio = Emissivity Ratio * EW+absorption Ratio * Optical Depth * Collisional to Recombination Ratio * Reddening

# Load in tables we'll need
path = os.getcwd()

hydrogen_emis_coeff = Table.read(path+'/tables/hydrogen_emissivity_HS1987', format='ascii', delimiter='\t')
helium_emis = Table.read(path+'/tables/helium_emissivity', format='ascii', delimiter='\t')
helium_emis_coeff = Table.read(path+'/tables/helium_emissivity_coeff', format='ascii', delimiter='\t')
helium_optical_depth = Table.read(path+'/tables/helium_optical_depth', format='ascii', delimiter='\t')
hydrogen_CR_coeff = Table.read(path+'/tables/hydrogen_CR_coeff', format='ascii', delimiter='\t')
helium_CR_coeff = Table.read(path+'/tables/helium_CR_coeff', format='ascii', delimiter='\t')

# Vacuum wavelengths of Balmer lines Ha, Hb, Hg, Hd, Heps, H8, H9, H10, H11, H12
#balmer_lines = np.array([6564.612, 4862.721, 4341.684, 4102.891, 3971.195, 3890.166, 3836.472, 3798.976, 3771.701, 3751.217])
# Vacuum wavelengths of Balmer lines Ha, Hb, Hg, Hd, H8 for MCMC
balmer_lines = np.array([6564.612, 4862.721, 4341.684, 4102.891, 3890.166])

# Vacuum wavelengths of Helium lines for MCMC
#helium_lines = np.array([7283.356, 7067.198, 6679.994, 5877.299, 5017.079, 4923.304, 4472.755, 4027.328, 3890.151])
helium_lines = np.array([10833.306, 7067.198, 6679.994, 5877.299, 5017.079, 4472.755, 4027.328, 3890.151])

##############
# Emissivity #
##############

# Hydrogen
# --------
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
    elif idx == 4:
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
        cij = np.array(hydrogen_emis_coeff[line]).reshape((3, 3))
        Xt = 0.

        for i in range(cij.shape[0]):
            for j in range(cij.shape[1]):
                # Balmer emissivity; from Equation in Section 3.1 Hydrogen emission of AOS 2010 Citation (3)
                Xt += cij[i][j] * (np.log10(T4) ** i) * (np.log10(dens) ** j)

    # print ('Emissivity for', line, 'is', Xt)

    return Xt


# Helium
# ------
# Interpolated Porter 2013 HeI emissivities
HeI_emis_3889 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['3889A'].reshape((14, 21)), kx=1, ky=1)
HeI_emis_4026 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4026A'].reshape((14, 21)), kx=1, ky=1)
HeI_emis_4471 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['4471A'].reshape((14, 21)), kx=1, ky=1)
HeI_emis_5016 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5016A'].reshape((14, 21)), kx=1, ky=1)
HeI_emis_5876 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['5876A'].reshape((14, 21)), kx=1, ky=1)
HeI_emis_6678 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['6678A'].reshape((14, 21)), kx=1, ky=1)
HeI_emis_7065 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['7065A'].reshape((14, 21)), kx=1, ky=1)
HeI_emis_10833 = interp.RectBivariateSpline(np.logspace(1, 14, num=14), np.linspace(5e3, 25e3, num=21), helium_emis['10830A'].reshape((14, 21)), kx=1, ky=1)


def helium_emissivity(wave, temp, dens, ratio=True):
    '''
    Calculate the emissivity of a HeI line
    using Porter et al.'s 2013 emissivities

    Defaults to returning the emissivity as
    a ratio relative to H(beta) and does not
    include the collisional correction for HeI

    Parameters
    ----------
    wavelength : float
        Wavelength of the HeI line (in Angstroms)
    temp : float
        Temperature of the gas (in Kelvin)
    dens : float
        Density of the gas (cm^-3), not in log!
    ratio : True/False (optional)
        If False, returns the HeI emissivity value

    Returns
    -------
    emissivity : float
        The E(HeI)/E(H(beta)) ratio (default).
        Units of emissivity are in ergs*cm^3*s^-1
    '''
    # Helium lines, in format of column names
    HeI_lines = np.array([3889, 4026, 4471, 5016, 5876, 6678, 7065, 10833])

    # Find column in Porter's 2013 emissivities corresponding to HeI wavelength of interest
    HeI_line = str(HeI_lines[np.where(np.abs(HeI_lines - wave) < 3.5)[0]][0])

    if HeI_line == '3889':
        HeI_emis = 10 ** HeI_emis_3889(dens, temp)[0][0]
    elif HeI_line == '4026':
        HeI_emis = 10 ** HeI_emis_4026(dens, temp)[0][0]
    elif HeI_line == '4471':
        HeI_emis = 10 ** HeI_emis_4471(dens, temp)[0][0]
    elif HeI_line == '5016':
        HeI_emis = 10 ** HeI_emis_5016(dens, temp)[0][0]
    elif HeI_line == '5876':
        HeI_emis = 10 ** HeI_emis_5876(dens, temp)[0][0]
    elif HeI_line == '6678':
        HeI_emis = 10 ** HeI_emis_6678(dens, temp)[0][0]
    elif HeI_line == '7065':
        HeI_emis = 10 ** HeI_emis_7065(dens, temp)[0][0]
    elif HeI_line == '10833':
        HeI_emis = 10 ** HeI_emis_10833(dens, temp)[0][0]

    # H beta emissivity; from Equation 3.1 of Citation (3) AOS 2010
    Hbeta_emis = (-2.6584e5 - (1420.9 * (np.log(temp) ** 2.)) + (35546 * np.log(temp)) + (6.5669e5 / np.log(temp))) \
                 * (1 / temp) * 1e-25

    return HeI_emis / Hbeta_emis

def helium_emissivity_2007(wave, temp):
    '''
    Calculate the emissivity of a HeI line
    using Porter et al.'s 2007 work

    Defaults to returning the emissivity as
    a ratio relative to H(beta) and does not
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
        The E(HeI)/E(H(beta)) ratio (default).
        Units of emissivity are in ergs*cm^3*s^-1
    '''
    # Find the row in Porter's 2007 emissivities corresponding to HeI wavelength of interest
    idx = np.abs(helium_emis_coeff['Wavelength'] - wave).argmin()

    a = helium_emis_coeff['a'][idx]
    b = helium_emis_coeff['b'][idx]
    c = helium_emis_coeff['c'][idx]
    d = helium_emis_coeff['d'][idx]
    
    # HeI emissivity; from Equation A1 of Citation (2) PFM 2007
    HeI_emis = ( a + ( b*(np.log(temp)**2.) ) + ( c*np.log(temp) ) + ( d / np.log(temp) ) ) * (1./temp) * 1e-25
    
    # H beta emissivity; from Equation 3.1 of Citation (3) AOS 2010
    Hbeta_emis = ( -2.6584e5 - ( 1420.9*(np.log(temp)**2.) ) + ( 35546*np.log(temp) ) + ( 6.5669e5 / np.log(temp) ) ) \
                                                                * (1./temp) * 1e-25
    
    return HeI_emis / Hbeta_emis

def helium_emissivity_2007_BSS(wave, temp, dens):
    '''
    Calculate the emissivity of a HeI line using
    Porter 2005/2007 data, but re-parameterized
    to the form of BSS as given in AOS 2010
    
    Returns the emissivity as a ratio relative
    to H(beta) and includes the collisional
    correction for helium
    
    Re-parameterized by AOS.
    
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
        E(HeI)/E(H(beta)) * (1+C/R(HeI))
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
fit_balmer_factor = np.polyfit(np.array([4862.721, 4341.684, 4102.891]), np.array([1.00, 0.959, 0.896]), deg=1) # Fit to Hb, Hg, Hd
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
    # Normalizations for stellar absorption given in Equations 5.1, 5.2 of AOS 2010
    # Underlying HI stellar absorption
    if ion in ['hydrogen', 'Hydrogen', 'H']:
        # Match to closest Balmer line
        H_idx = np.where(np.abs(balmer_lines - wave) < 3.5)[0]

        if H_idx == 0: # Ha
            a_at_wave = 0.942 * a_default
        elif H_idx == 1: # Hb
            a_at_wave = a_default
        elif H_idx == 2: # Hg
            a_at_wave = 0.959 * a_default
        elif H_idx == 3: # Hd
            a_at_wave = 0.896 * a_default
        elif H_idx == 4: # H8
            a_at_wave = a_HI_balmer_fit[4]*a_default
        else:
            print ('Cannot identify this hydrogen line')
            pdb.set_trace()

    # Underlying HeI stellar absorption
    elif ion in ['helium', 'Helium', 'He']:
        # Match to closest Helium line
        He_idx = np.where(np.abs(helium_lines - wave) < 3.5)[0]

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
            a_at_wave = a_He_fit[3] * a_default
        elif He_idx == 5:  # HeI4472
            a_at_wave = a_default
        elif He_idx == 6:  # HeI4027
            a_at_wave = 1.347 * a_default
        elif He_idx == 7:  # HeI3890
            a_at_wave = 1.400 * a_default
        else:
            print ('Cannot identify this helium line')
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

    # HeI5017 is not on the helium_optical_depth table, but its optical depth is 1
    if len(idx) == 0:
        if np.abs(5017.079 - wave) < 3:
            f_tau = 1
        else:
            print ('No expression for optical depth at this wavelength')

    # HeI7067 requires a different functional form, given in BSS 2002 Section 3.2, paragraph 4 (before Eq 5)
    elif idx == 9: # 9th index corresponds to the row in Table helium_optial_depth that HeI7067 is in
        f_tau = 1 + ((tau/2) * (0.359 + ((-3.46e-2 - (1.84e-4 * dens) + (3.039e-7 * dens ** 2)) * T4)))
    elif idx == 11: # 11th index corresponds to the row in Table helium_optical_depth that HeI10833 is in (I added it in the table; copied from Eq. 2.2 of AOS2015)
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
def hydrogen_collision_to_recomb(xi, wave, temp):
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
    xi : float
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
    for i in range(1, 9): # 1-9 here is to grab the 'Term1', 'Term2', etc. column names
        a, b, c = hydrogen_CR_coeff['Term ' + str(i)][rows]
        Keff_alphaeff += (a * np.exp(b / T4) * (T4 ** c))
    #    print (Keff_alphaeff)

    # Amount of collisional to recombination emission; from Equation 6.1 of AOS 2010
    hydrogen_CR = xi * Keff_alphaeff

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
    helium_CR = ((1 + (3552. * (T4 ** -0.55) / dens))**-1) * sum_upper_lev

    return helium_CR


#############
# Reddening #
#############
# Seaton 1979 extinction curve + interpolation over it
f_lambda_avg = Table.read(path+'/tables/average_extinction_curve', format='ascii', delimiter=' ')
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
    
    x = 1 / (wave * 0.0001) # change Angstrom wavelength to microns; 10000A = 1micron
    
    y = x - 1.82
    a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
    b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    
    Rv = 3.1
    
    Alambda_AV = a + (b/Rv)
    
    return Alambda_AV


################################
# Generate table of model flux #
################################
def generate_emission_line_ratio(filename, waves, EWs, EW_Hb, y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi):
#def generate_emission_line_ratio(filename, waves, EWs, EW_Hb, y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, n_HI):
    '''
    Generate the predicted flux ratio
    F(λ)/F(Hβ)
    
    Parameters
    ----------
    waves : array or float
        Wavelength of line(s) of interest (in Angstroms)
    EWs : array or float
        Equivalent width of line(s) of interest (in Angstroms)
    EW_Hb : float
        Equivalent width of Hβ (in Angstroms)
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
        absorption (in Angstroms) at Hβ
    a_He : float
        Underlying helium stellar
        absorption (in Angstroms) at HeI λ4472
    tau_He : float
        Optical depth at HeI λ3889
    n_HI : float
        Density of neutral hydrogen
        (cm^-3)
    xi : float
	Ratio of neutral to singly ionized
	hydrogen density
    
    Returns
    -------
    flux_ratio : float
    '''
    emis_lines = np.sort(np.concatenate((balmer_lines, helium_lines)))[1:] # 1: to remove the duplicate HeI+H8 3890 line

    species = []
    wavelength = []
    flux_ratio = []
    
    dens = 10**log_dens
    xi = 10**log_xi

    collisional_to_recomb_Hbeta = hydrogen_collision_to_recomb(xi, balmer_lines[1], temp)
    f_lambda_at_Hbeta = f_lambda_avg_interp(balmer_lines[1])

    for w in range(len(waves)):
        # Determine if working with hydrogen or helium line; within 3 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
        nearest_wave = emis_lines[np.where(np.abs(emis_lines - waves[w]) < 3)[0]][0]
        
        # Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166)
        if nearest_wave in balmer_lines and nearest_wave != 3890.166:
            line_species = 'hydrogen'
            
            emissivity_ratio = hydrogen_emissivity_HS1987(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp) # HS1987 does not include C/R 
            reddening_function = ( f_lambda_avg_interp(waves[w]) / f_lambda_at_Hbeta ) - 1.            

            flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
        # Any HeI line besides the blended HeI+H8 line (HeI at 3890.151)
        elif nearest_wave in helium_lines and nearest_wave != 3890.151:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity(waves[w], temp, dens)            
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            #collisional_to_recomb_ratio = helium_collision_to_recomb(waves[w], temp, dens) # C/R is included in Porter 2012 emissivities
            reddening_function = ( f_lambda_avg_interp(waves[w]) / f_lambda_at_Hbeta ) - 1.

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_He_at_wave)/(EWs[w]) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
        
        # The blended HeI+H8 line
        elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
            # Calculate fractional contribution of HeI and H8 to the blended line
#            frac_of_he = ( (y[0] * ( (EWs[w] + a_H_at_wave) / EWs[w]) ) - ( 0.104 * (temp/1e4)**0.046 * 10**-(reddening_function * c_Hb) ) ) / y[0]
#            EW_HeI = frac_of_he * EWs[w]
#            EW_H8 = (1-frac_of_he) * EWs[w]
            frac_of_he = 0.5
            EW_HeI = frac_of_he * EWs[w]
            EW_H8 = (1-frac_of_he) * EWs[w]

            reddening_function = ( f_lambda_avg_interp(waves[w]) / f_lambda_at_Hbeta ) - 1.

            # HeI 3890.151 contribution:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity(waves[w], temp, dens)
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            #collisional_to_recomb_ratio = helium_collision_to_recomb(waves[w], temp, dens) # C/R is included in Porter 2012 emissivities

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_HeI + a_He_at_wave)/(EW_HeI) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
            # H8 contribution:
            line_species = 'hydrogen'

            emissivity_ratio = hydrogen_emissivity_HS1987(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_factor = np.exp(( -13.6 * ((1/5**2)-(1/8**2)) ) / (8.6173303e-5 * temp)) # scale factor for C/R(Hg) to C/R(H8)
            collisional_to_recomb_ratio = collisional_to_recomb_factor * hydrogen_collision_to_recomb(xi, 4341.684, temp) # Calculate C/R(Hg) and multiply by above scale factor

            flux += emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_H8 + a_H_at_wave)/(EW_H8) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)

        flux_ratio.append(flux)
        species.append(line_species)
        wavelength.append(waves[w])
        
    flux_table = Table([wavelength, species, np.array(flux_ratio), EWs], names=('Wavelength', 'Species', 'Flux Ratio', 'EW'))
    flux_table.write(path+'/'+filename, format='ascii', overwrite=True)
    
    return


def generate_emission_line_ratio_CCMred(filename, waves, EWs, EW_Hb, y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi):
#def generate_emission_line_ratio(filename, waves, EWs, EW_Hb, y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, n_HI):
    '''
    Generate the predicted flux ratio
    F(λ)/F(Hβ)
    
    Parameters
    ----------
    waves : array or float
        Wavelength of line(s) of interest (in Angstroms)
    EWs : array or float
        Equivalent width of line(s) of interest (in Angstroms)
    EW_Hb : float
        Equivalent width of Hβ (in Angstroms)
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
        absorption (in Angstroms) at Hβ
    a_He : float
        Underlying helium stellar
        absorption (in Angstroms) at HeI λ4472
    tau_He : float
        Optical depth at HeI λ3889
    n_HI : float
        Density of neutral hydrogen
        (cm^-3)
    xi : float
    Ratio of neutral to singly ionized
    hydrogen density
    
    Returns
    -------
    flux_ratio : float
    '''
    emis_lines = np.sort(np.concatenate((balmer_lines, helium_lines)))[1:] # 1: to remove the duplicate HeI+H8 3890 line

    species = []
    wavelength = []
    flux_ratio = []
    
    dens = 10**log_dens
    xi = 10**log_xi

    collisional_to_recomb_Hbeta = hydrogen_collision_to_recomb(xi, balmer_lines[1], temp)
    AHbeta_AV = reddening_coefficient(balmer_lines[1])

    for w in range(len(waves)):
        # Determine if working with hydrogen or helium line; within 3 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
        nearest_wave = emis_lines[np.where(np.abs(emis_lines - waves[w]) < 3)[0]][0]
        
        # Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166)
        if nearest_wave in balmer_lines and nearest_wave != 3890.166:
            line_species = 'hydrogen'
            
            emissivity_ratio = hydrogen_emissivity(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp)            
            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.           

            flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
        # Any HeI line besides the blended HeI+H8 line (HeI at 3890.151)
        elif nearest_wave in helium_lines and nearest_wave != 3890.151:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity(waves[w], temp, dens)            
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            #collisional_to_recomb_ratio = helium_collision_to_recomb(waves[w], temp, dens); C/R included in Porter 2012 emissivities
            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.  

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_He_at_wave)/(EWs[w]) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
        
        # The blended HeI+H8 line
        elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
            # Calculate fractional contribution of HeI and H8 to the blended line
            frac_of_he = 0.5
            EW_HeI = frac_of_he * EWs[w]
            EW_H8 = (1-frac_of_he) * EWs[w]

            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.  

            # HeI 3890.151 contribution:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity(waves[w], temp, dens)
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            #collisional_to_recomb_ratio = helium_collision_to_recomb(waves[w], temp, dens); C/R included in Porter 2012 emissivities

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_HeI + a_He_at_wave)/(EW_HeI) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
            # H8 contribution:
            line_species = 'hydrogen'

            emissivity_ratio = hydrogen_emissivity(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_factor = np.exp(( -13.6 * ((1/5**2)-(1/8**2)) ) / (8.6173303e-5 * temp)) # scale factor for C/R(Hg) to C/R(H8)
            collisional_to_recomb_ratio = collisional_to_recomb_factor * hydrogen_collision_to_recomb(xi, 4341.684, temp) # Calculate C/R(Hg) and multiply by above scale factor

            flux += emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_H8 + a_H_at_wave)/(EW_H8) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)

        flux_ratio.append(flux)
        species.append(line_species)
        wavelength.append(waves[w])
        
    flux_table = Table([wavelength, species, np.array(flux_ratio), EWs], names=('Wavelength', 'Species', 'Flux Ratio', 'EW'))
    flux_table.write(path+'/'+filename, format='ascii', overwrite=True)
    
    return


def generate_emission_line_ratio_CCMred_2007emis(filename, waves, EWs, EW_Hb, y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi):
#def generate_emission_line_ratio(filename, waves, EWs, EW_Hb, y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, n_HI):
    '''
    Generate the predicted flux ratio
    F(λ)/F(Hβ)
    
    Parameters
    ----------
    waves : array or float
        Wavelength of line(s) of interest (in Angstroms)
    EWs : array or float
        Equivalent width of line(s) of interest (in Angstroms)
    EW_Hb : float
        Equivalent width of Hβ (in Angstroms)
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
        absorption (in Angstroms) at Hβ
    a_He : float
        Underlying helium stellar
        absorption (in Angstroms) at HeI λ4472
    tau_He : float
        Optical depth at HeI λ3889
    n_HI : float
        Density of neutral hydrogen
        (cm^-3)
    xi : float
    Ratio of neutral to singly ionized
    hydrogen density
    
    Returns
    -------
    flux_ratio : float
    '''
    emis_lines = np.sort(np.concatenate((balmer_lines, helium_lines)))[1:] # 1: to remove the duplicate HeI+H8 3890 line

    species = []
    wavelength = []
    flux_ratio = []
    
    dens = 10**log_dens
    xi = 10**log_xi

    collisional_to_recomb_Hbeta = hydrogen_collision_to_recomb(xi, balmer_lines[1], temp)
    AHbeta_AV = reddening_coefficient(balmer_lines[1])

    for w in range(len(waves)):
        # Determine if working with hydrogen or helium line; within 3 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
        nearest_wave = emis_lines[np.where(np.abs(emis_lines - waves[w]) < 3)[0]][0]
        
        # Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166)
        if nearest_wave in balmer_lines and nearest_wave != 3890.166:
            line_species = 'hydrogen'
            
            emissivity_ratio = hydrogen_emissivity(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp)            
            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.           

            flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
        # Any HeI line besides the blended HeI+H8 line (HeI at 3890.151)
        elif nearest_wave in helium_lines and nearest_wave != 3890.151:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity_2007(waves[w], temp)            
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            collisional_to_recomb_ratio = helium_collision_to_recomb(waves[w], temp, dens)            
            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.  

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_He_at_wave)/(EWs[w]) ) ) * \
                    optical_depth_at_wave * ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
        
        # The blended HeI+H8 line
        elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
            # Calculate fractional contribution of HeI and H8 to the blended line
            frac_of_he = 0.5
            EW_HeI = frac_of_he * EWs[w]
            EW_H8 = (1-frac_of_he) * EWs[w]

            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.  

            # HeI 3890.151 contribution:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity_2007(waves[w], temp)
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            collisional_to_recomb_ratio = helium_collision_to_recomb(waves[w], temp, dens)            

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_HeI + a_He_at_wave)/(EW_HeI) ) ) * \
                    optical_depth_at_wave * ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
            # H8 contribution:
            line_species = 'hydrogen'

            emissivity_ratio = hydrogen_emissivity(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_factor = np.exp(( -13.6 * ((1/5**2)-(1/8**2)) ) / (8.6173303e-5 * temp)) # scale factor for C/R(Hg) to C/R(H8)
            collisional_to_recomb_ratio = collisional_to_recomb_factor * hydrogen_collision_to_recomb(xi, 4341.684, temp) # Calculate C/R(Hg) and multiply by above scale factor

            flux += emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_H8 + a_H_at_wave)/(EW_H8) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)

        flux_ratio.append(flux)
        species.append(line_species)
        wavelength.append(waves[w])
        
    flux_table = Table([wavelength, species, np.array(flux_ratio), EWs], names=('Wavelength', 'Species', 'Flux Ratio', 'EW'))
    flux_table.write(path+'/'+filename, format='ascii', overwrite=True)
    
    return

def generate_emission_line_ratio_CCMred_BSSemis(filename, waves, EWs, EW_Hb, y_plus, temp, log_dens, c_Hb, a_H, a_He, tau_He, log_xi):
#def generate_emission_line_ratio(filename, waves, EWs, EW_Hb, y_plus, temp, dens, c_Hb, a_H, a_He, tau_He, n_HI):
    '''
    Generate the predicted flux ratio
    F(λ)/F(Hβ)
    
    Parameters
    ----------
    waves : array or float
        Wavelength of line(s) of interest (in Angstroms)
    EWs : array or float
        Equivalent width of line(s) of interest (in Angstroms)
    EW_Hb : float
        Equivalent width of Hβ (in Angstroms)
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
        absorption (in Angstroms) at Hβ
    a_He : float
        Underlying helium stellar
        absorption (in Angstroms) at HeI λ4472
    tau_He : float
        Optical depth at HeI λ3889
    n_HI : float
        Density of neutral hydrogen
        (cm^-3)
    xi : float
    Ratio of neutral to singly ionized
    hydrogen density
    
    Returns
    -------
    flux_ratio : float
    '''
    emis_lines = np.sort(np.concatenate((balmer_lines, helium_lines)))[1:] # 1: to remove the duplicate HeI+H8 3890 line

    species = []
    wavelength = []
    flux_ratio = []
    
    dens = 10**log_dens
    xi = 10**log_xi

    collisional_to_recomb_Hbeta = hydrogen_collision_to_recomb(xi, balmer_lines[1], temp)
    AHbeta_AV = reddening_coefficient(balmer_lines[1])

    for w in range(len(waves)):
        # Determine if working with hydrogen or helium line; within 3 Angstroms is arbitrary but should cover difference in vacuum vs air wavelength
        nearest_wave = emis_lines[np.where(np.abs(emis_lines - waves[w]) < 3)[0]][0]
        
        # Any Balmer line besides the blended HeI+H8 line (H8 at 3890.166)
        if nearest_wave in balmer_lines and nearest_wave != 3890.166:
            line_species = 'hydrogen'
            
            emissivity_ratio = hydrogen_emissivity(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_ratio = hydrogen_collision_to_recomb(xi, waves[w], temp)            
            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.           

            flux = emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_H_at_wave)/(EWs[w]) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
        # Any HeI line besides the blended HeI+H8 line (HeI at 3890.151)
        elif nearest_wave in helium_lines and nearest_wave != 3890.151:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity_2007_BSS(waves[w], temp, dens) # Includes correction for HeI collisional to recombination
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            
            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.  

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EWs[w] + a_He_at_wave)/(EWs[w]) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
        
        # The blended HeI+H8 line
        elif nearest_wave == 3890.151 or nearest_wave == 3890.166:
            # Calculate fractional contribution of HeI and H8 to the blended line
            frac_of_he = 0.5
            EW_HeI = frac_of_he * EWs[w]
            EW_H8 = (1-frac_of_he) * EWs[w]

            reddening_function = (reddening_coefficient(waves[w]) / AHbeta_AV ) - 1.  

            # HeI 3890.151 contribution:
            line_species = 'helium'
            
            emissivity_ratio = helium_emissivity_2007_BSS(waves[w], temp, dens)
            a_He_at_wave = stellar_absorption(waves[w], a_He, ion=line_species)            
            optical_depth_at_wave = optical_depth_function(waves[w], temp, dens, tau_He)            

            flux = y_plus * emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_HeI + a_He_at_wave)/(EW_HeI) ) ) * \
                    optical_depth_at_wave * ( 1 / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)
                    
            # H8 contribution:
            line_species = 'hydrogen'

            emissivity_ratio = hydrogen_emissivity(waves[w], temp, dens)
            a_H_at_wave = stellar_absorption(waves[w], a_H, ion=line_species)            
            collisional_to_recomb_factor = np.exp(( -13.6 * ((1/5**2)-(1/8**2)) ) / (8.6173303e-5 * temp)) # scale factor for C/R(Hg) to C/R(H8)
            collisional_to_recomb_ratio = collisional_to_recomb_factor * hydrogen_collision_to_recomb(xi, 4341.684, temp) # Calculate C/R(Hg) and multiply by above scale factor

            flux += emissivity_ratio * ( ( (EW_Hb + a_H)/(EW_Hb) ) / ( (EW_H8 + a_H_at_wave)/(EW_H8) ) ) * \
                    ( (1 + collisional_to_recomb_ratio) / (1 + collisional_to_recomb_Hbeta) ) * \
                    10**-(reddening_function * c_Hb)

        flux_ratio.append(flux)
        species.append(line_species)
        wavelength.append(waves[w])
        
    flux_table = Table([wavelength, species, np.array(flux_ratio), EWs], names=('Wavelength', 'Species', 'Flux Ratio', 'EW'))
    flux_table.write(path+'/'+filename, format='ascii', overwrite=True)
    
    return
