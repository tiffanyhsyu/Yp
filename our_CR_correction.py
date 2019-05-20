import numpy as np
import model_flux_ratio as mfr

## Effective collision strengths, Gamma, from Anderson et al. 2002: https://iopscience.iop.org/article/10.1088/0953-4075/35/6/701
# a correction of Anderson et al. 2000: https://iopscience.iop.org/article/10.1088/0953-4075/33/6/311/pdf
## Branching ratios, BR, from Omidvar 1983: https://www.sciencedirect.com/science/article/pii/0092640X83900116
## Recombination rates from Hummer & Storey 1987

kB = 8.61733e-5

def eV_to_K(eV):
    temp = eV * (1/kB)
    return temp

# Anderson et al. 2002 gives Gamma values at given electron temperatures...
electron_temp = np.array([0.5, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0])

# Working at T=11604.522 to avoid Gamma interpolations for now
T = eV_to_K(electron_temp[1])

# Effective collision strength, G_ij, from Table 1 of Anderson et al. 2002
# Halpha: 1-->3s, 1-->3p, 1-->3d, 1-->4s, 1-->4p, 1-->4d, 1-->4f
# In A2000's definition of i and j, these correspond to: 1-->4, 1-->5, 1-->6, 1-->7, 1-->8, 1-->9, 1-->10
G_Ha = np.array([6.96e-2, 1.26e-1, 6.58e-2, 2.55e-2, 4.79e-2, 3.19e-2, 1.14e-2])
# Hbeta: 1-->4s, 1-->4p, 1-->4d, 1-->4f, 1-->5s, 1-->5p, 1-->5d, 1-->5f 1-->5g
# 1-->7, 1-->8, 1-->9, 1-->10, 1-->11, 1-->12, 1-->13, 1-->14, 1-->15
G_Hb = np.array([2.55e-2, 4.79e-2, 3.19e-2, 1.14e-2, 1.72e-2, 3.15e-2, 2.22e-2, 9.14e-3, 4.03e-3])
# Hgamma: 1-->5s, 1-->5p, 1-->5d, 1-->5f 1-->5g
# 1-->11, 1-->12, 1-->13, 1-->14, 1-->15
G_Hg = np.array([1.72e-2, 3.15e-2, 2.22e-2, 9.14e-3, 4.03e-3])

# Branching Ratios from Table 2 of Omidvar 1983
BR_Ha = np.array([1.0, 1.0, 1.0, 4.16e-1, 4.20e-2, 2.54e-1, 1.0]) # last 4 are from 40,41,42,43 --> 3 
BR_Hb = np.array([1.0, 1.0, 1.0, 1.0, 2.27e-1, 2.20e-2, 1.07e-1, 3.63e-1, 1.0]) # last 5 are from 50,51,52,53,54 --> 4
BR_Hg = np.array([1.0, 1.0, 1.0, 1.0, 1.0])

def calculate_CR(gamma, branching_ratio, temp, eta, line):
    kB = 8.61733e-5
    
    # Recombination rate and scaling factors from H&S 1987
    # Taking T=10000K, n_e=100 for now to be close to T=11604.522, on p.59
    #### Eventually need to interpolate based on T, n_e
    recomb_42 = 3.022e-14

    if line == 'Ha':
        scale = 2.86
        i = np.array([3, 3, 3, 4, 4, 4, 4])
    elif line == 'Hb':
        scale = 1.00
        i = np.array([4, 4, 4, 4, 5, 5, 5, 5, 5])
    elif line == 'Hg':
        scale = 4.68e-1
        i = np.array([5, 5, 5, 5, 5])
    else:
        print ('Not ready for this hydrogen transition')
        
    K = 4.004e-8 * np.sqrt( 1/(kB*temp) ) * np.exp(-13.6*(1-(1/i**2))/(kB*temp)) * gamma
    numerator = K * branching_ratio
    CR = eta * np.sum(numerator) / (recomb_42*scale)
    
    return CR

eta = 1e-2
print ('For eta =', eta)
print ('Our C/R values:')
print ('C/R(Halpha): ', calculate_CR(G_Ha, BR_Ha, T, eta, 'Ha'))
print ('C/R(Hbeta): ', calculate_CR(G_Hb, BR_Hb, T, eta, 'Hb'))
print ('C/R(Hgamma): ', calculate_CR(G_Hg, BR_Hg, T, eta, 'Hg'))

print ('Compared to Erik:')

hydrogen_lines = np.array([10941.082, 6564.612, 4862.721, 4341.684, 4102.891, 3890.166])

print ('C/R(Halpha): ', mfr.hydrogen_collision_to_recomb(eta, hydrogen_lines[1], T))
print ('C/R(Hbeta): ', mfr.hydrogen_collision_to_recomb(eta, hydrogen_lines[2], T))
print ('C/R(Hgamma): ', mfr.hydrogen_collision_to_recomb(eta, hydrogen_lines[3], T))