import pdb
import numpy as np
import model_flux_ratio as mfr

## Effective collision strengths, Gamma, from Anderson et al. 2002: https://iopscience.iop.org/article/10.1088/0953-4075/35/6/701
# a correction of Anderson et al. 2000: https://iopscience.iop.org/article/10.1088/0953-4075/33/6/311/pdf
## Branching ratios, BR, from Omidvar 1983: https://www.sciencedirect.com/science/article/pii/0092640X83900116
## Recombination rates from Hummer & Storey 1987

kB = 8.61733e-5

def eV_to_K(eV):
    temp = eV * (1 / kB)
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
BR_Ha = np.array([1.0, 1.0, 1.0, 4.16e-1, 4.20e-2, 2.54e-1, 1.0])  # last 4 are from 40,41,42,43 --> 3
BR_Hb = np.array([1.0, 1.0, 1.0, 1.0, 2.27e-1, 2.20e-2, 1.07e-1, 3.63e-1, 1.0])  # last 5 are from 50,51,52,53,54 --> 4
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

    print("Erik", np.sum(gamma*branching_ratio))
    K = 4.004e-8 * np.sqrt(1 / (kB * temp)) * np.exp(-13.6 * (1 - (1 / i ** 2)) / (kB * temp)) * gamma
    numerator = K * branching_ratio
    CR = eta * np.sum(numerator) / (recomb_42 * scale)

    return CR


def colrate_raga(eta, T, line):
    kB = 8.61733e-5
    if line == 'Ha':
        i = 3
        acoeffs = np.array([0.2500, 0.2461, 0.3297, 0.3892, -0.0928, 0.0071])
        bcoeffs = np.array([-13.3377, -0.7161, -0.1435, -0.0386, 0.0077])
    elif line == 'Hb':
        i = 4
        acoeffs = np.array([0.1125, 0.1370, -0.1152, 0.1209, -0.0276, 0.0020])
        bcoeffs = np.array([-13.5225, -0.7928, -0.1749, -0.0412, 0.0154])
    elif line == 'Hg':
        i = 5
        acoeffs = np.array([0.0773, 0.0678, -0.0945, 0.0796, -0.0177, 0.0013])
        bcoeffs = np.array([-13.6820, -0.8629, -0.1957, -0.0375, 0.0199])
    else:
        print("Line not ready yet")
    ediff = -13.6 * (1 - (1 / i ** 2))
    tval = np.log10(T/1.0E4)
    omega1k = np.zeros(tval.size)
    for aa in range(acoeffs.size):
        omega1k += acoeffs[aa] * tval**aa
    print("Raga", omega1k)
    # Note 4.004E-8/np.sqrt(kB) = 0.5 * 8.629E-6  (i.e. Erik = Raga for the coefficient)
    q1k = 0.5 * 8.629E-6 * omega1k * np.exp(ediff/(kB*T)) / np.sqrt(T)
    alphak = np.zeros(tval.size)
    for aa in range(bcoeffs.size):
        alphak += bcoeffs[aa] * tval**aa
    alphak = 10.0**(alphak)
    CR = eta * q1k/alphak
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

plotit = True
if plotit:
    T = np.linspace(1.0E4, 2.5E4, 1000)
    from matplotlib import pyplot as plt
    plt.plot(T, mfr.hydrogen_collision_to_recomb(eta, hydrogen_lines[1], T), 'r-')
    plt.plot(T, mfr.hydrogen_collision_to_recomb(eta, hydrogen_lines[2], T), 'm-')
    plt.plot(T, mfr.hydrogen_collision_to_recomb(eta, hydrogen_lines[3], T), 'g-')
    plt.plot(T, colrate_raga(eta, T, 'Ha'), 'r-', ls='--')
    plt.plot(T, colrate_raga(eta, T, 'Hb'), 'm-', ls='--')
    plt.plot(T, colrate_raga(eta, T, 'Hg'), 'g-', ls='--')
    plt.show()

print ('Compared to Raga et al. (2015):')

hydrogen_lines = np.array([10941.082, 6564.612, 4862.721, 4341.684, 4102.891, 3890.166])

print ('C/R(Halpha): ', colrate_raga(eta, T, 'Ha'))
print ('C/R(Hbeta): ', colrate_raga(eta, T, 'Hb'))
print ('C/R(Hgamma): ', colrate_raga(eta, T, 'Hg'))