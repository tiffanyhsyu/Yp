# This script compares AOS2015 parameter values to our own to see which parameters are most different.
import numpy as np
import pdb
from astropy.table import Table
from matplotlib import pyplot as plt


def plotpar(num, pcen, ppos, pneg, shft=0.1, xplt=None):
    plt.subplot(4, 2, num)
    if xplt is None:
        xplt = np.arange(pcen.size)
    if shft > 0.0:
        color = 'b'
    else:
        color = 'r'
    plt.plot(xplt - shft, pcen, color+'x')
    plt.errorbar(xplt - shft, pcen, yerr=np.vstack((pneg, ppos)), linestyle='None', color=color)
    return

# Load in AOS2015 output results
aos2015 = Table.read('test_data/aos2015_mcmcoutput', format='ascii', delimiter=' ')

objnames = aos2015['Object']
yp, yp_p, yp_e = aos2015['y+'], aos2015['y+_p'], -aos2015['y+_m']
dens, dens_p, dens_e = aos2015['dens'], aos2015['dens_p'], -aos2015['dens_m']
aHe, aHe_p, aHe_e = aos2015['aHe'], aos2015['aHe_p'], -aos2015['aHe_m']
tauHe, tauHe_p, tauHe_e = aos2015['tauHe'], aos2015['tauHe_p'], -aos2015['tauHe_m']
temp, temp_p, temp_e = aos2015['temp'], aos2015['temp_p'], -aos2015['temp_m']
cHb, cHb_p, cHb_e = aos2015['cHb'], aos2015['cHb_p'], -aos2015['cHb_m']
aH, aH_p, aH_e = aos2015['aH'], aos2015['aH_p'], -aos2015['aH_m']
xi, xi_p, xi_e = aos2015['xi'], aos2015['xi_p'], -aos2015['xi_m']

aos_yp, aos_dens = yp.copy(),  dens.copy()

# Plot all of the data
plotpar(1, yp, yp_p, yp_e)
plotpar(2, dens, dens_p, dens_e)
plotpar(3, aHe, aHe_p, aHe_e)
plotpar(4, tauHe, tauHe_p, tauHe_e)
plotpar(5, temp, temp_p, temp_e)
plotpar(6, cHb, cHb_p, cHb_e)
plotpar(7, aH, aH_p, aH_e)
plotpar(8, xi, xi_p, xi_e)

# Load our output results
ourres = Table.read('all_output', format='ascii', delimiter=' ')

our_objnames = ourres['Object']

xplt = []
for ii in range(our_objnames.size):
    xplt.append(np.where(our_objnames[ii] == objnames.data)[0][0])
xplt = np.array(xplt)

yp, yp_p, yp_e = ourres['y+'], ourres['y+_p'], ourres['y+_m']
dens, dens_p, dens_e = ourres['dens'], ourres['dens_p'], ourres['dens_m']
aHe, aHe_p, aHe_e = ourres['aHe'], ourres['aHe_p'], ourres['aHe_m']
tauHe, tauHe_p, tauHe_e = ourres['tauHe'], ourres['tauHe_p'], ourres['tauHe_m']
temp, temp_p, temp_e = ourres['temp'], ourres['temp_p'], ourres['temp_m']
cHb, cHb_p, cHb_e = ourres['cHb'], ourres['cHb_p'], ourres['cHb_m']
aH, aH_p, aH_e = ourres['aH'], ourres['aH_p'], ourres['aH_m']
xi, xi_p, xi_e = ourres['xi'], ourres['xi_p'], ourres['xi_m']

# Plot all of the data
shft = -0.1
plotpar(1, yp, yp_p, yp_e, shft=shft, xplt=xplt)
plotpar(2, dens, dens_p, dens_e, shft=shft, xplt=xplt)
plotpar(3, aHe, aHe_p, aHe_e, shft=shft, xplt=xplt)
plotpar(4, tauHe, tauHe_p, tauHe_e, shft=shft, xplt=xplt)
plotpar(5, temp, temp_p, temp_e, shft=shft, xplt=xplt)
plotpar(6, cHb, cHb_p, cHb_e, shft=shft, xplt=xplt)
plotpar(7, aH, aH_p, aH_e, shft=shft, xplt=xplt)
plotpar(8, xi, xi_p, xi_e, shft=shft, xplt=xplt)

plt.show()

plt.clf()
plt.plot(aos_dens[(xplt,)], (aos_yp[(xplt,)]-yp)/yp, 'bx')
plt.show()