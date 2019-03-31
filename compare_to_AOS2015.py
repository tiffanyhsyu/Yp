# This script compares AOS2015 parameter values to our own to see which parameters are most different.
import fileinput
import numpy as np
import shutil

from astropy.table import Column
from astropy.table import Table
from matplotlib import pyplot as plt
import corner


def plotpar(num, pcen, ppos, pneg):
    plt.subplot(4, 2, num)
    objplot = np.arange(pcen.size)
    shft = 0.1
    plt.plot(objplot - shft, pcen, 'bx')
    plt.errorbar(objplot - shft, pcen, yerr=np.vstack((-pneg, ppos)), linestyle='None')
    return

# Load in AOS2015 output results
aos2015 = Table.read('test_data/aos2015_mcmcoutput', format='ascii', delimiter=' ')

objnames = aos2015['Object']
yp, yp_p, yp_e = aos2015['y+'], aos2015['y+_p'], aos2015['y+_m']
dens, dens_p, dens_e = aos2015['dens'], aos2015['dens_p'], aos2015['dens_m']
aHe, aHe_p, aHe_e = aos2015['aHe'], aos2015['aHe_p'], aos2015['aHe_m']
tauHe, tauHe_p, tauHe_e = aos2015['tauHe'], aos2015['tauHe_p'], aos2015['tauHe_m']
temp, temp_p, temp_e = aos2015['temp'], aos2015['temp_p'], aos2015['temp_m']
cHb, cHb_p, cHb_e = aos2015['cHb'], aos2015['cHb_p'], aos2015['cHb_m']
aH, aH_p, aH_e = aos2015['aH'], aos2015['aH_p'], aos2015['aH_m']
xi, xi_p, xi_e = aos2015['xi'], aos2015['xi_p'], aos2015['xi_m']

# Plot all of the data
plotpar(1, yp, yp_p, yp_e)
plotpar(2, dens, dens_p, dens_e)
plotpar(3, aHe, aHe_p, aHe_e)
plotpar(4, tauHe, tauHe_p, tauHe_e)
plotpar(5, temp, temp_p, temp_e)
plotpar(6, cHb, cHb_p, cHb_e)
plotpar(7, aH, aH_p, aH_e)
plotpar(8, xi, xi_p, xi_e)
plt.show()




