# This script compares AOS2015 parameter values to our own to see which parameters are most different.
import numpy as np
import pdb
from astropy.table import Table
from matplotlib import pyplot as plt
from mcmc_all import MCMCgal


def plotpar(num, xplt, pcen, ppos, pneg, color='k'):
    plt.subplot(4, 2, num)
    plt.plot(xplt, pcen, color+'x')
    plt.errorbar(xplt, pcen, yerr=np.vstack((pneg, ppos)), linestyle='None', color=color)
    return

names = ["IZw18SE1", "SBS0335-052E1", "SBS0335-052E3", "J0519+0007", "SBS0940+5442", "Tol65", "SBS1415+437No13",
         "SBS1415+437No2", "CGCG007-025No2", "Mrk209", "SBS1030+583", "Mrk71No1", "SBS1152+579", "Mrk59",
         "SBS1135+581", "Mrk450No1"]

# galname = "CGCG007-025No2"
# galname = "Mrk450No1"

# Load in AOS2015 output results
aos2015 = Table.read('test_data/aos2015_mcmcoutput', format='ascii', delimiter=' ')

objnames = aos2015['Object']
yp, yp_p, yp_m = aos2015['y+'], aos2015['y+_p'], -aos2015['y+_m']
dens, dens_p, dens_m = aos2015['dens'], aos2015['dens_p'], -aos2015['dens_m']
aHe, aHe_p, aHe_m = aos2015['aHe'], aos2015['aHe_p'], -aos2015['aHe_m']
tauHe, tauHe_p, tauHe_m = aos2015['tauHe'], aos2015['tauHe_p'], -aos2015['tauHe_m']
temp, temp_p, temp_m = aos2015['temp'], aos2015['temp_p'], -aos2015['temp_m']
cHb, cHb_p, cHb_m = aos2015['cHb'], aos2015['cHb_p'], -aos2015['cHb_m']
aH, aH_p, aH_m = aos2015['aH'], aos2015['aH_p'], -aos2015['aH_m']
xi, xi_p, xi_m = aos2015['xi'], aos2015['xi_p'], -aos2015['xi_m']

for ii in range(len(names)):
    idx = np.where(objnames == names[ii])[0][0]
    print(names[ii], xi[idx])

ncols = 2
cntx, cnty = 0, 0
fig, ax = plt.subplots(1+len(names)//ncols - 3, ncols)
for gg, galname in enumerate(names):
    idx = np.where(objnames == galname)[0][0]
    params = (yp[idx], temp[idx], np.log10(dens[idx]), cHb[idx], aH[idx], aHe[idx], tauHe[idx], -4.0,)
    try:
        mcmc_inst = MCMCgal(galname, run_mcmc=False)
    except:
        continue
    model_flux = mcmc_inst.get_model(params)

    # Plot the emission lines
    xplt = np.arange(mcmc_inst._y.size)
    ax[cntx, cnty].plot(xplt, (mcmc_inst._y-model_flux)/mcmc_inst._y, 'bx')
#    ax[gg].plot(xplt, model_flux, 'rx')
    ax[cntx, cnty].set_xticks(xplt)
    #ax[1].plot(xplt, (mcmc_inst._y-model_flux)/mcmc_inst._y, 'bx')
    #ax[1].set_xticks(xplt)
    cnty += 1
    if cnty == 2:
        cnty = 0
        cntx += 1
#mcmc_inst._y_error

# Draw the canvas and update the label names
fig.canvas.draw()
labels = [item.get_text() for item in ax[0,0].get_xticklabels()]
for ii in range(len(mcmc_inst._y_names)):
    labels[ii] = mcmc_inst._y_names[ii]

cntx, cnty = 0, 0
for gg in range(len(names)):
    try:
        ax[cntx, cnty].set_xticklabels(labels)
    except:
        break
    cnty += 1
    if cnty == 2:
        cnty = 0
        cntx += 1

plt.show()
plt.clf()

wrongidx = [3, 7]  # These indices are both wrong, and in the same direction
oppidx = [5, 8]  # These indices have He I 4027 and 4472 wrong (but in the opposite way)
cntr = 0
for gg, galname in enumerate(names):
    fail = False
    try:
        mcmc_inst = MCMCgal(galname, run_mcmc=False)
    except:
        fail = True
        continue
    col = 'k'
    if cntr in wrongidx:
        col = 'b'
    elif cntr == oppidx[0]:
        col = 'g'
    elif cntr == oppidx[1]:
        col = 'r'
    # Plot the data
    idx = np.where(objnames == galname)[0][0]
    plotpar(1, yp[idx], yp[idx], yp_p[idx], yp_m[idx], color=col)
    plotpar(2, yp[idx], dens[idx], dens_p[idx], dens_m[idx], color=col)
    plotpar(3, yp[idx], aHe[idx], aHe_p[idx], aHe_m[idx], color=col)
    plotpar(4, yp[idx], tauHe[idx], tauHe_p[idx], tauHe_m[idx], color=col)
    plotpar(5, yp[idx], temp[idx], temp_p[idx], temp_m[idx], color=col)
    plotpar(6, yp[idx], cHb[idx], cHb_p[idx], cHb_m[idx], color=col)
    plotpar(7, yp[idx], aH[idx], aH_p[idx], aH_m[idx], color=col)
    plotpar(8, yp[idx], xi[idx], xi_p[idx], xi_m[idx], color=col)
    cntr += 1

plt.show()