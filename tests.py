"""
Test various aspects of the code
"""
import os
import pdb
import numpy as np
import model_flux_ratio as mfr
from matplotlib import pyplot as plt
from astropy.table import Table
import scipy.interpolate as interp
from functools import reduce

# Load in tables we'll need
path = os.getcwd()

hydrogen_emis = Table.read(path+'/tables/hydrogen_emissivity_S2018', format='ascii', delimiter='\t')

hydrogen_lines = np.array([10941.082, 6564.612, 4862.721, 4341.684, 4102.891, 3890.166])  # Pa-g, Ha, Hb, Hg, Hd, H8
helium_lines = np.array([10833.306, 7067.198, 6679.994, 5877.299, 4472.755, 4027.328, 3890.151])
hel_colour = ['r', 'm', 'y', 'g', 'c', 'b', 'k']
emis_lines = np.sort(np.concatenate((hydrogen_lines, helium_lines)))[1:-1]


def test_He_emissivity():
    # Figure 1 from AOS2015
    temp = 18000.0
    dens = np.linspace(0.0, 300.0, 100)
    emissivity_ratio = np.zeros(dens.size)
    plt.subplot(121)
    for ii in range(helium_lines.size):
        for jj in range(dens.size):
            emissivity_ratio[jj] = mfr.helium_emissivity_PFSD2012(helium_lines[ii], temp, dens[jj])
        plt.plot(dens, emissivity_ratio/emissivity_ratio[0], hel_colour[ii]+'-')
    # Figure 2 from AOS2015
    dens = 100.0
    temp = np.linspace(1.2, 2.0, 100)
    emissivity_ratio = np.zeros(temp.size)
    plt.subplot(122)
    for ii in range(helium_lines.size):
        for jj in range(temp.size):
            emissivity_ratio[jj] = mfr.helium_emissivity_PFSD2012(helium_lines[ii], temp[jj]*1.0E4, dens)
        plt.plot(temp, emissivity_ratio/emissivity_ratio[0], hel_colour[ii]+'-')
    plt.show()
    return


def test_Hbeta_emissivity():
    # Calculate the porter result
    temp = np.linspace(1.0, 2.0, 100)*1.0E4
    Hbeta_emis_porter = (-2.6584e5 - (1420.9 * (np.log(temp) ** 2.)) + (35546 * np.log(temp)) + (6.5669e5 / np.log(temp))) * (1 / temp) * 1e-25
    # Calculate the storey result
    hb_RBS = np.zeros((21, 6))
    for t in range(len(np.arange(5000, 26000, 1000))):
        hb_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 4)[0], \
                                                                        np.where(hydrogen_emis['Nl'] == 2)[0], \
                                                                        np.where(hydrogen_emis['T'] ==
                                                                                 np.arange(5000, 26000, 1000)[t])))]
    S2018_hb_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0, 6), hb_RBS, kx=1, ky=1)
    Hbeta_emis_storey0 = S2018_hb_lin(temp, 0.0).flatten()
    Hbeta_emis_storey1 = S2018_hb_lin(temp, 1.0).flatten()
    Hbeta_emis_storey2 = S2018_hb_lin(temp, 2.0).flatten()
    Hbeta_emis_storey3 = S2018_hb_lin(temp, 3.0).flatten()
    tmp_pnts = np.arange(5000, 26000, 1000)
    Hbeta_storey_pnts = S2018_hb_lin(tmp_pnts, np.zeros(tmp_pnts.size))
    # Plot the comparisons
    plt.plot(temp, Hbeta_emis_porter, 'k-')
    plt.plot(temp, Hbeta_emis_storey0, 'r--')
    plt.plot(temp, Hbeta_emis_storey1, 'g--')
    plt.plot(temp, Hbeta_emis_storey2, 'b--')
    plt.plot(tmp_pnts, Hbeta_storey_pnts, 'ro')
    plt.show()
    plt.clf()
    plt.plot(temp, (Hbeta_emis_porter-Hbeta_emis_storey0)/Hbeta_emis_storey0, 'k-')
    plt.show()
    plt.clf()
    plt.plot(temp, (Hbeta_emis_storey1-Hbeta_emis_storey0)/Hbeta_emis_storey0, 'k-')
    plt.plot(temp, (Hbeta_emis_storey2-Hbeta_emis_storey0)/Hbeta_emis_storey0, 'r-')
    plt.plot(temp, (Hbeta_emis_storey3-Hbeta_emis_storey0)/Hbeta_emis_storey0, 'g-')
    plt.show()


def test_Pg_Hbeta_emissivity():
    # Calculate the storey result
    hb_RBS = np.zeros((21, 6))
    pg_RBS = np.zeros((21, 6))
    temp = np.linspace(1.0, 2.0, 100) * 1.0E4
    for t in range(len(np.arange(5000, 26000, 1000))):
        hb_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 4)[0], \
                                                                        np.where(hydrogen_emis['Nl'] == 2)[0], \
                                                                        np.where(hydrogen_emis['T'] ==
                                                                                 np.arange(5000, 26000, 1000)[t])))]
        pg_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 6)[0], \
                                                              np.where(hydrogen_emis['Nl'] == 3)[0], \
                                                              np.where(hydrogen_emis['T'] == np.arange(5000, 26000, 1000)[t])))]
    S2018_hb_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0, 6), hb_RBS, kx=1, ky=1)
    S2018_pg_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), np.arange(0, 6), pg_RBS, kx=1, ky=1)
    Hbeta_emis_storey = S2018_hb_lin(temp, 0.0).flatten()
    Pg_emis_storey = S2018_pg_lin(temp, 0.0).flatten()
    plt.plot(temp, Pg_emis_storey/Hbeta_emis_storey, 'k-')
    plt.show()


if __name__ == "__main__":
    #test_He_emissivity()
    test_Hbeta_emissivity()
    #test_Pg_Hbeta_emissivity()

