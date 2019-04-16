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
    emissivity_ratio = np.zeros(dens.size)
    for ii in range(helium_lines.size):
        for jj in range(dens.size):
            emissivity_ratio[jj] = mfr.helium_emissivity_PFSD2012(helium_lines[ii], temp, dens[jj], deg='cubic')
        plt.plot(dens, emissivity_ratio/emissivity_ratio[0], 'k--')
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


def test_HeH_emissivity():
    # Figure 1 from AOS2013
    lines = helium_lines[::-1][:-1]
    dens = 100.0
    temp = np.linspace(1.2, 2.0, 100)
    emissivity_ratio = np.zeros(temp.size)
    for ii in range(lines.size):
        plt.subplot(3, 2, ii+1)
        for jj in range(temp.size):
            emissivity_ratio[jj] = mfr.helium_emissivity_PFSD2012(lines[ii], temp[jj] * 1.0E4, dens)
        plt.plot(temp, emissivity_ratio, 'k-')
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


def test_H_Hbeta_emissivities():
    # Calculate the storey result
    ha_RBS = np.zeros((21, 6))
    hb_RBS = np.zeros((21, 6))
    hg_RBS = np.zeros((21, 6))
    hd_RBS = np.zeros((21, 6))
    h8_RBS = np.zeros((21, 6))
    temp = np.linspace(1.0, 2.0, 100) * 1.0E4
    logdens = 0.0
    dens = 10.0**logdens
    for t in range(len(np.arange(5000, 26000, 1000))):
        ha_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 3)[0],
                                                                        np.where(hydrogen_emis['Nl'] == 2)[0],
                                                                        np.where(hydrogen_emis['T'] ==
                                                                                 np.arange(5000, 26000, 1000)[t])))]
        hb_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 4)[0],
                                                                        np.where(hydrogen_emis['Nl'] == 2)[0],
                                                                        np.where(hydrogen_emis['T'] ==
                                                                                 np.arange(5000, 26000, 1000)[t])))]
        hg_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 5)[0],
                                                                        np.where(hydrogen_emis['Nl'] == 2)[0],
                                                                        np.where(hydrogen_emis['T'] ==
                                                                                 np.arange(5000, 26000, 1000)[t])))]
        hd_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 6)[0],
                                                                        np.where(hydrogen_emis['Nl'] == 2)[0],
                                                                        np.where(hydrogen_emis['T'] ==
                                                                                 np.arange(5000, 26000, 1000)[t])))]
        h8_RBS[t] = hydrogen_emis['emissivity'][reduce(np.intersect1d, (np.where(hydrogen_emis['Nu'] == 8)[0],
                                                                        np.where(hydrogen_emis['Nl'] == 2)[0],
                                                                        np.where(hydrogen_emis['T'] ==
                                                                                 np.arange(5000, 26000, 1000)[t])))]
    # Linear
    S2018_ha_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), 10.0 ** np.arange(0, 6), ha_RBS, kx=1, ky=1)
    S2018_hb_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), 10.0 ** np.arange(0, 6), hb_RBS, kx=1, ky=1)
    S2018_hg_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), 10.0 ** np.arange(0, 6), hg_RBS, kx=1, ky=1)
    S2018_hd_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), 10.0 ** np.arange(0, 6), hd_RBS, kx=1, ky=1)
    S2018_h8_lin = interp.RectBivariateSpline(np.arange(5000, 26000, 1000), 10.0 ** np.arange(0, 6), h8_RBS, kx=1, ky=1)
    # Calculate emissivities
    Ha_emis_storey = S2018_ha_lin(temp, dens).flatten()
    Hb_emis_storey = S2018_hb_lin(temp, dens).flatten()
    Hg_emis_storey = S2018_hg_lin(temp, dens).flatten()
    Hd_emis_storey = S2018_hd_lin(temp, dens).flatten()
    H8_emis_storey = S2018_h8_lin(temp, dens).flatten()
    plt.plot(temp, Ha_emis_storey/Hb_emis_storey, 'c-')
    plt.plot(temp, Hg_emis_storey/Hb_emis_storey, 'r-')
    plt.plot(temp, Hd_emis_storey/Hb_emis_storey, 'g-')
    plt.plot(temp, H8_emis_storey/Hb_emis_storey, 'b-')
    # Now plot the old HS87 emissivities
    lines = hydrogen_lines[1:]  # Ignore Pg
    res = np.zeros(temp.size)
    plt.clf()
    for ii in range(len(lines)):
        for t in range(temp.size):
            res[t] = mfr.hydrogen_emissivity_HS1987(lines[ii], temp[t], dens)
        if False:
            plt.plot(temp, res, 'k--')
        else:
            if ii == 0:
                comp, col = Ha_emis_storey/Hb_emis_storey, 'k'
            if ii == 2:
                comp, col = Hg_emis_storey/Hb_emis_storey, 'r'
            if ii == 3:
                comp, col = Hd_emis_storey/Hb_emis_storey, 'g'
            if ii == 4:
                comp, col = H8_emis_storey/Hb_emis_storey, 'b'
            plt.plot(temp, (comp-res)/res, col+'-')
    plt.show()


def test_xi():
    xi = 1.0e-4
    lines = hydrogen_lines[1:-1]
    temp = np.linspace(1.0, 2.5, 100)
    cols = ['r', 'm', 'g', 'b']
    for ii in range(lines.size):
        CR = mfr.hydrogen_collision_to_recomb(xi, lines[ii], temp * 1.0E4)
        plt.plot(temp, CR, cols[ii]+"-")
    plt.show()


if __name__ == "__main__":
    #test_He_emissivity()
    #test_HeH_emissivity()
    #test_Hbeta_emissivity()
    #test_Pg_Hbeta_emissivity()
    #test_xi()
    #test_od_func()
    test_H_Hbeta_emissivities()
