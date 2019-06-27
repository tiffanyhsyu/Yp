import os
from astropy.table import Table

#########################
# AOS2015 Optical + NIR #
#########################
def load_AOS2015(galaxyname):
    # T_OIII are taken from Column 4, Table 3, Izotov et al. (2014)
    outdict = dict()
    dir = '/test_data/optical+nir/'
    if galaxyname == 'IZw18SE1':
        gname = 'IZw18SE1_nir'
        T_OIII = 19000.0#17500.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS0335-052E1':
        gname = 'SBS0335052E1_nir'
        T_OIII = 20100  #20000.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS0335-052E3':
        gname = 'SBS0335052E3_nir'
        T_OIII = 20800  #20000.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0519+0007':
         gname = 'J0519+0007_nir'
         T_OIII = 20700
         full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS0940+5442':
        gname = 'SBS0940+5442_nir'
        T_OIII = 18700#18600.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Tol65':
        gname = 'Tol65_nir'
        T_OIII = 17200#17100.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1415+437No13':
        gname = 'SBS1415+437No13_nir'
        T_OIII = 17000#16800.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1415+437No2':
        gname = 'SBS1415+437No2_nir'
        T_OIII = 15900#16800.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'CGCG007-025No2':
        gname = 'CGCG007025No2_nir'
        T_OIII = 16500.0  # 16700.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk209':
        gname = 'Mrk209_nir'
        T_OIII = 16100#16400.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1030+583':
        gname = 'SBS1030+583_nir'
        T_OIII = 15500#15400.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk71No1':
        gname = 'Mrk71No1_nir'
        T_OIII = 15600#15700.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1152+579':
        gname = 'SBS1152+579_nir'
        T_OIII = 15400#15100.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk59':
        gname = 'Mrk59_nir'
        T_OIII = 13500#13300.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1135+581':
        gname = 'SBS1135+581_nir'
        T_OIII = 12600#12600.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk450No1':
        gname = 'Mrk450No1_nir'
        T_OIII = 11700#11600.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Test':
        #gname = 'test_output_flux_nir'
        gname = 'test_output_flux'
        T_OIII = 18000
        full_tbl = Table.read(os.getcwd() + '/' + gname, format='ascii', delimiter=' ')
    else:
        print('Galaxy not known: {0:s}'.format(galaxyname))
        return None
    outdict['full_tbl'] = full_tbl
    outdict['T_OIII'] = T_OIII
    return outdict

###################
# AOS2012 Optical #
###################
def load_AOS2012(galaxyname):
    # T_OIII are taken from Column 4, Table 3, Izotov et al. (2014)
    outdict = dict()
    dir = '/test_data/optical/'
    if galaxyname == 'IZw18SE1':
        gname = 'IZw18SE1'
        T_OIII = 19000.0#17500.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS0335-052E1':
        gname = 'SBS0335052E1'
        T_OIII = 20100  #20000.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS0335-052E3':
        gname = 'SBS0335052E3'
        T_OIII = 20800  #20000.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0519+0007':
         gname = 'J0519+0007'
         T_OIII = 20700
         full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS0940+5442':
        gname = 'SBS0940+5442'
        T_OIII = 18700#18600.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Tol65':
        gname = 'Tol65'
        T_OIII = 17200#17100.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1415+437No13':
        gname = 'SBS1415+437No13'
        T_OIII = 17000#16800.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1415+437No2':
        gname = 'SBS1415+437No2'
        T_OIII = 15900#16800.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'CGCG007-025No2':
        gname = 'CGCG007025No2'
        T_OIII = 16500.0  # 16700.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk209':
        gname = 'Mrk209'
        T_OIII = 16100#16400.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1030+583':
        gname = 'SBS1030+583'
        T_OIII = 15500#15400.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk71No1':
        gname = 'Mrk71No1'
        T_OIII = 15600#15700.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1152+579':
        gname = 'SBS1152+579'
        T_OIII = 15400#15100.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk59':
        gname = 'Mrk59'
        T_OIII = 13500#13300.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1135+581':
        gname = 'SBS1135+581'
        T_OIII = 12600#12600.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk450No1':
        gname = 'Mrk450No1'
        T_OIII = 11700#11600.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Test':
        gname = 'test_output_flux_opt'
        T_OIII = 18000
        full_tbl = Table.read(os.getcwd() + '/' + gname, format='ascii', delimiter=' ')
    else:
        print('Galaxy not known: {0:s}'.format(galaxyname))
        return None
    outdict['full_tbl'] = full_tbl
    outdict['T_OIII'] = T_OIII
    return outdict

###########################
# Erik's Synthetic Fluxes #
###########################
def load_synthetic(galaxyname):
    outdict = dict()
    dir = '/test_data/synthetic/'
    if galaxyname == 'synthetic1':
        gname = 'synthetic1_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'synthetic2':
        gname = 'synthetic2_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'synthetic3':
        gname = 'synthetic3_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'synthetic4':
        gname = 'synthetic4_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'synthetic5':
        gname = 'synthetic5_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'synthetic6':
        gname = 'synthetic6_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'synthetic7':
        gname = 'synthetic7_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'synthetic8':
        gname = 'synthetic8_erik'
        T_OIII = 16500.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    else:
        print('Galaxy not known: {0:s}'.format(galaxyname))
        return None
    outdict['full_tbl'] = full_tbl
    outdict['T_OIII'] = T_OIII
    return outdict

##########################
# Our Optical+NIR Sample #
##########################

def load_ours(galaxyname):
    outdict = dict()
    ##############################
    # Relative to Ha on red side #
    ##############################
    dir = '/test_data/no_HaHb/'
#    if galaxyname == 'LeoP':
#        gname = 'LeoP'
#        T_OIII = 17350.
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    if galaxyname == 'KJ2':
        gname = 'KJ2'
        T_OIII = 17413.404518841133
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'KJ5':
        gname = 'KJ5'
        T_OIII = 11501.089982360061
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'KJ29':
        gname = 'KJ29'
        T_OIII = 14212.458597864896
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'KJ5B':
        gname = 'KJ5B'
        T_OIII = 13979.130375718989
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'KJ97':
        gname = 'KJ97'
        T_OIII = 11822.063460125153
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0000p3052A':
        gname = 'J0000p3052A'
        T_OIII = 15034.028576533961
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0000p3052B':
        gname = 'J0000p3052B'
        T_OIII = 15110.068679551327
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0018p2345':
        gname = 'J0018p2345'
        T_OIII = 16740.044632787234
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0118p3512':
        gname = 'J0118p3512'
        T_OIII = 15414.42499005176
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0140p2951':
        gname = 'J0140p2951'
        T_OIII = 12154.766132309507
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0201p0919':
        gname = 'J0201p0919'
        T_OIII = 14672.599263221606
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0220p2044A':
        gname = 'J0220p2044A'
        T_OIII = 15798.807188110046
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0220p2044B':
        gname = 'J0220p2044B'
        T_OIII = 17523.586639214853
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0452m0541':
        gname = 'J0452m0541'
        T_OIII = 15399.829000035199
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0743p4807':
        gname = 'J0743p4807'
        T_OIII = 9390.403459203068
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0757p4750':
        gname = 'J0757p4750'
        T_OIII = 15846.837228485041
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'J0812p4836':
#        gname = 'J0812p4836'
#        T_OIII = 17864.510725965007
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'J0834p5905':
#        gname = 'J0834p5905'
#        T_OIII = 25900.258903826118
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J0943p3326':
        gname = 'J0943p3326'
        T_OIII = 16412.688803906996
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1044p6306':
        gname = 'J1044p6306'
        T_OIII = 18700.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1204p5259':
        gname = 'J1204p5259'
        T_OIII = 13295.120587868281
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1214p1245':
        gname = 'J1214p1245'
        T_OIII = 14896.286859605134
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1322p5425':
        gname = 'J1322p5425'
        T_OIII = 17158.63589425323
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1414m0208':
        gname = 'J1414m0208'
        T_OIII = 14675.038910218094
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1425p4441':
        gname = 'J1425p4441'
        T_OIII = 14986.196617398991
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1655p6337':
        gname = 'J1655p6337'
        T_OIII = 16513.50413647784
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1705p3527':
        gname = 'J1705p3527'
        T_OIII = 15428.266502581113
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1732p4452':
        gname = 'J1732p4452'
        T_OIII = 15131.086884967399
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J1757p6454':
        gname = 'J1757p6454'
        T_OIII = 14424.061263892125
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J2030m1343':
        gname = 'J2030m1343'
        T_OIII = 13817.727654937002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J2213p1722':
        gname = 'J2213p1722'
        T_OIII = 15328.109575371638
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J2230m0531':
        gname = 'J2230m0531'
        T_OIII = 14795.751597494424
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J2319p1616':
        gname = 'J2319p1616'
        T_OIII = 10576.43370005504
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'J2339p3230':
        gname = 'J2339p3230'
        T_OIII = 13931.435045664843
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    else:
        print('Galaxy not known: {0:s}'.format(galaxyname))
        return None
    outdict['full_tbl'] = full_tbl
    outdict['T_OIII'] = T_OIII
    return outdict