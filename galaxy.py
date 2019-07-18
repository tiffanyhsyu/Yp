import os
from astropy.table import Table

#########################
# AOS2015 Optical + NIR #
#########################
def load_AOS2015(galaxyname):
    # Used temperatures from Izotov, Thuan, and Stasinska 2007, Table 5; found at ~/Yp/test_data/HeBCD
    # Commented T_OIII are taken from Izotov, Thuan, and Guseva 2014, Table 3, Column 4 'te(OIII)'
    outdict = dict()
    dir = '/test_data/AOS/optical+nir/'
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
    dir = '/test_data/AOS/optical/'
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

#########################################
# All HeBCD galaxies w/optical+NIR data #
#########################################
def load_HeBCD_NIR(galaxyname):
    # Used temperatures from Izotov, Thuan, and Stasinska 2007, Table 5; found at ~/Yp/test_data/HeBCD
    # 12/21 systems are part of AOS2015 and commented out here
    outdict = dict()
    dir = '/test_data/AOS/full_HeBCD/'
#    if galaxyname == 'CGCG007025No2':
#        gname = 'CGCG007025No2'
#        T_OIII = 16500.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    if galaxyname == 'HS0837+4717':
        gname = 'HS0837+4717'
        T_OIII = 19400.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'IZw18SE1':
#        gname = 'IZw18SE1'
#        T_OIII = 19000.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'J0519+0007':
#        gname = 'J0519+0007'
#        T_OIII = 20700.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk162':
        gname = 'Mrk162'
        T_OIII = 11800.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk36':
        gname = 'Mrk36'
        T_OIII = 15100.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk930':
        gname = 'Mrk930'
        T_OIII = 12300.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk1315':
        gname = 'Mrk1315'
        T_OIII = 11000.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'Mrk1329':
        gname = 'Mrk1329'
        T_OIII = 10800.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'Mrk209':
#        gname = 'Mrk209'
#        T_OIII = 16100.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'Mrk450No1':
#        gname = 'Mrk450No1'
#        T_OIII = 11700.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'Mrk71No1':
#        gname = 'Mrk71No1'
#        T_OIII = 15600.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'SBS0940+5442':
#        gname = 'SBS0940+5442'
#        T_OIII = 18700.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'SBS1030+583':
#        gname = 'SBS1030+583'
#        T_OIII = 15500.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'SBS1135+581':
#        gname = 'SBS1135+581'
#        T_OIII = 12600.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'SBS1152+579':
#        gname = 'SBS1152+579'
#        T_OIII = 15400.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1222+614':
        gname = 'SBS1222+614'
        T_OIII = 14000.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'Mrk59':
#        gname = 'Mrk59'
#        T_OIII = 13500.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'SBS1437+370':
        gname = 'SBS1437+370'
        T_OIII = 14200.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
#    elif galaxyname == 'Tol65':
#        gname = 'Tol65'
#        T_OIII = 17200.0
#        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'UM311':
        gname = 'UM311'
        T_OIII = 9700.0
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
    dir = '/test_data/ours/fullspec/'
    if galaxyname == 'LeoP':
        gname = 'LeoP'
        T_OIII = 17350.
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'KJ2':
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
    elif galaxyname == 'J0214m0835':
        gname = 'J0214m0835'
        T_OIII = 17046.256263241805
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

def load_ours_noHaHb(galaxyname):
    outdict = dict()
    ##############################
    # Relative to Ha on red side #
    ##############################
    dir = '/test_data/ours/no_HaHb/'
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
    elif galaxyname == 'J0214m0835':
        gname = 'J0214m0835'
        T_OIII = 17046.256263241805
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

#################
# SDSS Galaxies #
#################
def load_SDSS(galaxyname):
    outdict = dict()
    dir = '/test_data/SDSS/'
    if galaxyname == 'spec-0266-51630-0407':
        gname = 'spec-0266-51630-0407'
        T_OIII = 11194.102994920871
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0284-51943-0408':
        gname = 'spec-0284-51943-0408'
        T_OIII = 12189.413194087507
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0283-51959-0572':
        gname = 'spec-0283-51959-0572'
        T_OIII = 16151.342198899576
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0284-51943-0007':
        gname = 'spec-0284-51943-0007'
        T_OIII = 10628.073500152708
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0278-51900-0392':
        gname = 'spec-0278-51900-0392'
        T_OIII = 16521.48652519121
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-0299-51671-0083':
    #	 gname = 'spec-0299-51671-0083'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0299-51671-0311':
        gname = 'spec-0299-51671-0311'
        T_OIII = 10938.354305551531
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0289-51990-0369':
        gname = 'spec-0289-51990-0369'
        T_OIII = 12639.6307982468
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0285-51930-0154':
        gname = 'spec-0285-51930-0154'
        T_OIII = 11270.454992047864
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0267-51608-0421':
        gname = 'spec-0267-51608-0421'
        T_OIII = 12697.635280234063
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0279-51984-0293':
        gname = 'spec-0279-51984-0293'
        T_OIII = 10194.079724756852
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0279-51984-0520':
        gname = 'spec-0279-51984-0520'
        T_OIII = 17627.59061443254
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0301-51942-0531':
        gname = 'spec-0301-51942-0531'
        T_OIII = 12418.20420184761
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0327-52294-0042':
        gname = 'spec-0327-52294-0042'
        T_OIII = 11244.419487957326
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0287-52023-0230':
        gname = 'spec-0287-52023-0230'
        T_OIII = 10479.810702407312
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-0278-51900-0081':
    #	 gname = 'spec-0278-51900-0081'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0304-51609-0020':
        gname = 'spec-0304-51609-0020'
        T_OIII = 12041.364818561056
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0304-51609-0034':
        gname = 'spec-0304-51609-0034'
        T_OIII = 13324.863162352554
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0380-51792-0253':
        gname = 'spec-0380-51792-0253'
        T_OIII = 9442.50295571747
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0394-51913-0493':
        gname = 'spec-0394-51913-0493'
        T_OIII = 9961.178300666343
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0394-51913-0599':
        gname = 'spec-0394-51913-0599'
        T_OIII = 13210.62507129942
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0280-51612-0192':
        gname = 'spec-0280-51612-0192'
        T_OIII = 11954.371659803212
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0281-51614-0129':
        gname = 'spec-0281-51614-0129'
        T_OIII = 11088.42493598055
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0358-51818-0311':
        gname = 'spec-0358-51818-0311'
        T_OIII = 11803.643189878701
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0358-51818-0364':
        gname = 'spec-0358-51818-0364'
        T_OIII = 14675.175232451153
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0358-51818-0403':
        gname = 'spec-0358-51818-0403'
        T_OIII = 11644.256719254627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0358-51818-0504':
        gname = 'spec-0358-51818-0504'
        T_OIII = 12497.240786234663
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0328-52282-0491':
        gname = 'spec-0328-52282-0491'
        T_OIII = 11754.681060842748
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0292-51609-0566':
        gname = 'spec-0292-51609-0566'
        T_OIII = 10998.171651487011
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0287-52023-0236':
        gname = 'spec-0287-52023-0236'
        T_OIII = 10962.995862035505
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0304-51609-0583':
        gname = 'spec-0304-51609-0583'
        T_OIII = 10748.662843508935
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0339-51692-0437':
        gname = 'spec-0339-51692-0437'
        T_OIII = 10405.237160936484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0281-51614-0499':
        gname = 'spec-0281-51614-0499'
        T_OIII = 10879.032297991607
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0282-51658-0296':
        gname = 'spec-0282-51658-0296'
        T_OIII = 11879.139612140378
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0296-51984-0416':
        gname = 'spec-0296-51984-0416'
        T_OIII = 10491.11666577711
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0306-51637-0515':
        gname = 'spec-0306-51637-0515'
        T_OIII = 10695.716664527605
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0381-51811-0370':
        gname = 'spec-0381-51811-0370'
        T_OIII = 9511.216994570748
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0305-51613-0594':
        gname = 'spec-0305-51613-0594'
        T_OIII = 12014.111667781179
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0339-51692-0622':
        gname = 'spec-0339-51692-0622'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0282-51658-0543':
        gname = 'spec-0282-51658-0543'
        T_OIII = 12128.991035239995
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0308-51662-0081':
        gname = 'spec-0308-51662-0081'
        T_OIII = 12783.64316161374
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0308-51662-0130':
        gname = 'spec-0308-51662-0130'
        T_OIII = 11497.438379600459
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0446-51899-0086':
        gname = 'spec-0446-51899-0086'
        T_OIII = 10584.81468233188
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0446-51899-0283':
        gname = 'spec-0446-51899-0283'
        T_OIII = 10795.543096281908
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0446-51899-0352':
        gname = 'spec-0446-51899-0352'
        T_OIII = 10918.544448136336
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0351-51695-0217':
        gname = 'spec-0351-51695-0217'
        T_OIII = 10371.147410700316
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0309-51994-0288':
        gname = 'spec-0309-51994-0288'
        T_OIII = 11964.650426800059
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0420-51871-0156':
        gname = 'spec-0420-51871-0156'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0420-51871-0474':
        gname = 'spec-0420-51871-0474'
        T_OIII = 12714.31359382933
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0423-51821-0402':
        gname = 'spec-0423-51821-0402'
        T_OIII = 12508.572716426379
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0383-51818-0266':
        gname = 'spec-0383-51818-0266'
        T_OIII = 13525.629567438933
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0384-51821-0281':
        gname = 'spec-0384-51821-0281'
        T_OIII = 11798.295322805352
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0467-51901-0628':
        gname = 'spec-0467-51901-0628'
        T_OIII = 12278.308640521916
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0309-51994-0489':
        gname = 'spec-0309-51994-0489'
        T_OIII = 11803.643189878701
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0310-51990-0136':
        gname = 'spec-0310-51990-0136'
        T_OIII = 10637.710554366235
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0443-51873-0542':
        gname = 'spec-0443-51873-0542'
        T_OIII = 12216.682187044686
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0456-51910-0195':
        gname = 'spec-0456-51910-0195'
        T_OIII = 12025.005518437003
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0456-51910-0306':
        gname = 'spec-0456-51910-0306'
        T_OIII = 15369.435600871382
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0456-51910-0546':
        gname = 'spec-0456-51910-0546'
        T_OIII = 10968.136502632646
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0395-51783-0570':
        gname = 'spec-0395-51783-0570'
        T_OIII = 13835.599784295056
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0364-52000-0066':
        gname = 'spec-0364-52000-0066'
        T_OIII = 11590.339716453223
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0364-52000-0187':
        gname = 'spec-0364-52000-0187'
        T_OIII = 10913.085976236822
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0275-51910-0015':
        gname = 'spec-0275-51910-0015'
        T_OIII = 23309.145018006344
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0275-51910-0445':
        gname = 'spec-0275-51910-0445'
        T_OIII = 11739.628458687444
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0501-52235-0361':
        gname = 'spec-0501-52235-0361'
        T_OIII = 12536.947516637038
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0367-51997-0561':
        gname = 'spec-0367-51997-0561'
        T_OIII = 10367.582572914509
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0312-51689-0508':
        gname = 'spec-0312-51689-0508'
        T_OIII = 10839.663134855515
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0331-52368-0006':
        gname = 'spec-0331-52368-0006'
        T_OIII = 13999.472075468642
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0345-51690-0572':
        gname = 'spec-0345-51690-0572'
        T_OIII = 17909.411166003403
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0375-52140-0348':
        gname = 'spec-0375-52140-0348'
        T_OIII = 11275.563602529657
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0465-51910-0524':
        gname = 'spec-0465-51910-0524'
        T_OIII = 11970.63486679113
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0353-51703-0486':
        gname = 'spec-0353-51703-0486'
        T_OIII = 10389.964108513297
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0335-52000-0604':
        gname = 'spec-0335-52000-0604'
        T_OIII = 10115.847071186598
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0336-51999-0033':
        gname = 'spec-0336-51999-0033'
        T_OIII = 10734.563020635222
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0505-52317-0228':
        gname = 'spec-0505-52317-0228'
        T_OIII = 14869.309006222302
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0507-52353-0087':
        gname = 'spec-0507-52353-0087'
        T_OIII = 12057.181122473297
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0412-52258-0441':
        gname = 'spec-0412-52258-0441'
        T_OIII = 11109.063795966624
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0412-52258-0518':
        gname = 'spec-0412-52258-0518'
        T_OIII = 13537.893993857477
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0472-51955-0235':
        gname = 'spec-0472-51955-0235'
        T_OIII = 12351.049523795427
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0414-51869-0459':
        gname = 'spec-0414-51869-0459'
        T_OIII = 10734.563020635222
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0414-51869-0524':
        gname = 'spec-0414-51869-0524'
        T_OIII = 10699.39432946481
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0375-52140-0118':
        gname = 'spec-0375-52140-0118'
        T_OIII = 13766.802451262882
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0390-51900-0430':
        gname = 'spec-0390-51900-0430'
        T_OIII = 13174.753759947303
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0390-51900-0445':
        gname = 'spec-0390-51900-0445'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0530-52026-0184':
        gname = 'spec-0530-52026-0184'
        T_OIII = 18948.906672371362
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0394-51913-0402':
        gname = 'spec-0394-51913-0402'
        T_OIII = 10766.229548623747
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0394-51913-0469':
        gname = 'spec-0394-51913-0469'
        T_OIII = 11347.32776582553
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0394-51913-0472':
        gname = 'spec-0394-51913-0472'
        T_OIII = 15675.418925537853
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0341-51690-0256':
        gname = 'spec-0341-51690-0256'
        T_OIII = 9831.125106270674
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0342-51691-0081':
        gname = 'spec-0342-51691-0081'
        T_OIII = 12485.919122013713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0481-51908-0483':
        gname = 'spec-0481-51908-0483'
        T_OIII = 10319.57702029705
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0415-51810-0285':
        gname = 'spec-0415-51810-0285'
        T_OIII = 11758.906763478532
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0530-52026-0525':
        gname = 'spec-0530-52026-0525'
        T_OIII = 10607.665178573705
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0348-51671-0331':
        gname = 'spec-0348-51671-0331'
        T_OIII = 12345.260727623503
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0487-51943-0583':
        gname = 'spec-0487-51943-0583'
        T_OIII = 11607.377298512543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0439-51877-0086':
        gname = 'spec-0439-51877-0086'
        T_OIII = 10571.26017735626
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0450-51908-0377':
        gname = 'spec-0450-51908-0377'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0546-52205-0419':
        gname = 'spec-0546-52205-0419'
        T_OIII = 10651.683161756944
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0550-51959-0092':
        gname = 'spec-0550-51959-0092'
        T_OIII = 18353.083866145415
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0454-51908-0570':
        gname = 'spec-0454-51908-0570'
        T_OIII = 13587.06300866288
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0455-51909-0047':
        gname = 'spec-0455-51909-0047'
        T_OIII = 11617.902340920347
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0349-51699-0337':
        gname = 'spec-0349-51699-0337'
        T_OIII = 23489.41032682246
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0492-51955-0449':
        gname = 'spec-0492-51955-0449'
        T_OIII = 11349.278624789184
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0493-51957-0219':
        gname = 'spec-0493-51957-0219'
        T_OIII = 13174.753759947303
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0497-51989-0198':
        gname = 'spec-0497-51989-0198'
        T_OIII = 12150.997043837559
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0450-51908-0520':
        gname = 'spec-0450-51908-0520'
        T_OIII = 15011.489873224315
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0450-51908-0545':
        gname = 'spec-0450-51908-0545'
        T_OIII = 11291.079771779818
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0436-51883-0600':
        gname = 'spec-0436-51883-0600'
        T_OIII = 12520.306217341384
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0564-52224-0216':
        gname = 'spec-0564-52224-0216'
        T_OIII = 10908.652977199827
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0564-52224-0223':
        gname = 'spec-0564-52224-0223'
        T_OIII = 20007.893234940188
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0569-52264-0609':
        gname = 'spec-0569-52264-0609'
        T_OIII = 11954.371659803212
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0570-52266-0061':
        gname = 'spec-0570-52266-0061'
        T_OIII = 11398.867355156632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0351-51695-0049':
        gname = 'spec-0351-51695-0049'
        T_OIII = 11466.218949055712
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0452-51911-0487':
        gname = 'spec-0452-51911-0487'
        T_OIII = 10657.01088604571
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0502-51957-0007':
        gname = 'spec-0502-51957-0007'
        T_OIII = 11981.489294658955
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0494-51915-0537':
        gname = 'spec-0494-51915-0537'
        T_OIII = 22448.6955852049
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0409-51871-0033':
        gname = 'spec-0409-51871-0033'
        T_OIII = 11347.32776582553
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0409-51871-0050':
        gname = 'spec-0409-51871-0050'
        T_OIII = 12841.706498841864
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0437-51876-0460':
        gname = 'spec-0437-51876-0460'
        T_OIII = 10643.031289881616
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0575-52319-0521':
        gname = 'spec-0575-52319-0521'
        T_OIII = 15066.010921192563
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0488-51914-0439':
        gname = 'spec-0488-51914-0439'
        T_OIII = 10070.1086130086
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0524-52027-0016':
        gname = 'spec-0524-52027-0016'
        T_OIII = 11632.43530276649
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0507-52353-0377':
        gname = 'spec-0507-52353-0377'
        T_OIII = 11083.054724220938
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0441-51868-0634':
        gname = 'spec-0441-51868-0634'
        T_OIII = 12211.71965486194
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0442-51882-0156':
        gname = 'spec-0442-51882-0156'
        T_OIII = 11497.438379600459
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0442-51882-0413':
        gname = 'spec-0442-51882-0413'
        T_OIII = 11461.56126215855
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0516-52017-0315':
        gname = 'spec-0516-52017-0315'
        T_OIII = 11981.489294658955
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0516-52017-0403':
        gname = 'spec-0516-52017-0403'
        T_OIII = 15196.293350535128
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0516-52017-0459':
        gname = 'spec-0516-52017-0459'
        T_OIII = 13501.134036786234
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0518-52282-0335':
        gname = 'spec-0518-52282-0335'
        T_OIII = 10338.300152334337
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0425-51898-0634':
        gname = 'spec-0425-51898-0634'
        T_OIII = 10690.870764829013
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0447-51877-0172':
        gname = 'spec-0447-51877-0172'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0447-51877-0361':
        gname = 'spec-0447-51877-0361'
        T_OIII = 11846.513484698893
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0597-52059-0586':
        gname = 'spec-0597-52059-0586'
        T_OIII = 11440.803853135543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0406-51817-0490':
        gname = 'spec-0406-51817-0490'
        T_OIII = 11681.253314758873
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0469-51913-0009':
        gname = 'spec-0469-51913-0009'
        T_OIII = 10888.896910241861
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0469-51913-0110':
        gname = 'spec-0469-51913-0110'
        T_OIII = 11476.615995485094
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0536-52024-0326':
        gname = 'spec-0536-52024-0326'
        T_OIII = 12812.44172273529
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0525-52295-0626':
        gname = 'spec-0525-52295-0626'
        T_OIII = 11492.22924449282
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0535-51999-0052':
        gname = 'spec-0535-51999-0052'
        T_OIII = 11093.104334738524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0519-52283-0124':
        gname = 'spec-0519-52283-0124'
        T_OIII = 11992.353564828874
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0498-51984-0521':
        gname = 'spec-0498-51984-0521'
        T_OIII = 10575.22554824027
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0564-52224-0283':
        gname = 'spec-0564-52224-0283'
        T_OIII = 10323.770633852138
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0483-51924-0594':
        gname = 'spec-0483-51924-0594'
        T_OIII = 11296.02121098684
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0431-51877-0191':
        gname = 'spec-0431-51877-0191'
        T_OIII = 9799.988194065465
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0463-51908-0299':
        gname = 'spec-0463-51908-0299'
        T_OIII = 12340.245965617973
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0499-51988-0124':
        gname = 'spec-0499-51988-0124'
        T_OIII = 12041.364818561056
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0512-51992-0524':
        gname = 'spec-0512-51992-0524'
        T_OIII = 11093.104334738524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0519-52283-0615':
        gname = 'spec-0519-52283-0615'
        T_OIII = 12311.739172898178
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0568-52254-0244':
        gname = 'spec-0568-52254-0244'
        T_OIII = 10993.01693375783
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0463-51908-0470':
        gname = 'spec-0463-51908-0470'
        T_OIII = 11296.02121098684
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0490-51929-0279':
        gname = 'spec-0490-51929-0279'
        T_OIII = 16521.48652519121
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0490-51929-0363':
        gname = 'spec-0490-51929-0363'
        T_OIII = 12139.989053287203
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0490-51929-0396':
        gname = 'spec-0490-51929-0396'
        T_OIII = 12112.512686851956
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0490-51929-0401':
        gname = 'spec-0490-51929-0401'
        T_OIII = 10498.8245514363
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0490-51929-0519':
        gname = 'spec-0490-51929-0519'
        T_OIII = 11760.928034558468
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0494-51915-0007':
        gname = 'spec-0494-51915-0007'
        T_OIII = 13300.731228795645
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0630-52050-0395':
        gname = 'spec-0630-52050-0395'
        T_OIII = 11873.386451018981
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0424-51893-0279':
        gname = 'spec-0424-51893-0279'
        T_OIII = 20624.699696250675
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0531-52028-0218':
        gname = 'spec-0531-52028-0218'
        T_OIII = 11487.022469483838
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0469-51913-0498':
        gname = 'spec-0469-51913-0498'
        T_OIII = 12046.44636806152
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0472-51955-0546':
        gname = 'spec-0472-51955-0546'
        T_OIII = 10268.744448852292
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0496-51988-0541':
        gname = 'spec-0496-51988-0541'
        T_OIII = 12519.914921897609
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0483-51924-0495':
        gname = 'spec-0483-51924-0495'
        T_OIII = 10824.93645684694
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0495-51988-0466':
        gname = 'spec-0495-51988-0466'
        T_OIII = 11148.540300786497
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0541-51959-0025':
        gname = 'spec-0541-51959-0025'
        T_OIII = 11739.628458687444
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0552-51992-0027':
        gname = 'spec-0552-51992-0027'
        T_OIII = 11628.436926958453
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0603-52056-0486':
        gname = 'spec-0603-52056-0486'
        T_OIII = 10604.019043372815
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0607-52368-0239':
        gname = 'spec-0607-52368-0239'
        T_OIII = 10958.200104611962
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0544-52201-0067':
        gname = 'spec-0544-52201-0067'
        T_OIII = 11239.851898105111
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0555-52266-0517':
        gname = 'spec-0555-52266-0517'
        T_OIII = 11586.355810012077
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0585-52027-0255':
        gname = 'spec-0585-52027-0255'
        T_OIII = 11296.02121098684
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0473-51929-0422':
        gname = 'spec-0473-51929-0422'
        T_OIII = 11607.377298512543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0503-51999-0123':
        gname = 'spec-0503-51999-0123'
        T_OIII = 14056.796001509865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0503-51999-0273':
        gname = 'spec-0503-51999-0273'
        T_OIII = 9949.044260494311
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0498-51984-0236':
        gname = 'spec-0498-51984-0236'
        T_OIII = 10869.176622420213
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0664-52174-0447':
        gname = 'spec-0664-52174-0447'
        T_OIII = 12990.964491645
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0459-51924-0253':
        gname = 'spec-0459-51924-0253'
        T_OIII = 10604.019043372815
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0553-51999-0342':
        gname = 'spec-0553-51999-0342'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0648-52559-0288':
        gname = 'spec-0648-52559-0288'
        T_OIII = 11445.453104749757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0649-52201-0019':
        gname = 'spec-0649-52201-0019'
        T_OIII = 15369.435600871382
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0649-52201-0117':
        gname = 'spec-0649-52201-0117'
        T_OIII = 10908.652977199827
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0656-52148-0243':
        gname = 'spec-0656-52148-0243'
        T_OIII = 14815.499785332746
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0566-52238-0497':
        gname = 'spec-0566-52238-0497'
        T_OIII = 15251.485596321565
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0570-52266-0586':
        gname = 'spec-0570-52266-0586'
        T_OIII = 14275.058065697667
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0587-52026-0495':
        gname = 'spec-0587-52026-0495'
        T_OIII = 16770.416620452357
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0587-52026-0546':
        gname = 'spec-0587-52026-0546'
        T_OIII = 11782.26625491757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0527-52342-0103':
        gname = 'spec-0527-52342-0103'
        T_OIII = 12497.240786234663
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0515-52051-0378':
        gname = 'spec-0515-52051-0378'
        T_OIII = 11670.67088060627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0519-52283-0507':
        gname = 'spec-0519-52283-0507'
        T_OIII = 11209.331901904156
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0519-52283-0515':
        gname = 'spec-0519-52283-0515'
        T_OIII = 17990.755753531575
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0460-51924-0128':
        gname = 'spec-0460-51924-0128'
        T_OIII = 10540.57888021582
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0520-52288-0551':
        gname = 'spec-0520-52288-0551'
        T_OIII = 22479.235691644666
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0545-52202-0339':
        gname = 'spec-0545-52202-0339'
        T_OIII = 10475.062623048356
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0625-52145-0432':
        gname = 'spec-0625-52145-0432'
        T_OIII = 11750.273420439087
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0637-52174-0526':
        gname = 'spec-0637-52174-0526'
        T_OIII = 10319.57702029705
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0614-53437-0055':
        gname = 'spec-0614-53437-0055'
        T_OIII = 10456.091807953026
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0523-52026-0114':
        gname = 'spec-0523-52026-0114'
        T_OIII = 11771.592309797865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0524-52027-0260':
        gname = 'spec-0524-52027-0260'
        T_OIII = 10903.710602618274
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0529-52025-0257':
        gname = 'spec-0529-52025-0257'
        T_OIII = 14181.564284025484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0682-52525-0172':
        gname = 'spec-0682-52525-0172'
        T_OIII = 16861.863612637975
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0584-52049-0421':
        gname = 'spec-0584-52049-0421'
        T_OIII = 11565.372392383904
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0691-52199-0252':
        gname = 'spec-0691-52199-0252'
        T_OIII = 12311.739172898178
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0704-52205-0494':
        gname = 'spec-0704-52205-0494'
        T_OIII = 13446.181170901009
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0633-52079-0191':
        gname = 'spec-0633-52079-0191'
        T_OIII = 15293.011289960772
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0541-51959-0600':
        gname = 'spec-0541-51959-0600'
        T_OIII = 14294.478427611624
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0542-51993-0013':
        gname = 'spec-0542-51993-0013'
        T_OIII = 21806.967099901325
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0539-52017-0155':
        gname = 'spec-0539-52017-0155'
        T_OIII = 11178.894777806525
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0539-52017-0356':
        gname = 'spec-0539-52017-0356'
        T_OIII = 12085.098507516668
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0615-52347-0272':
        gname = 'spec-0615-52347-0272'
        T_OIII = 11117.747035927965
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0622-52054-0445':
        gname = 'spec-0622-52054-0445'
        T_OIII = 9654.522994844869
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0641-52176-0213':
        gname = 'spec-0641-52176-0213'
        T_OIII = 11409.20333022049
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0566-52238-0338':
        gname = 'spec-0566-52238-0338'
        T_OIII = 12289.442052990076
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0544-52201-0610':
        gname = 'spec-0544-52201-0610'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0553-51999-0134':
        gname = 'spec-0553-51999-0134'
        T_OIII = 10884.813922643169
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0542-51993-0269':
        gname = 'spec-0542-51993-0269'
        T_OIII = 10489.313318651464
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0621-52055-0618':
        gname = 'spec-0621-52055-0618'
        T_OIII = 11316.515936324256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0621-52055-0624':
        gname = 'spec-0621-52055-0624'
        T_OIII = 10744.296634708051
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0664-52174-0112':
        gname = 'spec-0664-52174-0112'
        T_OIII = 11445.453104749757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0639-52146-0143':
        gname = 'spec-0639-52146-0143'
        T_OIII = 11948.955502303877
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0639-52146-0242':
        gname = 'spec-0639-52146-0242'
        T_OIII = 10291.555881379016
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0658-52146-0017':
        gname = 'spec-0658-52146-0017'
        T_OIII = 11905.71448917048
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0590-52057-0200':
        gname = 'spec-0590-52057-0200'
        T_OIII = 10978.081910527662
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0660-52177-0414':
        gname = 'spec-0660-52177-0414'
        T_OIII = 13258.605460273811
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0666-52149-0492':
        gname = 'spec-0666-52149-0492'
        T_OIII = 10780.370965406275
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0556-51991-0224':
        gname = 'spec-0556-51991-0224'
        T_OIII = 10908.652977199827
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0629-52051-0497':
        gname = 'spec-0629-52051-0497'
        T_OIII = 10399.385255947842
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0686-52519-0185':
        gname = 'spec-0686-52519-0185'
        T_OIII = 12041.364818561056
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0650-52143-0330':
        gname = 'spec-0650-52143-0330'
        T_OIII = 11507.86373439115
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0677-52606-0533':
        gname = 'spec-0677-52606-0533'
        T_OIII = 17271.75640239755
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0608-52081-0139':
        gname = 'spec-0608-52081-0139'
        T_OIII = 20897.538059098697
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0615-52347-0590':
        gname = 'spec-0615-52347-0590'
        T_OIII = 10618.113322221714
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0569-52264-0472':
        gname = 'spec-0569-52264-0472'
        T_OIII = 12559.693693465491
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0653-52145-0056':
        gname = 'spec-0653-52145-0056'
        T_OIII = 10918.544448136336
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0668-52162-0089':
        gname = 'spec-0668-52162-0089'
        T_OIII = 11921.352616699249
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0786-52319-0443':
        gname = 'spec-0786-52319-0443'
        T_OIII = 10408.814946050896
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0795-52378-0405':
        gname = 'spec-0795-52378-0405'
        T_OIII = 10254.312684436489
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0662-52147-0333':
        gname = 'spec-0662-52147-0333'
        T_OIII = 11691.845344582025
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0692-52201-0465':
        gname = 'spec-0692-52201-0465'
        T_OIII = 10666.674179366131
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0665-52168-0324':
        gname = 'spec-0665-52168-0324'
        T_OIII = 18966.088680087632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0620-52375-0428':
        gname = 'spec-0620-52375-0428'
        T_OIII = 11728.993140559112
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0690-52261-0362':
        gname = 'spec-0690-52261-0362'
        T_OIII = 17326.632913683254
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0582-52045-0440':
        gname = 'spec-0582-52045-0440'
        T_OIII = 11260.244713697612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0582-52045-0445':
        gname = 'spec-0582-52045-0445'
        T_OIII = 12559.693693465491
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0679-52177-0605':
        gname = 'spec-0679-52177-0605'
        T_OIII = 12228.141882422768
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0712-52199-0158':
        gname = 'spec-0712-52199-0158'
        T_OIII = 9979.251178673212
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0673-52162-0073':
        gname = 'spec-0673-52162-0073'
        T_OIII = 11260.244713697612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0721-52228-0555':
        gname = 'spec-0721-52228-0555'
        T_OIII = 11695.865523035109
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0624-52377-0422':
        gname = 'spec-0624-52377-0422'
        T_OIII = 11378.223487652196
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0630-52050-0463':
        gname = 'spec-0630-52050-0463'
        T_OIII = 10300.88779670391
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0640-52200-0270':
        gname = 'spec-0640-52200-0270'
        T_OIII = 10456.091807953026
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0805-52586-0482':
        gname = 'spec-0805-52586-0482'
        T_OIII = 11586.355810012077
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0707-52177-0525':
        gname = 'spec-0707-52177-0525'
        T_OIII = 12548.315451087996
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0708-52175-0155':
        gname = 'spec-0708-52175-0155'
        T_OIII = 14418.089075810085
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0602-52072-0019':
        gname = 'spec-0602-52072-0019'
        T_OIII = 11560.132478468087
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0684-52523-0560':
        gname = 'spec-0684-52523-0560'
        T_OIII = 12325.406755766515
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0686-52519-0406':
        gname = 'spec-0686-52519-0406'
        T_OIII = 11848.550164831573
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0695-52202-0261':
        gname = 'spec-0695-52202-0261'
        T_OIII = 11728.993140559112
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0663-52145-0406':
        gname = 'spec-0663-52145-0406'
        T_OIII = 12929.29636094603
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0834-52316-0157':
        gname = 'spec-0834-52316-0157'
        T_OIII = 11174.70305029536
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0602-52072-0369':
        gname = 'spec-0602-52072-0369'
        T_OIII = 13667.343822670808
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0602-52072-0500':
        gname = 'spec-0602-52072-0500'
        T_OIII = 11306.263929839473
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0831-52294-0526':
        gname = 'spec-0831-52294-0526'
        T_OIII = 13097.366378029848
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0831-52294-0572':
        gname = 'spec-0831-52294-0572'
        T_OIII = 11803.643189878701
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0840-52374-0570':
        gname = 'spec-0840-52374-0570'
        T_OIII = 13483.002198146052
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0749-52226-0500':
        gname = 'spec-0749-52226-0500'
        T_OIII = 11058.834543527928
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0755-52235-0416':
        gname = 'spec-0755-52235-0416'
        T_OIII = 11424.724868332823
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0742-52263-0179':
        gname = 'spec-0742-52263-0179'
        T_OIII = 13061.80260230364
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0675-52590-0094':
        gname = 'spec-0675-52590-0094'
        T_OIII = 11419.370229971493
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0687-52518-0274':
        gname = 'spec-0687-52518-0274'
        T_OIII = 12091.143192837555
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0708-52175-0404':
        gname = 'spec-0708-52175-0404'
        T_OIII = 12691.287393325083
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0606-52365-0604':
        gname = 'spec-0606-52365-0604'
        T_OIII = 11178.894777806525
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0750-52235-0312':
        gname = 'spec-0750-52235-0312'
        T_OIII = 12311.739172898178
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0769-54530-0086':
        gname = 'spec-0769-54530-0086'
        T_OIII = 12289.442052990076
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0771-52370-0490':
        gname = 'spec-0771-52370-0490'
        T_OIII = 11651.537579641648
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0634-52164-0230':
        gname = 'spec-0634-52164-0230'
        T_OIII = 16976.873896417892
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0879-52365-0014':
        gname = 'spec-0879-52365-0014'
        T_OIII = 12005.666589400918
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0879-52365-0097':
        gname = 'spec-0879-52365-0097'
        T_OIII = 22560.879252529303
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0904-52381-0501':
        gname = 'spec-0904-52381-0501'
        T_OIII = 10958.200104611962
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0762-52232-0217':
        gname = 'spec-0762-52232-0217'
        T_OIII = 11113.230901627261
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0762-52232-0575':
        gname = 'spec-0762-52232-0575'
        T_OIII = 12952.754388580213
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0764-52238-0477':
        gname = 'spec-0764-52238-0477'
        T_OIII = 10953.577612439196
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0881-52368-0567':
        gname = 'spec-0881-52368-0567'
        T_OIII = 10164.652123175232
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0894-52615-0185':
        gname = 'spec-0894-52615-0185'
        T_OIII = 11868.006985422713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0755-52235-0300':
        gname = 'spec-0755-52235-0300'
        T_OIII = 11535.230225304234
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0730-52466-0376':
        gname = 'spec-0730-52466-0376'
        T_OIII = 10963.509818993989
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0726-52226-0532':
        gname = 'spec-0726-52226-0532'
        T_OIII = 12311.739172898178
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0737-52518-0544':
        gname = 'spec-0737-52518-0544'
        T_OIII = 16604.049939221743
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0651-52141-0452':
        gname = 'spec-0651-52141-0452'
        T_OIII = 11133.046037033146
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0653-52145-0462':
        gname = 'spec-0653-52145-0462'
        T_OIII = 11739.628458687444
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0663-52145-0288':
        gname = 'spec-0663-52145-0288'
        T_OIII = 11239.851898105111
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0782-52320-0022':
        gname = 'spec-0782-52320-0022'
        T_OIII = 10277.573857418769
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0662-52147-0466':
        gname = 'spec-0662-52147-0466'
        T_OIII = 10083.808392923831
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0758-52253-0176':
        gname = 'spec-0758-52253-0176'
        T_OIII = 12725.842355781295
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0903-52400-0631':
        gname = 'spec-0903-52400-0631'
        T_OIII = 10824.93645684694
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0761-54524-0477':
        gname = 'spec-0761-54524-0477'
        T_OIII = 14018.627050134624
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0804-52286-0085':
        gname = 'spec-0804-52286-0085'
        T_OIII = 10883.453259938144
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0735-52519-0461':
        gname = 'spec-0735-52519-0461'
        T_OIII = 11728.993140559112
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0740-52263-0259':
        gname = 'spec-0740-52263-0259'
        T_OIII = 9576.09047420312
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0740-52263-0267':
        gname = 'spec-0740-52263-0267'
        T_OIII = 10189.93879287919
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0676-52178-0192':
        gname = 'spec-0676-52178-0192'
        T_OIII = 15732.351333605193
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0922-52426-0212':
        gname = 'spec-0922-52426-0212'
        T_OIII = 12272.745717156748
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-0933-52642-0080':
    #	 gname = 'spec-0933-52642-0080'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0802-52289-0423':
        gname = 'spec-0802-52289-0423'
        T_OIII = 13109.652189038374
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0785-52339-0079':
        gname = 'spec-0785-52339-0079'
        T_OIII = 9827.131610001748
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0767-52252-0463':
        gname = 'spec-0767-52252-0463'
        T_OIII = 12830.072771726329
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0911-52426-0253':
        gname = 'spec-0911-52426-0253'
        T_OIII = 10483.414129234881
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0766-52247-0492':
        gname = 'spec-0766-52247-0492'
        T_OIII = 10070.1086130086
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0775-52295-0029':
        gname = 'spec-0775-52295-0029'
        T_OIII = 11927.31540013229
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0679-52177-0055':
        gname = 'spec-0679-52177-0055'
        T_OIII = 11798.295322805352
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0841-52375-0333':
        gname = 'spec-0841-52375-0333'
        T_OIII = 12239.80359412302
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0703-52209-0285':
        gname = 'spec-0703-52209-0285'
        T_OIII = 11146.624276640516
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0709-52205-0584':
        gname = 'spec-0709-52205-0584'
        T_OIII = 9888.749910497594
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0775-52295-0291':
        gname = 'spec-0775-52295-0291'
        T_OIII = 12222.792687463036
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0838-52378-0595':
        gname = 'spec-0838-52378-0595'
        T_OIII = 14444.248266730488
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0792-52353-0110':
        gname = 'spec-0792-52353-0110'
        T_OIII = 11575.85934665115
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0764-52238-0633':
        gname = 'spec-0764-52238-0633'
        T_OIII = 14320.413347858608
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0779-52342-0218':
        gname = 'spec-0779-52342-0218'
        T_OIII = 10060.985768187482
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0783-52325-0594':
        gname = 'spec-0783-52325-0594'
        T_OIII = 13501.134036786234
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0964-52646-0570':
        gname = 'spec-0964-52646-0570'
        T_OIII = 11665.930143395257
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0971-52644-0599':
        gname = 'spec-0971-52644-0599'
        T_OIII = 11275.563602529657
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0713-52178-0481':
        gname = 'spec-0713-52178-0481'
        T_OIII = 10617.283727448379
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0951-52398-0017':
        gname = 'spec-0951-52398-0017'
        T_OIII = 14191.207589965601
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0849-52439-0598':
        gname = 'spec-0849-52439-0598'
        T_OIII = 11367.915578253484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0850-52338-0010':
        gname = 'spec-0850-52338-0010'
        T_OIII = 11194.102994920871
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0864-52320-0124':
        gname = 'spec-0864-52320-0124'
        T_OIII = 12294.436169578046
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0780-52370-0398':
        gname = 'spec-0780-52370-0398'
        T_OIII = 19017.72823894171
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0884-52374-0216':
        gname = 'spec-0884-52374-0216'
        T_OIII = 10734.563020635222
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-0865-52323-0080':
    #	 gname = 'spec-0865-52323-0080'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0844-52378-0299':
        gname = 'spec-0844-52378-0299'
        T_OIII = 10551.290701921758
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0875-52354-0226':
        gname = 'spec-0875-52354-0226'
        T_OIII = 9584.773638000595
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0908-52373-0027':
        gname = 'spec-0908-52373-0027'
        T_OIII = 11461.56126215855
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0909-52379-0616':
        gname = 'spec-0909-52379-0616'
        T_OIII = 13587.06300866288
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1002-52646-0525':
        gname = 'spec-1002-52646-0525'
        T_OIII = 9799.988194065465
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0847-52426-0522':
        gname = 'spec-0847-52426-0522'
        T_OIII = 11296.02121098684
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0854-52373-0514':
        gname = 'spec-0854-52373-0514'
        T_OIII = 11782.26625491757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0778-54525-0497':
        gname = 'spec-0778-54525-0497'
        T_OIII = 9812.247107106516
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0989-52468-0143':
        gname = 'spec-0989-52468-0143'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0741-52261-0279':
        gname = 'spec-0741-52261-0279'
        T_OIII = 10973.622509525656
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0898-52606-0045':
        gname = 'spec-0898-52606-0045'
        T_OIII = 11539.196552486676
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0843-52378-0001':
        gname = 'spec-0843-52378-0001'
        T_OIII = 11560.132478468087
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0845-52381-0001':
        gname = 'spec-0845-52381-0001'
        T_OIII = 10551.290701921758
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0846-52407-0267':
        gname = 'spec-0846-52407-0267'
        T_OIII = 11285.78777135236
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1007-52706-0478':
        gname = 'spec-1007-52706-0478'
        T_OIII = 10948.272708303248
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0890-52583-0314':
        gname = 'spec-0890-52583-0314'
        T_OIII = 12244.968883338277
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0816-52379-0277':
        gname = 'spec-0816-52379-0277'
        T_OIII = 11316.515936324256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0940-52670-0594':
        gname = 'spec-0940-52670-0594'
        T_OIII = 11884.152698441407
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0949-52427-0568':
        gname = 'spec-0949-52427-0568'
        T_OIII = 11219.496014829632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0951-52398-0142':
        gname = 'spec-0951-52398-0142'
        T_OIII = 12334.07674716982
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0945-52652-0423':
        gname = 'spec-0945-52652-0423'
        T_OIII = 11798.295322805352
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0952-52409-0439':
        gname = 'spec-0952-52409-0439'
        T_OIII = 12593.890361666976
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0955-52409-0383':
        gname = 'spec-0955-52409-0383'
        T_OIII = 11771.592309797865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0957-52398-0530':
        gname = 'spec-0957-52398-0530'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0892-52378-0229':
        gname = 'spec-0892-52378-0229'
        T_OIII = 12485.919122013713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0854-52373-0373':
        gname = 'spec-0854-52373-0373'
        T_OIII = 10465.572916981222
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0832-52312-0473':
        gname = 'spec-0832-52312-0473'
        T_OIII = 11560.132478468087
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1039-52707-0119':
        gname = 'spec-1039-52707-0119'
        T_OIII = 10217.204263589458
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0960-52425-0581':
        gname = 'spec-0960-52425-0581'
        T_OIII = 11681.253314758873
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0885-52379-0301':
        gname = 'spec-0885-52379-0301'
        T_OIII = 16350.178172330667
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0967-52636-0244':
        gname = 'spec-0967-52636-0244'
        T_OIII = 10680.68483040906
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0967-52636-0339':
        gname = 'spec-0967-52636-0339'
        T_OIII = 11138.440468991766
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0829-52296-0210':
        gname = 'spec-0829-52296-0210'
        T_OIII = 12656.2329230936
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0897-52605-0553':
        gname = 'spec-0897-52605-0553'
        T_OIII = 14665.318044598422
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1014-52707-0393':
        gname = 'spec-1014-52707-0393'
        T_OIII = 11362.232453727618
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1015-52734-0003':
        gname = 'spec-1015-52734-0003'
        T_OIII = 10198.222339403897
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0960-52425-0368':
        gname = 'spec-0960-52425-0368'
        T_OIII = 11814.346197280407
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0963-52643-0395':
        gname = 'spec-0963-52643-0395'
        T_OIII = 15636.399057801034
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0877-52353-0253':
        gname = 'spec-0877-52353-0253'
        T_OIII = 11649.534764548354
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0836-52376-0378':
        gname = 'spec-0836-52376-0378'
        T_OIII = 11047.952729466202
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0991-52707-0004':
        gname = 'spec-0991-52707-0004'
        T_OIII = 11083.054724220938
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0992-52644-0524':
        gname = 'spec-0992-52644-0524'
        T_OIII = 11792.949878685811
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1017-52706-0153':
        gname = 'spec-1017-52706-0153'
        T_OIII = 11771.592309797865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0848-52669-0528':
        gname = 'spec-0848-52669-0528'
        T_OIII = 10766.229548623747
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1003-52641-0327':
        gname = 'spec-1003-52641-0327'
        T_OIII = 9619.585099600687
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1008-52707-0430':
        gname = 'spec-1008-52707-0430'
        T_OIII = 12474.607714463225
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0952-52409-0264':
        gname = 'spec-0952-52409-0264'
        T_OIII = 10456.091807953026
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0991-52707-0488':
        gname = 'spec-0991-52707-0488'
        T_OIII = 11602.118354884411
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0890-52583-0049':
        gname = 'spec-0890-52583-0049'
        T_OIII = 10556.073333529317
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0890-52583-0065':
        gname = 'spec-0890-52583-0065'
        T_OIII = 11528.742812258617
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0900-52637-0388':
        gname = 'spec-0900-52637-0388'
        T_OIII = 10775.991876465538
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0902-52409-0135':
        gname = 'spec-0902-52409-0135'
        T_OIII = 11083.054724220938
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0902-52409-0276':
        gname = 'spec-0902-52409-0276'
        T_OIII = 10815.129788489168
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0884-52374-0404':
        gname = 'spec-0884-52374-0404'
        T_OIII = 15107.031617171797
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0850-52338-0397':
        gname = 'spec-0850-52338-0397'
        T_OIII = 14509.854076241529
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0853-52374-0577':
        gname = 'spec-0853-52374-0577'
        T_OIII = 11122.091203466698
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0962-52620-0227':
        gname = 'spec-0962-52620-0227'
        T_OIII = 10037.11681453645
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1061-52641-0393':
        gname = 'spec-1061-52641-0393'
        T_OIII = 12795.23478846811
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0962-52620-0028':
        gname = 'spec-0962-52620-0028'
        T_OIII = 12278.308640521916
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1018-52672-0311':
        gname = 'spec-1018-52672-0311'
        T_OIII = 12463.306554291392
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0914-52721-0052':
        gname = 'spec-0914-52721-0052'
        T_OIII = 16078.314456919477
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0893-52589-0177':
        gname = 'spec-0893-52589-0177'
        T_OIII = 12463.306554291392
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-1031-53172-0544':
    #	 gname = 'spec-1031-53172-0544'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0859-52317-0070':
        gname = 'spec-0859-52317-0070'
        T_OIII = 11435.084289814917
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0972-52435-0370':
        gname = 'spec-0972-52435-0370'
        T_OIII = 11713.058226051426
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0926-52413-0279':
        gname = 'spec-0926-52413-0279'
        T_OIII = 10761.351701709767
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0950-52378-0214':
        gname = 'spec-0950-52378-0214'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0922-52426-0547':
        gname = 'spec-0922-52426-0547'
        T_OIII = 11194.102994920871
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0900-52637-0234':
        gname = 'spec-0900-52637-0234'
        T_OIII = 14056.796001509865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0917-52400-0336':
        gname = 'spec-0917-52400-0336'
        T_OIII = 11691.845344582025
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0896-52592-0125':
        gname = 'spec-0896-52592-0125'
        T_OIII = 11398.867355156632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0992-52644-0019':
        gname = 'spec-0992-52644-0019'
        T_OIII = 13692.140912947623
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1096-52974-0102':
        gname = 'spec-1096-52974-0102'
        T_OIII = 11275.563602529657
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1005-52703-0452':
        gname = 'spec-1005-52703-0452'
        T_OIII = 12517.078608827946
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0957-52398-0210':
        gname = 'spec-0957-52398-0210'
        T_OIII = 11052.96048533375
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0959-52411-0137':
        gname = 'spec-0959-52411-0137'
        T_OIII = 11750.273420439087
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0937-52707-0275':
        gname = 'spec-0937-52707-0275'
        T_OIII = 11445.453104749757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0944-52614-0487':
        gname = 'spec-0944-52614-0487'
        T_OIII = 11260.244713697612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1177-52824-0556':
        gname = 'spec-1177-52824-0556'
        T_OIII = 14204.075532297544
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1177-52824-0616':
        gname = 'spec-1177-52824-0616'
        T_OIII = 11522.979285904317
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1184-52641-0581':
        gname = 'spec-1184-52641-0581'
        T_OIII = 10489.313318651464
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1195-52724-0060':
        gname = 'spec-1195-52724-0060'
        T_OIII = 12030.456147473606
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1054-52516-0499':
        gname = 'spec-1054-52516-0499'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0899-52620-0594':
        gname = 'spec-0899-52620-0594'
        T_OIII = 12651.091841232617
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0900-52637-0085':
        gname = 'spec-0900-52637-0085'
        T_OIII = 13544.03037674597
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0999-52636-0150':
        gname = 'spec-0999-52636-0150'
        T_OIII = 12882.507679404664
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1020-52721-0287':
        gname = 'spec-1020-52721-0287'
        T_OIII = 11466.218949055712
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0930-52618-0274':
        gname = 'spec-0930-52618-0274'
        T_OIII = 10394.186325590279
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1040-52722-0318':
        gname = 'spec-1040-52722-0318'
        T_OIII = 13002.744108062774
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1040-52722-0358':
        gname = 'spec-1040-52722-0358'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1040-52722-0494':
        gname = 'spec-1040-52722-0494'
        T_OIII = 13848.145277624182
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1025-53239-0103':
        gname = 'spec-1025-53239-0103'
        T_OIII = 12074.150216632315
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1207-52672-0506':
        gname = 'spec-1207-52672-0506'
        T_OIII = 10888.896910241861
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1173-52790-0240':
        gname = 'spec-1173-52790-0240'
        T_OIII = 10254.312684436489
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1176-52791-0591':
        gname = 'spec-1176-52791-0591'
        T_OIII = 12406.954139485417
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1108-53227-0397':
        gname = 'spec-1108-53227-0397'
        T_OIII = 11681.253314758873
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1012-52649-0392':
        gname = 'spec-1012-52649-0392'
        T_OIII = 10628.073500152708
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0938-52708-0627':
        gname = 'spec-0938-52708-0627'
        T_OIII = 14056.796001509865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1215-52725-0501':
        gname = 'spec-1215-52725-0501'
        T_OIII = 12983.607679774886
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1185-52642-0580':
        gname = 'spec-1185-52642-0580'
        T_OIII = 12395.714268927213
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1185-52642-0630':
        gname = 'spec-1185-52642-0630'
        T_OIII = 11771.592309797865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1040-52722-0021':
        gname = 'spec-1040-52722-0021'
        T_OIII = 12272.745717156748
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1028-52884-0562':
        gname = 'spec-1028-52884-0562'
        T_OIII = 13150.893681349155
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-1034-52813-0083':
    #	 gname = 'spec-1034-52813-0083'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1034-52813-0521':
        gname = 'spec-1034-52813-0521'
        T_OIII = 11750.273420439087
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1240-52734-0340':
        gname = 'spec-1240-52734-0340'
        T_OIII = 11073.014217978294
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0949-52427-0359':
        gname = 'spec-0949-52427-0359'
        T_OIII = 12662.563276571509
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1048-52736-0424':
        gname = 'spec-1048-52736-0424'
        T_OIII = 12806.836926083124
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1052-52466-0115':
        gname = 'spec-1052-52466-0115'
        T_OIII = 12003.227686225393
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0963-52643-0165':
        gname = 'spec-0963-52643-0165'
        T_OIII = 12485.919122013713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0960-52425-0063':
        gname = 'spec-0960-52425-0063'
        T_OIII = 11483.074075914252
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1096-52974-0288':
        gname = 'spec-1096-52974-0288'
        T_OIII = 12139.989053287203
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1156-52641-0423':
        gname = 'spec-1156-52641-0423'
        T_OIII = 11229.669344110953
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-1159-52669-0305':
    #	 gname = 'spec-1159-52669-0305'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1050-52721-0274':
        gname = 'spec-1050-52721-0274'
        T_OIII = 10948.272708303248
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1050-52721-0402':
        gname = 'spec-1050-52721-0402'
        T_OIII = 10460.34089764467
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0976-52413-0238':
        gname = 'spec-0976-52413-0238'
        T_OIII = 14171.927530970173
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0987-52523-0394':
        gname = 'spec-0987-52523-0394'
        T_OIII = 11570.614681419422
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1158-52668-0120':
        gname = 'spec-1158-52668-0120'
        T_OIII = 11103.163057786347
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1163-52669-0465':
        gname = 'spec-1163-52669-0465'
        T_OIII = 10389.964108513297
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1233-52734-0335':
        gname = 'spec-1233-52734-0335'
        T_OIII = 12345.260727623503
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1116-52932-0089':
        gname = 'spec-1116-52932-0089'
        T_OIII = 16107.485846785732
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1005-52703-0268':
        gname = 'spec-1005-52703-0268'
        T_OIII = 14181.564284025484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1283-52762-0315':
        gname = 'spec-1283-52762-0315'
        T_OIII = 11148.540300786497
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-0999-52636-0517':
        gname = 'spec-0999-52636-0517'
        T_OIII = 10107.788382357448
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1279-52736-0147':
        gname = 'spec-1279-52736-0147'
        T_OIII = 11644.256719254627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1269-52937-0177':
        gname = 'spec-1269-52937-0177'
        T_OIII = 10938.354305551531
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1212-52703-0388':
        gname = 'spec-1212-52703-0388'
        T_OIII = 11670.67088060627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1212-52703-0606':
        gname = 'spec-1212-52703-0606'
        T_OIII = 12074.150216632315
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1073-52649-0409':
        gname = 'spec-1073-52649-0409'
        T_OIII = 13823.0656563437
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1073-52649-0419':
        gname = 'spec-1073-52649-0419'
        T_OIII = 14242.749410016355
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1104-52912-0511':
        gname = 'spec-1104-52912-0511'
        T_OIII = 12662.563276571509
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1021-52460-0446':
        gname = 'spec-1021-52460-0446'
        T_OIII = 11825.058909682957
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1217-52672-0271':
        gname = 'spec-1217-52672-0271'
        T_OIII = 10853.2225887174
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1209-52674-0598':
        gname = 'spec-1209-52674-0598'
        T_OIII = 11219.496014829632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1284-52736-0256':
        gname = 'spec-1284-52736-0256'
        T_OIII = 10480.302005849706
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1293-52765-0589':
        gname = 'spec-1293-52765-0589'
        T_OIII = 21173.98572406204
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1214-52731-0142':
        gname = 'spec-1214-52731-0142'
        T_OIII = 11857.255364936307
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1214-52731-0339':
        gname = 'spec-1214-52731-0339'
        T_OIII = 10319.57702029705
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1238-52761-0515':
        gname = 'spec-1238-52761-0515'
        T_OIII = 11528.742812258617
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1028-52884-0386':
        gname = 'spec-1028-52884-0386'
        T_OIII = 13300.731228795645
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1028-52884-0401':
        gname = 'spec-1028-52884-0401'
        T_OIII = 11713.058226051426
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1030-52914-0077':
        gname = 'spec-1030-52914-0077'
        T_OIII = 10968.136502632646
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1030-52914-0107':
        gname = 'spec-1030-52914-0107'
        T_OIII = 10051.87118805355
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1170-52756-0485':
        gname = 'spec-1170-52756-0485'
        T_OIII = 11022.947962474924
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1324-53088-0524':
        gname = 'spec-1324-53088-0524'
        T_OIII = 11260.244713697612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1163-52669-0284':
        gname = 'spec-1163-52669-0284'
        T_OIII = 12156.694748140775
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1171-52753-0376':
        gname = 'spec-1171-52753-0376'
        T_OIII = 11296.02121098684
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1346-52822-0185':
        gname = 'spec-1346-52822-0185'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1338-52765-0329':
        gname = 'spec-1338-52765-0329'
        T_OIII = 11219.496014829632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1313-52790-0219':
        gname = 'spec-1313-52790-0219'
        T_OIII = 17587.69419138251
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1164-52674-0418':
        gname = 'spec-1164-52674-0418'
        T_OIII = 12031.020146641786
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1053-52468-0249':
        gname = 'spec-1053-52468-0249'
        T_OIII = 15453.24321812379
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1355-52823-0098':
        gname = 'spec-1355-52823-0098'
        T_OIII = 10310.22817378569
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1363-53053-0599':
        gname = 'spec-1363-53053-0599'
        T_OIII = 11429.36758585671
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1365-53062-0338':
        gname = 'spec-1365-53062-0338'
        T_OIII = 10489.313318651464
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1279-52736-0546':
        gname = 'spec-1279-52736-0546'
        T_OIII = 11649.534764548354
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1074-52937-0186':
        gname = 'spec-1074-52937-0186'
        T_OIII = 10475.062623048356
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1199-52703-0600':
        gname = 'spec-1199-52703-0600'
        T_OIII = 12052.283381135889
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1267-52932-0384':
        gname = 'spec-1267-52932-0384'
        T_OIII = 10527.410027201355
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1331-52766-0554':
        gname = 'spec-1331-52766-0554'
        T_OIII = 13020.433563005136
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1213-52972-0088':
        gname = 'spec-1213-52972-0088'
        T_OIII = 10015.495364746293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1059-52618-0309':
        gname = 'spec-1059-52618-0309'
        T_OIII = 17318.784528456297
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1059-52618-0359':
        gname = 'spec-1059-52618-0359'
        T_OIII = 11927.31540013229
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1369-53089-0085':
        gname = 'spec-1369-53089-0085'
        T_OIII = 12418.20420184761
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1286-52725-0579':
        gname = 'spec-1286-52725-0579'
        T_OIII = 11660.098033431159
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1093-52591-0261':
        gname = 'spec-1093-52591-0261'
        T_OIII = 12345.260727623503
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1213-52972-0500':
        gname = 'spec-1213-52972-0500'
        T_OIII = 13079.368014936108
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1233-52734-0136':
        gname = 'spec-1233-52734-0136'
        T_OIII = 10456.091807953026
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1203-52669-0570':
        gname = 'spec-1203-52669-0570'
        T_OIII = 10218.960829093914
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-1073-52649-0188':
    #	 gname = 'spec-1073-52649-0188'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1074-52937-0573':
        gname = 'spec-1074-52937-0573'
        T_OIII = 11803.643189878701
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1306-52996-0230':
        gname = 'spec-1306-52996-0230'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1306-52996-0292':
        gname = 'spec-1306-52996-0292'
        T_OIII = 14204.075532297544
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1215-52725-0273':
        gname = 'spec-1215-52725-0273'
        T_OIII = 11660.098033431159
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1394-53108-0294':
        gname = 'spec-1394-53108-0294'
        T_OIII = 12222.792687463036
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1373-53063-0205':
        gname = 'spec-1373-53063-0205'
        T_OIII = 11644.256719254627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1399-53172-0299':
        gname = 'spec-1399-53172-0299'
        T_OIII = 15803.80771026763
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1299-52972-0486':
        gname = 'spec-1299-52972-0486'
        T_OIII = 10106.209000506708
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1286-52725-0150':
        gname = 'spec-1286-52725-0150'
        T_OIII = 12457.07581456919
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1333-52782-0172':
        gname = 'spec-1333-52782-0172'
        T_OIII = 13038.147083391987
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1349-52797-0175':
        gname = 'spec-1349-52797-0175'
        T_OIII = 12282.91435200761
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1158-52668-0062':
        gname = 'spec-1158-52668-0062'
        T_OIII = 13766.802451262882
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1319-52791-0443':
        gname = 'spec-1319-52791-0443'
        T_OIII = 10779.697143604779
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1424-52912-0189':
        gname = 'spec-1424-52912-0189'
        T_OIII = 11316.515936324256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1401-53144-0397':
        gname = 'spec-1401-53144-0397'
        T_OIII = 10666.674179366131
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1350-52786-0300':
        gname = 'spec-1350-52786-0300'
        T_OIII = 11570.614681419422
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1467-53115-0032':
        gname = 'spec-1467-53115-0032'
        T_OIII = 11580.201562660724
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1357-53034-0540':
        gname = 'spec-1357-53034-0540'
        T_OIII = 11008.488339349717
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1345-52814-0364':
        gname = 'spec-1345-52814-0364'
        T_OIII = 10864.252134685259
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1457-53116-0394':
        gname = 'spec-1457-53116-0394'
        T_OIII = 12497.240786234663
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1425-52913-0518':
        gname = 'spec-1425-52913-0518'
        T_OIII = 11316.515936324256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1304-52993-0210':
        gname = 'spec-1304-52993-0210'
        T_OIII = 10958.200104611962
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1314-53050-0214':
        gname = 'spec-1314-53050-0214'
        T_OIII = 12074.150216632315
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1268-52933-0318':
        gname = 'spec-1268-52933-0318'
        T_OIII = 11270.454992047864
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1275-52996-0232':
        gname = 'spec-1275-52996-0232'
        T_OIII = 10371.147410700316
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1277-52765-0031':
        gname = 'spec-1277-52765-0031'
        T_OIII = 9193.382658496088
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1162-52668-0224':
        gname = 'spec-1162-52668-0224'
        T_OIII = 11199.176996985096
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1313-52790-0408':
        gname = 'spec-1313-52790-0408'
        T_OIII = 10894.173044276355
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1356-53033-0248':
        gname = 'spec-1356-53033-0248'
        T_OIII = 10399.385255947842
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1370-53090-0254':
        gname = 'spec-1370-53090-0254'
        T_OIII = 11057.452127270644
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1377-53050-0362':
        gname = 'spec-1377-53050-0362'
        T_OIII = 12322.90289865434
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1461-53062-0068':
        gname = 'spec-1461-53062-0068'
        T_OIII = 10188.98343211608
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1464-53091-0442':
        gname = 'spec-1464-53091-0442'
        T_OIII = 13026.335394246644
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1467-53115-0579':
        gname = 'spec-1467-53115-0579'
        T_OIII = 11586.355810012077
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1455-53089-0569':
        gname = 'spec-1455-53089-0569'
        T_OIII = 15509.368694134177
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1322-52791-0470':
        gname = 'spec-1322-52791-0470'
        T_OIII = 10371.147410700316
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1326-52764-0028':
        gname = 'spec-1326-52764-0028'
        T_OIII = 9870.22396233294
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1376-53089-0637':
        gname = 'spec-1376-53089-0637'
        T_OIII = 10819.524782081646
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1310-53033-0508':
        gname = 'spec-1310-53033-0508'
        T_OIII = 11250.04368517661
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1223-52781-0287':
        gname = 'spec-1223-52781-0287'
        T_OIII = 12406.954139485417
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1385-53108-0015':
        gname = 'spec-1385-53108-0015'
        T_OIII = 11012.961908937923
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1385-53108-0559':
        gname = 'spec-1385-53108-0559'
        T_OIII = 10671.008844375976
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1466-53083-0092':
        gname = 'spec-1466-53083-0092'
        T_OIII = 11132.872065200108
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1316-52790-0581':
        gname = 'spec-1316-52790-0581'
        T_OIII = 13246.594050585314
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1590-52974-0628':
        gname = 'spec-1590-52974-0628'
        T_OIII = 14496.709119002047
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1595-52999-0624':
        gname = 'spec-1595-52999-0624'
        T_OIII = 10249.987112748755
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1487-52964-0510':
        gname = 'spec-1487-52964-0510'
        T_OIII = 13115.184560905886
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1362-53050-0617':
        gname = 'spec-1362-53050-0617'
        T_OIII = 12035.909247141666
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1342-52793-0537':
        gname = 'spec-1342-52793-0537'
        T_OIII = 9400.249673124554
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1607-53083-0205':
        gname = 'spec-1607-53083-0205'
        T_OIII = 11649.534764548354
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1405-52826-0395':
        gname = 'spec-1405-52826-0395'
        T_OIII = 10854.070618869435
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1408-52822-0221':
        gname = 'spec-1408-52822-0221'
        T_OIII = 10701.06641335077
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1350-52786-0090':
        gname = 'spec-1350-52786-0090'
        T_OIII = 10389.964108513297
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1372-53062-0517':
        gname = 'spec-1372-53062-0517'
        T_OIII = 10002.98252386497
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1215-52725-0629':
        gname = 'spec-1215-52725-0629'
        T_OIII = 12452.015632214747
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1624-53386-0636':
        gname = 'spec-1624-53386-0636'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1292-52736-0544':
        gname = 'spec-1292-52736-0544'
        T_OIII = 22663.35081901792
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1301-52976-0356':
        gname = 'spec-1301-52976-0356'
        T_OIII = 10958.200104611962
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1387-53118-0458':
        gname = 'spec-1387-53118-0458'
        T_OIII = 9336.129745079415
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1576-53496-0446':
        gname = 'spec-1576-53496-0446'
        T_OIII = 9764.523887739684
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1428-52998-0568':
        gname = 'spec-1428-52998-0568'
        T_OIII = 10484.560933949853
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1430-53002-0389':
        gname = 'spec-1430-53002-0389'
        T_OIII = 10724.838224565667
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1442-53050-0548':
        gname = 'spec-1442-53050-0548'
        T_OIII = 10943.825416081932
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1357-53034-0202':
        gname = 'spec-1357-53034-0202'
        T_OIII = 11435.084289814917
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1457-53116-0533':
        gname = 'spec-1457-53116-0533'
        T_OIII = 10633.721753606136
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1596-52998-0479':
        gname = 'spec-1596-52998-0479'
        T_OIII = 10296.703478374322
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1597-52999-0050':
        gname = 'spec-1597-52999-0050'
        T_OIII = 11636.980752400705
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1455-53089-0287':
        gname = 'spec-1455-53089-0287'
        T_OIII = 12029.892174745002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1417-53141-0306':
        gname = 'spec-1417-53141-0306'
        T_OIII = 11103.163057786347
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1651-53442-0255':
        gname = 'spec-1651-53442-0255'
        T_OIII = 14181.564284025484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1450-53120-0201':
        gname = 'spec-1450-53120-0201'
        T_OIII = 10600.374161445068
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1442-53050-0396':
        gname = 'spec-1442-53050-0396'
        T_OIII = 9598.863075529807
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1682-53173-0277':
        gname = 'spec-1682-53173-0277'
        T_OIII = 10461.321701993696
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1492-52932-0010':
        gname = 'spec-1492-52932-0010'
        T_OIII = 14470.454919090893
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1573-53226-0512':
        gname = 'spec-1573-53226-0512'
        T_OIII = 16424.44071994387
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1584-52943-0372':
        gname = 'spec-1584-52943-0372'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1620-53137-0470':
        gname = 'spec-1620-53137-0470'
        T_OIII = 10647.356347022687
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1464-53091-0232':
        gname = 'spec-1464-53091-0232'
        T_OIII = 10022.854104760241
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1310-53033-0488':
        gname = 'spec-1310-53033-0488'
        T_OIII = 12748.931250436992
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1689-53177-0521':
        gname = 'spec-1689-53177-0521'
        T_OIII = 12085.098507516668
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1363-53053-0411':
        gname = 'spec-1363-53053-0411'
        T_OIII = 11959.790272308908
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1372-53062-0402':
        gname = 'spec-1372-53062-0402'
        T_OIII = 16195.31796006577
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1615-53166-0206':
        gname = 'spec-1615-53166-0206'
        T_OIII = 11326.777238862745
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1495-52944-0627':
        gname = 'spec-1495-52944-0627'
        T_OIII = 13336.945542915813
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1579-53473-0103':
        gname = 'spec-1579-53473-0103'
        T_OIII = 10184.844570400828
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1434-53053-0021':
        gname = 'spec-1434-53053-0021'
        T_OIII = 14909.794135786862
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1314-53050-0582':
        gname = 'spec-1314-53050-0582'
        T_OIII = 12870.836989148727
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1325-52762-0412':
        gname = 'spec-1325-52762-0412'
        T_OIII = 11575.85934665115
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1705-53848-0450':
        gname = 'spec-1705-53848-0450'
        T_OIII = 10740.603535660415
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1602-53117-0460':
        gname = 'spec-1602-53117-0460'
        T_OIII = 10508.344408574927
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1676-53147-0294':
        gname = 'spec-1676-53147-0294'
        T_OIII = 10766.229548623747
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1685-53463-0234':
        gname = 'spec-1685-53463-0234'
        T_OIII = 9693.97982653396
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1685-53463-0440':
        gname = 'spec-1685-53463-0440'
        T_OIII = 12877.274682591395
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1622-53385-0193':
        gname = 'spec-1622-53385-0193'
        T_OIII = 10498.8245514363
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1458-53119-0543':
        gname = 'spec-1458-53119-0543'
        T_OIII = 10815.129788489168
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1462-53112-0184':
        gname = 'spec-1462-53112-0184'
        T_OIII = 11022.947962474924
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1502-53741-0087':
        gname = 'spec-1502-53741-0087'
        T_OIII = 13361.143181237905
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1712-53531-0569':
        gname = 'spec-1712-53531-0569'
        T_OIII = 10690.870764829013
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1721-53857-0206':
        gname = 'spec-1721-53857-0206'
        T_OIII = 10775.991876465538
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1578-53496-0438':
        gname = 'spec-1578-53496-0438'
        T_OIII = 11032.943070911602
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1627-53473-0426':
        gname = 'spec-1627-53473-0426'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1695-53473-0627':
        gname = 'spec-1695-53473-0627'
        T_OIII = 9857.892616013474
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1576-53496-0378':
        gname = 'spec-1576-53496-0378'
        T_OIII = 11681.253314758873
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1623-53089-0454':
        gname = 'spec-1623-53089-0454'
        T_OIII = 11378.223487652196
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1344-52792-0396':
        gname = 'spec-1344-52792-0396'
        T_OIII = 12457.07581456919
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1351-52790-0437':
        gname = 'spec-1351-52790-0437'
        T_OIII = 15182.526525228563
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1351-52790-0474':
        gname = 'spec-1351-52790-0474'
        T_OIII = 13464.47389542078
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1374-53083-0560':
        gname = 'spec-1374-53083-0560'
        T_OIII = 12150.997043837559
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1741-53052-0219':
        gname = 'spec-1741-53052-0219'
        T_OIII = 12063.211844167205
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1588-52965-0598':
        gname = 'spec-1588-52965-0598'
        T_OIII = 12644.767222631524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1707-53885-0025':
        gname = 'spec-1707-53885-0025'
        T_OIII = 9952.154139421937
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1708-53503-0492':
        gname = 'spec-1708-53503-0492'
        T_OIII = 12800.434446355654
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1708-53503-0512':
        gname = 'spec-1708-53503-0512'
        T_OIII = 10010.957652361838
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1724-53859-0322':
        gname = 'spec-1724-53859-0322'
        T_OIII = 11103.163057786347
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1388-53119-0189':
        gname = 'spec-1388-53119-0189'
        T_OIII = 10556.073333529317
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1744-53055-0450':
        gname = 'spec-1744-53055-0450'
        T_OIII = 11083.054724220938
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1597-52999-0358':
        gname = 'spec-1597-52999-0358'
        T_OIII = 10456.091807953026
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1605-53062-0061':
        gname = 'spec-1605-53062-0061'
        T_OIII = 11549.659771691455
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1733-53047-0240':
        gname = 'spec-1733-53047-0240'
        T_OIII = 12870.836989148727
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1609-53142-0238':
        gname = 'spec-1609-53142-0238'
        T_OIII = 12362.05571578742
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1418-53142-0442':
        gname = 'spec-1418-53142-0442'
        T_OIII = 11954.371659803212
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1771-53498-0429':
        gname = 'spec-1771-53498-0429'
        T_OIII = 11148.540300786497
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1682-53173-0003':
        gname = 'spec-1682-53173-0003'
        T_OIII = 16641.71503582341
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1761-53376-0601':
        gname = 'spec-1761-53376-0601'
        T_OIII = 11398.867355156632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1722-53852-0358':
        gname = 'spec-1722-53852-0358'
        T_OIII = 10508.344408574927
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1431-52992-0451':
        gname = 'spec-1431-52992-0451'
        T_OIII = 15439.243613424378
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1433-53035-0376':
        gname = 'spec-1433-53035-0376'
        T_OIII = 12162.015015933623
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1782-53299-0624':
        gname = 'spec-1782-53299-0624'
        T_OIII = 12345.260727623503
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1801-54156-0583':
        gname = 'spec-1801-54156-0583'
        T_OIII = 11158.649290649884
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1803-54152-0448':
        gname = 'spec-1803-54152-0448'
        T_OIII = 10175.617782476917
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1814-54555-0293':
        gname = 'spec-1814-54555-0293'
        T_OIII = 11798.295322805352
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1814-54555-0395':
        gname = 'spec-1814-54555-0395'
        T_OIII = 11209.331901904156
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1644-53144-0564':
        gname = 'spec-1644-53144-0564'
        T_OIII = 14523.010952731023
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1693-53446-0541':
        gname = 'spec-1693-53446-0541'
        T_OIII = 11378.223487652196
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1708-53503-0263':
        gname = 'spec-1708-53503-0263'
        T_OIII = 11032.943070911602
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1616-53169-0205':
        gname = 'spec-1616-53169-0205'
        T_OIII = 11575.85934665115
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1626-53472-0539':
        gname = 'spec-1626-53472-0539'
        T_OIII = 11012.961908937923
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1795-54507-0315':
        gname = 'spec-1795-54507-0315'
        T_OIII = 11959.790272308908
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1799-53556-0055':
        gname = 'spec-1799-53556-0055'
        T_OIII = 10952.721807790953
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1800-53884-0520':
        gname = 'spec-1800-53884-0520'
        T_OIII = 11445.453104749757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1647-53531-0034':
        gname = 'spec-1647-53531-0034'
        T_OIII = 12300.585560730844
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1649-53149-0344':
        gname = 'spec-1649-53149-0344'
        T_OIII = 11644.256719254627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1714-53521-0042':
        gname = 'spec-1714-53521-0042'
        T_OIII = 11713.058226051426
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1827-53531-0228':
        gname = 'spec-1827-53531-0228'
        T_OIII = 13513.376251770394
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1840-53472-0032':
        gname = 'spec-1840-53472-0032'
        T_OIII = 11073.014217978294
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1745-53061-0463':
        gname = 'spec-1745-53061-0463'
        T_OIII = 14181.564284025484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1745-53061-0475':
        gname = 'spec-1745-53061-0475'
        T_OIII = 13716.982993362524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1752-53379-0530':
        gname = 'spec-1752-53379-0530'
        T_OIII = 10835.259960873123
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1753-53383-0146':
        gname = 'spec-1753-53383-0146'
        T_OIII = 11721.022975334794
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1666-52991-0310':
        gname = 'spec-1666-52991-0310'
        T_OIII = 12052.283381135889
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1713-53827-0432':
        gname = 'spec-1713-53827-0432'
        T_OIII = 12691.287393325083
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1718-53850-0484':
        gname = 'spec-1718-53850-0484'
        T_OIII = 13210.62507129942
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1718-53850-0542':
        gname = 'spec-1718-53850-0542'
        T_OIII = 11702.4469787766
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1570-53149-0146':
        gname = 'spec-1570-53149-0146'
        T_OIII = 10705.415054489516
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1580-53145-0419':
        gname = 'spec-1580-53145-0419'
        T_OIII = 12917.583283274278
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1647-53531-0344':
        gname = 'spec-1647-53531-0344'
        T_OIII = 12492.164288487753
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1455-53089-0076':
        gname = 'spec-1455-53089-0076'
        T_OIII = 11660.098033431159
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1719-53876-0008':
        gname = 'spec-1719-53876-0008'
        T_OIII = 11062.982807762703
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1719-53876-0196':
        gname = 'spec-1719-53876-0196'
        T_OIII = 10338.300152334337
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1679-53149-0534':
        gname = 'spec-1679-53149-0534'
        T_OIII = 11782.26625491757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1734-53034-0490':
        gname = 'spec-1734-53034-0490'
        T_OIII = 11709.032138011396
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1650-53174-0305':
        gname = 'spec-1650-53174-0305'
        T_OIII = 9831.125106270674
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1461-53062-0582':
        gname = 'spec-1461-53062-0582'
        T_OIII = 12067.359715003811
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1820-54208-0585':
        gname = 'spec-1820-54208-0585'
        T_OIII = 11435.084289814917
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1702-53144-0062':
        gname = 'spec-1702-53144-0062'
        T_OIII = 11857.255364936307
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1715-54212-0211':
        gname = 'spec-1715-54212-0211'
        T_OIII = 10676.346234921906
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1771-53498-0627':
        gname = 'spec-1771-53498-0627'
        T_OIII = 12639.6307982468
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1801-54156-0525':
        gname = 'spec-1801-54156-0525'
        T_OIII = 13459.004498410439
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1747-53075-0182':
        gname = 'spec-1747-53075-0182'
        T_OIII = 10446.619288175401
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1749-53357-0173':
        gname = 'spec-1749-53357-0173'
        T_OIII = 12317.897218610982
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1870-53383-0072':
        gname = 'spec-1870-53383-0072'
        T_OIII = 12128.991035239995
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1724-53859-0627':
        gname = 'spec-1724-53859-0627'
        T_OIII = 15775.186295326017
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1725-54266-0068':
        gname = 'spec-1725-54266-0068'
        T_OIII = 15629.313082271625
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1571-53174-0155':
        gname = 'spec-1571-53174-0155'
        T_OIII = 10319.57702029705
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1817-53851-0358':
        gname = 'spec-1817-53851-0358'
        T_OIII = 10006.42199587506
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1819-54540-0638':
        gname = 'spec-1819-54540-0638'
        T_OIII = 12435.098425241405
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1792-54270-0119':
        gname = 'spec-1792-54270-0119'
        T_OIII = 12830.072771726329
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1697-53142-0048':
        gname = 'spec-1697-53142-0048'
        T_OIII = 11239.851898105111
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1706-53442-0099':
        gname = 'spec-1706-53442-0099'
        T_OIII = 13452.275982700403
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1733-53047-0326':
        gname = 'spec-1733-53047-0326'
        T_OIII = 13300.731228795645
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1733-53047-0528':
        gname = 'spec-1733-53047-0528'
        T_OIII = 10705.415054489516
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1778-53883-0143':
        gname = 'spec-1778-53883-0143'
        T_OIII = 12737.381571490869
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1926-53317-0369':
        gname = 'spec-1926-53317-0369'
        T_OIII = 10233.982509902751
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1823-53886-0464':
        gname = 'spec-1823-53886-0464'
        T_OIII = 12311.739172898178
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1839-53471-0483':
        gname = 'spec-1839-53471-0483'
        T_OIII = 13986.8987732756
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1763-53463-0225':
        gname = 'spec-1763-53463-0225'
        T_OIII = 10259.441653250507
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1774-53759-0060':
        gname = 'spec-1774-53759-0060'
        T_OIII = 10213.053938304782
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1806-53559-0414':
        gname = 'spec-1806-53559-0414'
        T_OIII = 12979.195546771743
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1601-53115-0168':
        gname = 'spec-1601-53115-0168'
        T_OIII = 12173.042978626236
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1610-53144-0379':
        gname = 'spec-1610-53144-0379'
        T_OIII = 11158.649290649884
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1758-53084-0338':
        gname = 'spec-1758-53084-0338'
        T_OIII = 11409.20333022049
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1946-53432-0528':
        gname = 'spec-1946-53432-0528'
        T_OIII = 12685.933833791785
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1953-53358-0618':
        gname = 'spec-1953-53358-0618'
        T_OIII = 10613.634286097713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1785-54439-0341':
        gname = 'spec-1785-54439-0341'
        T_OIII = 13324.863162352554
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1594-52992-0343':
        gname = 'spec-1594-52992-0343'
        T_OIII = 11306.263929839473
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1622-53385-0403':
        gname = 'spec-1622-53385-0403'
        T_OIII = 11487.022469483838
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1983-53442-0308':
        gname = 'spec-1983-53442-0308'
        T_OIII = 12334.07674716982
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1872-53386-0526':
        gname = 'spec-1872-53386-0526'
        T_OIII = 10730.202539316151
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1799-53556-0392':
        gname = 'spec-1799-53556-0392'
        T_OIII = 11771.592309797865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1850-53786-0209':
        gname = 'spec-1850-53786-0209'
        T_OIII = 13525.629567438933
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1690-53475-0108':
        gname = 'spec-1690-53475-0108'
        T_OIII = 11660.098033431159
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1702-53144-0234':
        gname = 'spec-1702-53144-0234'
        T_OIII = 9533.83532230403
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-1775-53847-0305':
    #	 gname = 'spec-1775-53847-0305'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1938-53379-0406':
        gname = 'spec-1938-53379-0406'
        T_OIII = 11235.286163651193
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1802-53885-0346':
        gname = 'spec-1802-53885-0346'
        T_OIII = 10993.01693375783
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1824-53491-0063':
        gname = 'spec-1824-53491-0063'
        T_OIII = 11164.928445988717
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1836-54567-0353':
        gname = 'spec-1836-54567-0353'
        T_OIII = 9799.988194065465
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1836-54567-0390':
        gname = 'spec-1836-54567-0390'
        T_OIII = 11363.65297399599
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1834-54562-0290':
        gname = 'spec-1834-54562-0290'
        T_OIII = 10291.555881379016
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1834-54562-0306':
        gname = 'spec-1834-54562-0306'
        T_OIII = 10849.49204893719
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2037-53446-0293':
        gname = 'spec-2037-53446-0293'
        T_OIII = 11113.230901627261
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1936-53330-0020':
        gname = 'spec-1936-53330-0020'
        T_OIII = 13367.199447480294
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1937-53388-0488':
        gname = 'spec-1937-53388-0488'
        T_OIII = 12272.745717156748
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-1950-53436-0189':
    #	 gname = 'spec-1950-53436-0189'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1843-53816-0039':
        gname = 'spec-1843-53816-0039'
        T_OIII = 10938.354305551531
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1848-54180-0463':
        gname = 'spec-1848-54180-0463'
        T_OIII = 10831.027826489903
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1996-53436-0475':
        gname = 'spec-1996-53436-0475'
        T_OIII = 15244.575618077715
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2003-53442-0450':
        gname = 'spec-2003-53442-0450'
        T_OIII = 14923.313659752961
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1842-53501-0210':
        gname = 'spec-1842-53501-0210'
        T_OIII = 13385.384722093637
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1678-53433-0425':
        gname = 'spec-1678-53433-0425'
        T_OIII = 17217.053694730042
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1864-53313-0228':
        gname = 'spec-1864-53313-0228'
        T_OIII = 11194.102994920871
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1957-53415-0258':
        gname = 'spec-1957-53415-0258'
        T_OIII = 10523.791477819543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1749-53357-0026':
        gname = 'spec-1749-53357-0026'
        T_OIII = 11018.470336457072
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1763-53463-0094':
        gname = 'spec-1763-53463-0094'
        T_OIII = 11269.926647222295
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2031-53848-0009':
        gname = 'spec-2031-53848-0009'
        T_OIII = 19286.780960935914
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2032-53815-0531':
        gname = 'spec-2032-53815-0531'
        T_OIII = 11894.928708207906
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1847-54176-0496':
        gname = 'spec-1847-54176-0496'
        T_OIII = 10613.634286097713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1851-53524-0225':
        gname = 'spec-1851-53524-0225'
        T_OIII = 9906.072166518683
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1710-53504-0632':
        gname = 'spec-1710-53504-0632'
        T_OIII = 10338.300152334337
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1876-54464-0379':
        gname = 'spec-1876-54464-0379'
        T_OIII = 10666.674179366131
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1974-53430-0155':
        gname = 'spec-1974-53430-0155'
        T_OIII = 14346.395322644312
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1929-53349-0593':
        gname = 'spec-1929-53349-0593'
        T_OIII = 11528.742812258617
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1785-54439-0201':
        gname = 'spec-1785-54439-0201'
        T_OIII = 15134.440778138678
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1684-53239-0484':
        gname = 'spec-1684-53239-0484'
        T_OIII = 12058.31165309344
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1881-53261-0528':
        gname = 'spec-1881-53261-0528'
        T_OIII = 11164.928445988717
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1706-53442-0015':
        gname = 'spec-1706-53442-0015'
        T_OIII = 11414.018102661712
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1744-53055-0385':
        gname = 'spec-1744-53055-0385'
        T_OIII = 14984.303381941287
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1922-53315-0588':
        gname = 'spec-1922-53315-0588'
        T_OIII = 10226.468760245178
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1924-53330-0612':
        gname = 'spec-1924-53330-0612'
        T_OIII = 14229.846431738068
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1927-53321-0006':
        gname = 'spec-1927-53321-0006'
        T_OIII = 10932.202575500773
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2098-53460-0291':
        gname = 'spec-2098-53460-0291'
        T_OIII = 11507.86373439115
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1949-53433-0592':
        gname = 'spec-1949-53433-0592'
        T_OIII = 14082.299687081457
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1748-53112-0217':
        gname = 'spec-1748-53112-0217'
        T_OIII = 10347.674453240526
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1748-53112-0218':
        gname = 'spec-1748-53112-0218'
        T_OIII = 11798.295322805352
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2100-53713-0211':
        gname = 'spec-2100-53713-0211'
        T_OIII = 10175.617782476917
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1942-53415-0577':
        gname = 'spec-1942-53415-0577'
        T_OIII = 10695.716664527605
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1767-53436-0514':
        gname = 'spec-1767-53436-0514'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1774-53759-0638':
        gname = 'spec-1774-53759-0638'
        T_OIII = 12882.507679404664
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1957-53415-0535':
        gname = 'spec-1957-53415-0535'
        T_OIII = 12184.080940974438
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2112-53534-0555':
        gname = 'spec-2112-53534-0555'
        T_OIII = 17046.256263241805
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1875-54453-0549':
        gname = 'spec-1875-54453-0549'
        T_OIII = 11285.78777135236
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1754-53385-0131':
        gname = 'spec-1754-53385-0131'
        T_OIII = 12457.659827760282
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1959-53440-0545':
        gname = 'spec-1959-53440-0545'
        T_OIII = 13234.593522430534
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2123-53793-0288':
        gname = 'spec-2123-53793-0288'
        T_OIII = 11073.014217978294
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2132-53493-0241':
        gname = 'spec-2132-53493-0241'
        T_OIII = 16499.040568222572
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1981-53463-0426':
        gname = 'spec-1981-53463-0426'
        T_OIII = 13061.80260230364
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1981-53463-0438':
        gname = 'spec-1981-53463-0438'
        T_OIII = 10582.995541891527
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1786-54450-0448':
        gname = 'spec-1786-54450-0448'
        T_OIII = 11518.29854241994
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2132-53493-0445':
        gname = 'spec-2132-53493-0445'
        T_OIII = 10839.663134855515
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2154-54539-0217':
        gname = 'spec-2154-54539-0217'
        T_OIII = 11590.339716453223
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1993-53762-0564':
        gname = 'spec-1993-53762-0564'
        T_OIII = 10628.073500152708
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1999-53503-0119':
        gname = 'spec-1999-53503-0119'
        T_OIII = 10681.519379107072
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2029-53819-0505':
        gname = 'spec-2029-53819-0505'
        T_OIII = 10399.385255947842
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1997-53442-0006':
        gname = 'spec-1997-53442-0006'
        T_OIII = 15675.418925537853
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1795-54507-0461':
        gname = 'spec-1795-54507-0461'
        T_OIII = 10604.019043372815
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2156-54525-0081':
        gname = 'spec-2156-54525-0081'
        T_OIII = 12139.989053287203
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2156-54525-0140':
        gname = 'spec-2156-54525-0140'
        T_OIII = 13642.591641052915
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2165-53917-0022':
        gname = 'spec-2165-53917-0022'
        T_OIII = 11239.851898105111
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2010-53495-0338':
        gname = 'spec-2010-53495-0338'
        T_OIII = 23468.130492255994
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2010-53495-0127':
        gname = 'spec-2010-53495-0127'
        T_OIII = 11158.649290649884
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1809-53792-0586':
        gname = 'spec-1809-53792-0586'
        T_OIII = 13324.863162352554
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1824-53491-0320':
        gname = 'spec-1824-53491-0320'
        T_OIII = 10498.8245514363
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2028-53818-0494':
        gname = 'spec-2028-53818-0494'
        T_OIII = 10432.915678364892
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1864-53313-0359':
        gname = 'spec-1864-53313-0359'
        T_OIII = 22135.533805400744
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-2114-53848-0381':
    #	 gname = 'spec-2114-53848-0381'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-2136-53494-0081':
    #	 gname = 'spec-2136-53494-0081'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1935-53387-0204':
        gname = 'spec-1935-53387-0204'
        T_OIII = 11042.947242458555
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2016-53799-0476':
        gname = 'spec-2016-53799-0476'
        T_OIII = 12795.23478846811
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2167-53889-0429':
        gname = 'spec-2167-53889-0429'
        T_OIII = 10993.01693375783
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2173-53874-0127':
        gname = 'spec-2173-53874-0127'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2101-53858-0387':
        gname = 'spec-2101-53858-0387'
        T_OIII = 10993.01693375783
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1854-53566-0551':
        gname = 'spec-1854-53566-0551'
        T_OIII = 10948.272708303248
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2095-53474-0074':
        gname = 'spec-2095-53474-0074'
        T_OIII = 13791.779992446884
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2231-53816-0443':
        gname = 'spec-2231-53816-0443'
        T_OIII = 11654.997328988069
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2166-54232-0594':
        gname = 'spec-2166-54232-0594'
        T_OIII = 11954.371659803212
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2173-53874-0538':
        gname = 'spec-2173-53874-0538'
        T_OIII = 9848.962024690649
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2117-54115-0484':
        gname = 'spec-2117-54115-0484'
        T_OIII = 10724.838224565667
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2130-53881-0238':
        gname = 'spec-2130-53881-0238'
        T_OIII = 10517.872897887539
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1986-53475-0260':
        gname = 'spec-1986-53475-0260'
        T_OIII = 11194.102994920871
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1868-53318-0353':
        gname = 'spec-1868-53318-0353'
        T_OIII = 10613.634286097713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2107-53786-0443':
        gname = 'spec-2107-53786-0443'
        T_OIII = 10637.710554366235
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1995-53415-0191':
        gname = 'spec-1995-53415-0191'
        T_OIII = 10019.721981598661
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2098-53460-0503':
        gname = 'spec-2098-53460-0503'
        T_OIII = 10575.22554824027
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2088-53493-0585':
        gname = 'spec-2088-53493-0585'
        T_OIII = 11857.255364936307
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2090-53463-0559':
        gname = 'spec-2090-53463-0559'
        T_OIII = 12941.020059504513
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2092-53460-0557':
        gname = 'spec-2092-53460-0557'
        T_OIII = 12173.042978626236
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2155-53820-0479':
        gname = 'spec-2155-53820-0479'
        T_OIII = 13716.982993362524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2088-53493-0353':
        gname = 'spec-2088-53493-0353'
        T_OIII = 12014.111667781179
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2108-53473-0077':
        gname = 'spec-2108-53473-0077'
        T_OIII = 9610.870398947525
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1932-53350-0565':
        gname = 'spec-1932-53350-0565'
        T_OIII = 10446.619288175401
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2138-53757-0152':
        gname = 'spec-2138-53757-0152'
        T_OIII = 16387.26737903439
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2148-54526-0555':
        gname = 'spec-2148-54526-0555'
        T_OIII = 13860.702146636713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2270-53714-0470':
        gname = 'spec-2270-53714-0470'
        T_OIII = 11042.947242458555
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2283-53729-0182':
        gname = 'spec-2283-53729-0182'
        T_OIII = 12435.098425241405
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2110-53467-0466':
        gname = 'spec-2110-53467-0466'
        T_OIII = 12485.919122013713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2110-53467-0496':
        gname = 'spec-2110-53467-0496'
        T_OIII = 11260.244713697612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2036-53446-0598':
        gname = 'spec-2036-53446-0598'
        T_OIII = 10709.765462805386
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2103-53467-0109':
        gname = 'spec-2103-53467-0109'
        T_OIII = 12719.480367402108
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2106-53714-0487':
        gname = 'spec-2106-53714-0487'
        T_OIII = 11868.006985422713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2115-53535-0414':
        gname = 'spec-2115-53535-0414'
        T_OIII = 12166.957349415672
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2150-54510-0397':
        gname = 'spec-2150-54510-0397'
        T_OIII = 12256.072064829234
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2154-54539-0318':
        gname = 'spec-2154-54539-0318'
        T_OIII = 10888.896910241861
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2154-54539-0338':
        gname = 'spec-2154-54539-0338'
        T_OIII = 10724.838224565667
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2240-53823-0220':
        gname = 'spec-2240-53823-0220'
        T_OIII = 12020.120847135122
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2123-53793-0003':
        gname = 'spec-2123-53793-0003'
        T_OIII = 23951.62293663108
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2126-53794-0436':
        gname = 'spec-2126-53794-0436'
        T_OIII = 12440.734938958238
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2196-53534-0331':
        gname = 'spec-2196-53534-0331'
        T_OIII = 10884.30365220158
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2153-54212-0357':
        gname = 'spec-2153-54212-0357'
        T_OIII = 10903.710602618274
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1927-53321-0560':
        gname = 'spec-1927-53321-0560'
        T_OIII = 11586.355810012077
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1935-53387-0141':
        gname = 'spec-1935-53387-0141'
        T_OIII = 10083.808392923831
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2137-54206-0012':
        gname = 'spec-2137-54206-0012'
        T_OIII = 12107.024880525423
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1974-53430-0425':
        gname = 'spec-1974-53430-0425'
        T_OIII = 13061.80260230364
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2143-54184-0161':
        gname = 'spec-2143-54184-0161'
        T_OIII = 11916.510050189228
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2171-53557-0183':
        gname = 'spec-2171-53557-0183'
        T_OIII = 9782.239969481387
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2171-53557-0324':
        gname = 'spec-2171-53557-0324'
        T_OIII = 14216.955142694998
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2147-53491-0514':
        gname = 'spec-2147-53491-0514'
        T_OIII = 11052.96048533375
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2154-54539-0639':
        gname = 'spec-2154-54539-0639'
        T_OIII = 10908.652977199827
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1994-53845-0417':
        gname = 'spec-1994-53845-0417'
        T_OIII = 10958.200104611962
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2008-53473-0505':
        gname = 'spec-2008-53473-0505'
        T_OIII = 12162.015015933623
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2167-53889-0210':
        gname = 'spec-2167-53889-0210'
        T_OIII = 12139.989053287203
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2169-53556-0041':
        gname = 'spec-2169-53556-0041'
        T_OIII = 12882.507679404664
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2268-53682-0518':
        gname = 'spec-2268-53682-0518'
        T_OIII = 12760.491402107278
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2152-53874-0110':
        gname = 'spec-2152-53874-0110'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2152-53874-0157':
        gname = 'spec-2152-53874-0157'
        T_OIII = 12870.836989148727
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1982-53436-0532':
        gname = 'spec-1982-53436-0532'
        T_OIII = 12025.005518437003
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1987-53765-0573':
        gname = 'spec-1987-53765-0573'
        T_OIII = 13312.791727617876
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2163-53823-0546':
        gname = 'spec-2163-53823-0546'
        T_OIII = 13067.110586568806
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2018-53800-0096':
        gname = 'spec-2018-53800-0096'
        T_OIII = 12616.739851683944
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2018-53800-0202':
        gname = 'spec-2018-53800-0202'
        T_OIII = 12074.150216632315
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2218-53816-0522':
        gname = 'spec-2218-53816-0522'
        T_OIII = 12435.098425241405
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2226-53819-0157':
        gname = 'spec-2226-53819-0157'
        T_OIII = 13723.20055273777
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2230-53799-0014':
        gname = 'spec-2230-53799-0014'
        T_OIII = 10913.085976236822
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2101-53858-0304':
        gname = 'spec-2101-53858-0304'
        T_OIII = 13409.670245136727
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2214-53794-0230':
        gname = 'spec-2214-53794-0230'
        T_OIII = 10114.740595279956
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-1995-53415-0045':
        gname = 'spec-1995-53415-0045'
        T_OIII = 11022.947962474924
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2241-54169-0451':
        gname = 'spec-2241-54169-0451'
        T_OIII = 11560.132478468087
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2246-53767-0295':
        gname = 'spec-2246-53767-0295'
        T_OIII = 10356.733569150587
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2364-53737-0324':
        gname = 'spec-2364-53737-0324'
        T_OIII = 15615.153970052514
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2289-53708-0204':
        gname = 'spec-2289-53708-0204'
        T_OIII = 12639.6307982468
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2187-54270-0155':
        gname = 'spec-2187-54270-0155'
        T_OIII = 12702.795276164543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2234-53823-0301':
        gname = 'spec-2234-53823-0301'
        T_OIII = 10628.073500152708
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2009-53904-0640':
        gname = 'spec-2009-53904-0640'
        T_OIII = 11388.540743794836
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2292-53713-0389':
        gname = 'spec-2292-53713-0389'
        T_OIII = 11938.130547875764
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2135-53827-0450':
        gname = 'spec-2135-53827-0450'
        T_OIII = 11644.256719254627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2245-54208-0594':
        gname = 'spec-2245-54208-0594'
        T_OIII = 10584.81468233188
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2025-53431-0571':
        gname = 'spec-2025-53431-0571'
        T_OIII = 11026.738144142759
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2353-53794-0526':
        gname = 'spec-2353-53794-0526'
        T_OIII = 17429.00943121612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2144-53770-0359':
        gname = 'spec-2144-53770-0359'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2144-53770-0550':
        gname = 'spec-2144-53770-0550'
        T_OIII = 13061.80260230364
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2148-54526-0212':
        gname = 'spec-2148-54526-0212'
        T_OIII = 13501.134036786234
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2034-53466-0618':
        gname = 'spec-2034-53466-0618'
        T_OIII = 13258.605460273811
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2237-53828-0259':
        gname = 'spec-2237-53828-0259'
        T_OIII = 9827.745888439898
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2237-53828-0382':
        gname = 'spec-2237-53828-0382'
        T_OIII = 11032.943070911602
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2093-53818-0320':
        gname = 'spec-2093-53818-0320'
        T_OIII = 10352.364790783382
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2224-53815-0337':
        gname = 'spec-2224-53815-0337'
        T_OIII = 10556.073333529317
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2347-53757-0501':
        gname = 'spec-2347-53757-0501'
        T_OIII = 11788.159470585906
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-2350-53765-0146':
    #	 gname = 'spec-2350-53765-0146'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2227-53820-0445':
        gname = 'spec-2227-53820-0445'
        T_OIII = 14815.499785332746
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2109-53468-0261':
        gname = 'spec-2109-53468-0261'
        T_OIII = 11123.307874531605
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-2267-53713-0082':
    #	 gname = 'spec-2267-53713-0082'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2275-53709-0316':
        gname = 'spec-2275-53709-0316'
        T_OIII = 14701.800862843473
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2275-53709-0349':
        gname = 'spec-2275-53709-0349'
        T_OIII = 14346.395322644312
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2280-53680-0176':
        gname = 'spec-2280-53680-0176'
        T_OIII = 12644.767222631524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2125-53795-0471':
        gname = 'spec-2125-53795-0471'
        T_OIII = 13679.73674913785
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2429-53799-0069':
        gname = 'spec-2429-53799-0069'
        T_OIII = 12009.23142166479
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2371-53762-0356':
        gname = 'spec-2371-53762-0356'
        T_OIII = 10938.354305551531
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2356-53786-0172':
        gname = 'spec-2356-53786-0172'
        T_OIII = 11239.851898105111
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2202-53566-0001':
        gname = 'spec-2202-53566-0001'
        T_OIII = 12035.345018778575
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2205-53793-0297':
        gname = 'spec-2205-53793-0297'
        T_OIII = 12565.386684581365
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2097-53491-0355':
        gname = 'spec-2097-53491-0355'
        T_OIII = 11798.295322805352
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2419-54139-0611':
        gname = 'spec-2419-54139-0611'
        T_OIII = 11510.201985913074
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2427-53815-0117':
        gname = 'spec-2427-53815-0117'
        T_OIII = 10754.039074779943
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2433-53820-0275':
        gname = 'spec-2433-53820-0275'
        T_OIII = 10418.253186568521
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2290-53727-0527':
        gname = 'spec-2290-53727-0527'
        T_OIII = 16907.773932530574
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2293-53730-0005':
        gname = 'spec-2293-53730-0005'
        T_OIII = 12063.211844167205
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2359-53826-0222':
        gname = 'spec-2359-53826-0222'
        T_OIII = 10993.01693375783
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2212-53789-0156':
        gname = 'spec-2212-53789-0156'
        T_OIII = 10006.42199587506
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2109-53468-0607':
        gname = 'spec-2109-53468-0607'
        T_OIII = 10427.699985253888
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2356-53786-0466':
        gname = 'spec-2356-53786-0466'
        T_OIII = 12435.098425241405
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2229-53823-0031':
        gname = 'spec-2229-53823-0031'
        T_OIII = 11189.03129173066
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2119-53792-0135':
        gname = 'spec-2119-53792-0135'
        T_OIII = 13020.433563005136
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2122-54178-0393':
        gname = 'spec-2122-54178-0393'
        T_OIII = 12519.914921897609
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2504-54179-0371':
        gname = 'spec-2504-54179-0371'
        T_OIII = 11022.947962474924
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2372-53768-0508':
        gname = 'spec-2372-53768-0508'
        T_OIII = 15148.164001474333
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2483-53852-0254':
        gname = 'spec-2483-53852-0254'
        T_OIII = 14216.955142694998
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2519-54570-0212':
        gname = 'spec-2519-54570-0212'
        T_OIII = 16785.62327235894
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2479-54174-0510':
        gname = 'spec-2479-54174-0510'
        T_OIII = 10184.844570400828
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2482-54175-0389':
        gname = 'spec-2482-54175-0389'
        T_OIII = 9459.634784789014
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2365-53739-0497':
        gname = 'spec-2365-53739-0497'
        T_OIII = 10676.346234921906
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2426-53795-0015':
        gname = 'spec-2426-53795-0015'
        T_OIII = 11518.29854241994
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2167-53889-0598':
        gname = 'spec-2167-53889-0598'
        T_OIII = 13630.23236556919
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2498-54169-0212':
        gname = 'spec-2498-54169-0212'
        T_OIII = 9743.489643856263
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2500-54178-0063':
        gname = 'spec-2500-54178-0063'
        T_OIII = 9720.373926563034
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2217-53794-0352':
        gname = 'spec-2217-53794-0352'
        T_OIII = 13642.591641052915
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2364-53737-0435':
        gname = 'spec-2364-53737-0435'
        T_OIII = 11470.878528718307
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2364-53737-0441':
        gname = 'spec-2364-53737-0441'
        T_OIII = 11175.052295877631
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2364-53737-0546':
        gname = 'spec-2364-53737-0546'
        T_OIII = 10042.764865119534
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2513-54141-0309':
        gname = 'spec-2513-54141-0309'
        T_OIII = 10347.674453240526
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2273-53709-0464':
        gname = 'spec-2273-53709-0464'
        T_OIII = 15926.029193516453
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2407-53771-0068':
        gname = 'spec-2407-53771-0068'
        T_OIII = 10042.764865119534
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2549-54523-0248':
        gname = 'spec-2549-54523-0248'
        T_OIII = 10380.551496008155
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2232-53827-0062':
        gname = 'spec-2232-53827-0062'
        T_OIII = 10928.444888209293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2571-54055-0543':
        gname = 'spec-2571-54055-0543'
        T_OIII = 10452.49777253023
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2522-54570-0387':
        gname = 'spec-2522-54570-0387'
        T_OIII = 15011.489873224315
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2208-53880-0306':
        gname = 'spec-2208-53880-0306'
        T_OIII = 10203.323253131312
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2425-54139-0203':
        gname = 'spec-2425-54139-0203'
        T_OIII = 9938.633223902725
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2429-53799-0383':
        gname = 'spec-2429-53799-0383'
        T_OIII = 10649.186865273103
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2579-54068-0388':
        gname = 'spec-2579-54068-0388'
        T_OIII = 9961.178300666343
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2527-54569-0571':
        gname = 'spec-2527-54569-0571'
        T_OIII = 11530.724860518709
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2529-54585-0291':
        gname = 'spec-2529-54585-0291'
        T_OIII = 9725.843740717663
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2529-54585-0414':
        gname = 'spec-2529-54585-0414'
        T_OIII = 11649.534764548354
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2483-53852-0399':
        gname = 'spec-2483-53852-0399'
        T_OIII = 11128.349786969071
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2483-53852-0572':
        gname = 'spec-2483-53852-0572'
        T_OIII = 11492.22924449282
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2430-53815-0117':
        gname = 'spec-2430-53815-0117'
        T_OIII = 12418.20420184761
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2613-54481-0507':
        gname = 'spec-2613-54481-0507'
        T_OIII = 11108.71660182317
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2493-54115-0164':
        gname = 'spec-2493-54115-0164'
        T_OIII = 12345.260727623503
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2234-53823-0182':
        gname = 'spec-2234-53823-0182'
        T_OIII = 12041.364818561056
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2587-54138-0565':
        gname = 'spec-2587-54138-0565'
        T_OIII = 12251.093532104507
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == 'spec-2291-53714-0567':
    #	 gname = 'spec-2291-53714-0567'
    #	 T_OIII = nan
    #	 full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2480-53851-0136':
        gname = 'spec-2480-53851-0136'
        T_OIII = 14242.749410016355
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2548-54152-0200':
        gname = 'spec-2548-54152-0200'
        T_OIII = 10227.90720842272
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2522-54570-0226':
        gname = 'spec-2522-54570-0226'
        T_OIII = 11414.37483179368
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2344-53740-0444':
        gname = 'spec-2344-53740-0444'
        T_OIII = 12272.745717156748
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2478-54097-0585':
        gname = 'spec-2478-54097-0585'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2489-53857-0209':
        gname = 'spec-2489-53857-0209'
        T_OIII = 10805.332004316992
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2496-54178-0241':
        gname = 'spec-2496-54178-0241'
        T_OIII = 11638.981065280612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2505-53856-0213':
        gname = 'spec-2505-53856-0213'
        T_OIII = 10565.645101262739
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2660-54504-0419':
        gname = 'spec-2660-54504-0419'
        T_OIII = 11103.163057786347
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2577-54086-0315':
        gname = 'spec-2577-54086-0315'
        T_OIII = 13222.603865951582
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2499-54176-0080':
        gname = 'spec-2499-54176-0080'
        T_OIII = 10508.344408574927
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2495-54175-0057':
        gname = 'spec-2495-54175-0057'
        T_OIII = 11073.014217978294
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2271-53726-0443':
        gname = 'spec-2271-53726-0443'
        T_OIII = 10775.991876465538
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2278-53711-0411':
        gname = 'spec-2278-53711-0411'
        T_OIII = 14159.088712577077
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2510-53877-0560':
        gname = 'spec-2510-53877-0560'
        T_OIII = 14038.137285016983
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2516-54241-0245':
        gname = 'spec-2516-54241-0245'
        T_OIII = 12112.512686851956
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2505-53856-0072':
        gname = 'spec-2505-53856-0072'
        T_OIII = 10637.710554366235
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2515-54180-0507':
        gname = 'spec-2515-54180-0507'
        T_OIII = 12571.082253116232
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2527-54569-0630':
        gname = 'spec-2527-54569-0630'
        T_OIII = 23214.271186179794
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2586-54169-0585':
        gname = 'spec-2586-54169-0585'
        T_OIII = 11429.010392348184
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2524-54568-0146':
        gname = 'spec-2524-54568-0146'
        T_OIII = 11342.186641212293
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2286-53700-0472':
        gname = 'spec-2286-53700-0472'
        T_OIII = 10559.70298285448
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2642-54232-0571':
        gname = 'spec-2642-54232-0571'
        T_OIII = 11213.887089206391
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2322-53727-0365':
        gname = 'spec-2322-53727-0365'
        T_OIII = 11062.982807762703
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2606-54154-0474':
        gname = 'spec-2606-54154-0474'
        T_OIII = 11047.952729466202
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2608-54474-0124':
        gname = 'spec-2608-54474-0124'
        T_OIII = 10918.544448136336
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2611-54477-0634':
        gname = 'spec-2611-54477-0634'
        T_OIII = 10456.091807953026
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2520-54584-0205':
        gname = 'spec-2520-54584-0205'
        T_OIII = 11440.803853135543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2582-54139-0510':
        gname = 'spec-2582-54139-0510'
        T_OIII = 11219.496014829632
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2428-53801-0393':
        gname = 'spec-2428-53801-0393'
        T_OIII = 10437.15534986711
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2528-54571-0513':
        gname = 'spec-2528-54571-0513'
        T_OIII = 10038.999149938383
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2361-53762-0304':
        gname = 'spec-2361-53762-0304'
        T_OIII = 11728.993140559112
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2567-54179-0190':
        gname = 'spec-2567-54179-0190'
        T_OIII = 10116.95366566284
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2770-54510-0583':
        gname = 'spec-2770-54510-0583'
        T_OIII = 10418.253186568521
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2588-54174-0369':
        gname = 'spec-2588-54174-0369'
        T_OIII = 14255.664088129159
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2574-54084-0323':
        gname = 'spec-2574-54084-0323'
        T_OIII = 11012.961908937923
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2583-54095-0384':
        gname = 'spec-2583-54095-0384'
        T_OIII = 10446.619288175401
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2585-54097-0334':
        gname = 'spec-2585-54097-0334'
        T_OIII = 10194.079724756852
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2592-54178-0330':
        gname = 'spec-2592-54178-0330'
        T_OIII = 11103.163057786347
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2601-54144-0342':
        gname = 'spec-2601-54144-0342'
        T_OIII = 12003.227686225393
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2776-54554-0372':
        gname = 'spec-2776-54554-0372'
        T_OIII = 10983.057995730782
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2780-54557-0245':
        gname = 'spec-2780-54557-0245'
        T_OIII = 11209.331901904156
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2618-54506-0623':
        gname = 'spec-2618-54506-0623'
        T_OIII = 11073.014217978294
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2481-54086-0081':
        gname = 'spec-2481-54086-0081'
        T_OIII = 11513.619700297988
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2376-53770-0030':
        gname = 'spec-2376-53770-0030'
        T_OIII = 10215.448301863054
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2745-54231-0004':
        gname = 'spec-2745-54231-0004'
        T_OIII = 12479.677077665141
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2587-54138-0509':
        gname = 'spec-2587-54138-0509'
        T_OIII = 11455.83132165492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2398-53768-0241':
        gname = 'spec-2398-53768-0241'
        T_OIII = 13061.80260230364
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2592-54178-0449':
        gname = 'spec-2592-54178-0449'
        T_OIII = 12025.005518437003
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2574-54084-0620':
        gname = 'spec-2574-54084-0620'
        T_OIII = 11686.000278916908
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2576-54086-0384':
        gname = 'spec-2576-54086-0384'
        T_OIII = 11250.04368517661
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2421-54153-0494':
        gname = 'spec-2421-54153-0494'
        T_OIII = 12112.512686851956
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2779-54540-0261':
        gname = 'spec-2779-54540-0261'
        T_OIII = 10565.645101262739
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2742-54233-0048':
        gname = 'spec-2742-54233-0048'
        T_OIII = 12322.90289865434
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2604-54484-0439':
        gname = 'spec-2604-54484-0439'
        T_OIII = 12373.265066298753
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2606-54154-0383':
        gname = 'spec-2606-54154-0383'
        T_OIII = 12096.056725813787
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2612-54480-0516':
        gname = 'spec-2612-54480-0516'
        T_OIII = 12452.015632214747
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2762-54533-0008':
        gname = 'spec-2762-54533-0008'
        T_OIII = 11868.006985422713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2661-54505-0036':
        gname = 'spec-2661-54505-0036'
        T_OIII = 11760.928034558468
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2646-54479-0622':
        gname = 'spec-2646-54479-0622'
        T_OIII = 11083.054724220938
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2864-54467-0595':
        gname = 'spec-2864-54467-0595'
        T_OIII = 11788.52790017836
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2769-54527-0135':
        gname = 'spec-2769-54527-0135'
        T_OIII = 11814.346197280407
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2785-54537-0284':
        gname = 'spec-2785-54537-0284'
        T_OIII = 13936.282828356027
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2645-54477-0042':
        gname = 'spec-2645-54477-0042'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2645-54477-0070':
        gname = 'spec-2645-54477-0070'
        T_OIII = 11209.331901904156
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2654-54231-0527':
        gname = 'spec-2654-54231-0527'
        T_OIII = 21561.303729122017
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2580-54092-0541':
        gname = 'spec-2580-54092-0541'
        T_OIII = 19648.480029200535
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2643-54208-0236':
        gname = 'spec-2643-54208-0236'
        T_OIII = 10795.543096281908
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2489-53857-0154':
        gname = 'spec-2489-53857-0154'
        T_OIII = 10934.594503747658
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2480-53851-0364':
        gname = 'spec-2480-53851-0364'
        T_OIII = 11575.85934665115
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2488-54149-0396':
        gname = 'spec-2488-54149-0396'
        T_OIII = 9866.83130519079
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2488-54149-0535':
        gname = 'spec-2488-54149-0535'
        T_OIII = 21737.900577737124
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2653-54230-0380':
        gname = 'spec-2653-54230-0380'
        T_OIII = 12124.82198008786
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2748-54234-0552':
        gname = 'spec-2748-54234-0552'
        T_OIII = 11113.230901627261
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2729-54419-0548':
        gname = 'spec-2729-54419-0548'
        T_OIII = 10092.95193226424
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2499-54176-0604':
        gname = 'spec-2499-54176-0604'
        T_OIII = 10578.861782951948
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2501-54084-0470':
        gname = 'spec-2501-54084-0470'
        T_OIII = 10338.300152334337
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2661-54505-0382':
        gname = 'spec-2661-54505-0382'
        T_OIII = 12014.111667781179
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2663-54234-0324':
        gname = 'spec-2663-54234-0324'
        T_OIII = 11721.022975334794
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2881-54502-0550':
        gname = 'spec-2881-54502-0550'
        T_OIII = 13192.677223723627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2742-54233-0385':
        gname = 'spec-2742-54233-0385'
        T_OIII = 11549.659771691455
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2742-54233-0565':
        gname = 'spec-2742-54233-0565'
        T_OIII = 22963.157925311705
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2508-53875-0615':
        gname = 'spec-2508-53875-0615'
        T_OIII = 14483.5760702146
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2510-53877-0330':
        gname = 'spec-2510-53877-0330'
        T_OIII = 10879.032297991607
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2744-54272-0182':
        gname = 'spec-2744-54272-0182'
        T_OIII = 10770.772981133998
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2926-54625-0136':
        gname = 'spec-2926-54625-0136'
        T_OIII = 13716.982993362524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2748-54234-0343':
        gname = 'spec-2748-54234-0343'
        T_OIII = 11168.767446886002
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2760-54506-0362':
        gname = 'spec-2760-54506-0362'
        T_OIII = 10333.13176005107
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2756-54508-0185':
        gname = 'spec-2756-54508-0185'
        T_OIII = 12139.989053287203
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2782-54592-0293':
        gname = 'spec-2782-54592-0293'
        T_OIII = 15251.485596321565
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2772-54529-0566':
        gname = 'spec-2772-54529-0566'
        T_OIII = 10748.326923558492
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2970-54589-0280':
        gname = 'spec-2970-54589-0280'
        T_OIII = 11093.104334738524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2777-54554-0633':
        gname = 'spec-2777-54554-0633'
        T_OIII = 11118.268246429308
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2763-54507-0169':
        gname = 'spec-2763-54507-0169'
        T_OIII = 10766.229548623747
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2788-54553-0205':
        gname = 'spec-2788-54553-0205'
        T_OIII = 11660.098033431159
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2768-54265-0449':
        gname = 'spec-2768-54265-0449'
        T_OIII = 11846.513484698893
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2770-54510-0330':
        gname = 'spec-2770-54510-0330'
        T_OIII = 11378.223487652196
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2773-54533-0362':
        gname = 'spec-2773-54533-0362'
        T_OIII = 11122.43881562274
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2794-54537-0201':
        gname = 'spec-2794-54537-0201'
        T_OIII = 10686.027060658256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2759-54534-0475':
        gname = 'spec-2759-54534-0475'
        T_OIII = 10869.176622420213
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2647-54495-0242':
        gname = 'spec-2647-54495-0242'
        T_OIII = 10998.343518375543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2603-54479-0322':
        gname = 'spec-2603-54479-0322'
        T_OIII = 12395.714268927213
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2767-54243-0015':
        gname = 'spec-2767-54243-0015'
        T_OIII = 13138.979851327573
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2865-54503-0304':
        gname = 'spec-2865-54503-0304'
        T_OIII = 10918.544448136336
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2954-54561-0037':
        gname = 'spec-2954-54561-0037'
        T_OIII = 10983.057995730782
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2918-54554-0391':
        gname = 'spec-2918-54554-0391'
        T_OIII = 11539.196552486676
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2927-54621-0442':
        gname = 'spec-2927-54621-0442'
        T_OIII = 12014.111667781179
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2877-54523-0344':
        gname = 'spec-2877-54523-0344'
        T_OIII = 12063.211844167205
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2747-54233-0293':
        gname = 'spec-2747-54233-0293'
        T_OIII = 15572.753550174768
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2750-54242-0045':
        gname = 'spec-2750-54242-0045'
        T_OIII = 12474.607714463225
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2882-54498-0236':
        gname = 'spec-2882-54498-0236'
        T_OIII = 10879.032297991607
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2788-54553-0451':
        gname = 'spec-2788-54553-0451'
        T_OIII = 11776.007947121634
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2926-54625-0538':
        gname = 'spec-2926-54625-0538'
        T_OIII = 10968.136502632646
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2947-54533-0059':
        gname = 'spec-2947-54533-0059'
        T_OIII = 10715.122238510845
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2974-54592-0526':
        gname = 'spec-2974-54592-0526'
        T_OIII = 10584.81468233188
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2880-54509-0095':
        gname = 'spec-2880-54509-0095'
        T_OIII = 14181.564284025484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2911-54631-0344':
        gname = 'spec-2911-54631-0344'
        T_OIII = 10305.073815434032
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2960-54561-0379':
        gname = 'spec-2960-54561-0379'
        T_OIII = 10993.01693375783
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3239-54888-0207':
        gname = 'spec-3239-54888-0207'
        T_OIII = 18039.739743126604
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2757-54509-0068':
        gname = 'spec-2757-54509-0068'
        T_OIII = 11388.540743794836
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2764-54535-0172':
        gname = 'spec-2764-54535-0172'
        T_OIII = 11638.981065280612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2969-54586-0054':
        gname = 'spec-2969-54586-0054'
        T_OIII = 10671.509112671123
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3242-54889-0544':
        gname = 'spec-3242-54889-0544'
        T_OIII = 10948.272708303248
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2974-54592-0360':
        gname = 'spec-2974-54592-0360'
        T_OIII = 10498.8245514363
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2881-54502-0455':
        gname = 'spec-2881-54502-0455'
        T_OIII = 11679.245736833644
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2956-54525-0457':
        gname = 'spec-2956-54525-0457'
        T_OIII = 18039.739743126604
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-2956-54525-0606':
        gname = 'spec-2956-54525-0606'
        T_OIII = 10385.74360654094
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3783-55246-0123':
        gname = 'spec-3783-55246-0123'
        T_OIII = 11586.355810012077
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3317-54908-0565':
        gname = 'spec-3317-54908-0565'
        T_OIII = 20089.662999872773
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3587-55182-0148':
        gname = 'spec-3587-55182-0148'
        T_OIII = 19121.429541360983
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3767-55214-0459':
        gname = 'spec-3767-55214-0459'
        T_OIII = 14181.564284025484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3691-55274-0169':
        gname = 'spec-3691-55274-0169'
        T_OIII = 14463.898801784284
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3925-55338-0501':
        gname = 'spec-3925-55338-0501'
        T_OIII = 24678.823566307234
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3615-56544-0376':
        gname = 'spec-3615-56544-0376'
        T_OIII = 15467.255517028703
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3848-55647-0395':
        gname = 'spec-3848-55647-0395'
        T_OIII = 13791.779992446884
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3952-55330-0899':
        gname = 'spec-3952-55330-0899'
        T_OIII = 10316.029908618977
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3961-55654-0019':
        gname = 'spec-3961-55654-0019'
        T_OIII = 11905.71448917048
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3825-55533-0627':
        gname = 'spec-3825-55533-0627'
        T_OIII = 13835.599784295056
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3801-55509-0539':
        gname = 'spec-3801-55509-0539'
        T_OIII = 14761.885289850534
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3828-55539-0803':
        gname = 'spec-3828-55539-0803'
        T_OIII = 10556.073333529317
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3973-55323-0647':
        gname = 'spec-3973-55323-0647'
        T_OIII = 12536.947516637038
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3830-55574-0381':
        gname = 'spec-3830-55574-0381'
        T_OIII = 15134.440778138678
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3840-55574-0619':
        gname = 'spec-3840-55574-0619'
        T_OIII = 13162.818314290524
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3663-55176-0387':
        gname = 'spec-3663-55176-0387'
        T_OIII = 12952.754388580213
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3932-55337-0561':
        gname = 'spec-3932-55337-0561'
        T_OIII = 12477.142396064184
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3857-55272-0747':
        gname = 'spec-3857-55272-0747'
        T_OIII = 20494.262438468515
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3864-55649-0861':
        gname = 'spec-3864-55649-0861'
        T_OIII = 12702.795276164543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3942-55338-0939':
        gname = 'spec-3942-55338-0939'
        T_OIII = 11067.478522516425
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3857-55272-0595':
        gname = 'spec-3857-55272-0595'
        T_OIII = 13397.521980857759
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4181-55685-0071':
        gname = 'spec-4181-55685-0071'
        T_OIII = 15107.031617171797
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3968-55590-0761':
        gname = 'spec-3968-55590-0761'
        T_OIII = 11713.058226051426
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4020-55332-0027':
        gname = 'spec-4020-55332-0027'
        T_OIII = 9645.776642795578
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3926-55327-0774':
        gname = 'spec-3926-55327-0774'
        T_OIII = 9637.038214356204
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4006-55328-0874':
        gname = 'spec-4006-55328-0874'
        T_OIII = 10775.991876465538
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4296-55499-0913':
        gname = 'spec-4296-55499-0913'
        T_OIII = 16122.091383183888
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4321-55504-0583':
        gname = 'spec-4321-55504-0583'
        T_OIII = 13138.979851327573
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4093-55475-0275':
        gname = 'spec-4093-55475-0275'
        T_OIII = 9884.733006486616
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4189-55679-0519':
        gname = 'spec-4189-55679-0519'
        T_OIII = 14444.248266730488
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4317-55480-0057':
        gname = 'spec-4317-55480-0057'
        T_OIII = 11981.489294658955
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4371-55830-0779':
        gname = 'spec-4371-55830-0779'
        T_OIII = 10647.356347022687
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4454-55536-0235':
        gname = 'spec-4454-55536-0235'
        T_OIII = 12194.55721966369
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4297-55806-0323':
        gname = 'spec-4297-55806-0323'
        T_OIII = 14869.309006222302
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4490-55629-0573':
        gname = 'spec-4490-55629-0573'
        T_OIII = 12261.052620703023
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4478-55600-0939':
        gname = 'spec-4478-55600-0939'
        T_OIII = 15750.185485418875
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4276-55505-0279':
        gname = 'spec-4276-55505-0279'
        T_OIII = 14483.5760702146
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4351-55538-0807':
        gname = 'spec-4351-55538-0807'
        T_OIII = 11782.26625491757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4379-55881-0577':
        gname = 'spec-4379-55881-0577'
        T_OIII = 9319.221589255085
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4424-55532-0766':
        gname = 'spec-4424-55532-0766'
        T_OIII = 13464.47389542078
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4485-55836-0407':
        gname = 'spec-4485-55836-0407'
        T_OIII = 17925.650606419076
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4232-55447-0200':
        gname = 'spec-4232-55447-0200'
        T_OIII = 14181.564284025484
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4570-55623-0245':
        gname = 'spec-4570-55623-0245'
        T_OIII = 16151.342198899576
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4532-55559-0651':
        gname = 'spec-4532-55559-0651'
        T_OIII = 12311.739172898178
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4548-55565-0823':
        gname = 'spec-4548-55565-0823'
        T_OIII = 17643.5745127752
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4320-55894-0831':
        gname = 'spec-4320-55894-0831'
        T_OIII = 13885.850052981183
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4445-55869-0993':
        gname = 'spec-4445-55869-0993'
        T_OIII = 10749.334725405466
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4446-55589-0041':
        gname = 'spec-4446-55589-0041'
        T_OIII = 13234.593522430534
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4698-55623-0355':
        gname = 'spec-4698-55623-0355'
        T_OIII = 11617.902340920347
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4531-55563-0115':
        gname = 'spec-4531-55563-0115'
        T_OIII = 23973.341180180905
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4557-55588-0819':
        gname = 'spec-4557-55588-0819'
        T_OIII = 9631.016324844102
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4724-55742-0563':
        gname = 'spec-4724-55742-0563'
        T_OIII = 13312.791727617876
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4602-55644-0075':
        gname = 'spec-4602-55644-0075'
        T_OIII = 11209.331901904156
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4537-55806-0899':
        gname = 'spec-4537-55806-0899'
        T_OIII = 11905.71448917048
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4483-55587-0301':
        gname = 'spec-4483-55587-0301'
        T_OIII = 10565.645101262739
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4706-55705-0119':
        gname = 'spec-4706-55705-0119'
        T_OIII = 12356.454849202595
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4538-55860-0881':
        gname = 'spec-4538-55860-0881'
        T_OIII = 9898.025920130092
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4658-55592-0859':
        gname = 'spec-4658-55592-0859'
        T_OIII = 14307.44001124752
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4655-55620-0387':
        gname = 'spec-4655-55620-0387'
        T_OIII = 9734.054233655163
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4627-55626-0187':
        gname = 'spec-4627-55626-0187'
        T_OIII = 17310.93614322934
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4724-55742-0101':
        gname = 'spec-4724-55742-0101'
        T_OIII = 11002.984902097443
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4775-55708-0842':
        gname = 'spec-4775-55708-0842'
        T_OIII = 18562.193691014814
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4733-55649-0430':
        gname = 'spec-4733-55649-0430'
        T_OIII = 10705.415054489516
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4770-55928-0857':
        gname = 'spec-4770-55928-0857'
        T_OIII = 11026.738144142759
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4862-55685-0045':
        gname = 'spec-4862-55685-0045'
        T_OIII = 9766.202627381766
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4895-55708-0827':
        gname = 'spec-4895-55708-0827'
        T_OIII = 15320.757879677238
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4977-56047-0378':
        gname = 'spec-4977-56047-0378'
        T_OIII = 16165.98750224542
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4764-55646-0337':
        gname = 'spec-4764-55646-0337'
        T_OIII = 13587.06300866288
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4760-55656-0232':
        gname = 'spec-4760-55656-0232'
        T_OIII = 23140.74748304121
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4763-55869-0183':
        gname = 'spec-4763-55869-0183'
        T_OIII = 12222.792687463036
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4862-55685-0295':
        gname = 'spec-4862-55685-0295'
        T_OIII = 14516.2041829868
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4770-55928-0407':
        gname = 'spec-4770-55928-0407'
        T_OIII = 10815.129788489168
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4691-55651-0745':
        gname = 'spec-4691-55651-0745'
        T_OIII = 11445.453104749757
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4801-55653-0952':
        gname = 'spec-4801-55653-0952'
        T_OIII = 9085.69756257142
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4903-55927-0800':
        gname = 'spec-4903-55927-0800'
        T_OIII = 12812.842159201124
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4966-55712-0793':
        gname = 'spec-4966-55712-0793'
        T_OIII = 11628.436926958453
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5165-56063-0345':
        gname = 'spec-5165-56063-0345'
        T_OIII = 10903.710602618274
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4783-55652-0053':
        gname = 'spec-4783-55652-0053'
        T_OIII = 13312.791727617876
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5177-56245-0697':
        gname = 'spec-5177-56245-0697'
        T_OIII = 13127.076814439111
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5200-56091-0681':
        gname = 'spec-5200-56091-0681'
        T_OIII = 9628.3077023485
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5045-56181-0997':
        gname = 'spec-5045-56181-0997'
        T_OIII = 13873.27040164763
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5184-56352-0557':
        gname = 'spec-5184-56352-0557'
        T_OIII = 11503.189130987197
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5184-56352-0343':
        gname = 'spec-5184-56352-0343'
        T_OIII = 11702.4469787766
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5128-55912-0299':
        gname = 'spec-5128-55912-0299'
        T_OIII = 9645.776642795578
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5301-55987-0625':
        gname = 'spec-5301-55987-0625'
        T_OIII = 11260.244713697612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5391-56000-0177':
        gname = 'spec-5391-56000-0177'
        T_OIII = 11178.894777806525
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5316-55955-0133':
        gname = 'spec-5316-55955-0133'
        T_OIII = 11194.102994920871
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5398-56011-0062':
        gname = 'spec-5398-56011-0062'
        T_OIII = 14171.927530970173
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5363-55956-0305':
        gname = 'spec-5363-55956-0305'
        T_OIII = 15251.485596321565
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5371-55976-0568':
        gname = 'spec-5371-55976-0568'
        T_OIII = 13741.87014554258
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5398-56011-0229':
        gname = 'spec-5398-56011-0229'
        T_OIII = 12497.240786234663
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5466-56033-0441':
        gname = 'spec-5466-56033-0441'
        T_OIII = 12639.6307982468
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5649-55912-0445':
        gname = 'spec-5649-55912-0445'
        T_OIII = 16195.31796006577
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5488-56013-0339':
        gname = 'spec-5488-56013-0339'
        T_OIII = 16305.78191335606
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5177-56245-0742':
        gname = 'spec-5177-56245-0742'
        T_OIII = 12859.17687175749
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5362-56017-0341':
        gname = 'spec-5362-56017-0341'
        T_OIII = 10551.290701921758
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5444-56038-0595':
        gname = 'spec-5444-56038-0595'
        T_OIII = 12853.350774892599
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5450-55986-0047':
        gname = 'spec-5450-55986-0047'
        T_OIII = 9929.629486989834
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5410-56016-0303':
        gname = 'spec-5410-56016-0303'
        T_OIII = 11739.628458687444
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5466-56033-0198':
        gname = 'spec-5466-56033-0198'
        T_OIII = 10277.573857418769
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5427-56001-0289':
        gname = 'spec-5427-56001-0289'
        T_OIII = 11087.558595690873
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5436-56015-0463':
        gname = 'spec-5436-56015-0463'
        T_OIII = 11644.256719254627
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5321-55945-0361':
        gname = 'spec-5321-55945-0361'
        T_OIII = 11199.176996985096
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5495-55896-0529':
        gname = 'spec-5495-55896-0529'
        T_OIII = 10135.309171796935
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5777-56280-0893':
        gname = 'spec-5777-56280-0893'
        T_OIII = 11260.244713697612
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5866-56035-0797':
        gname = 'spec-5866-56035-0797'
        T_OIII = 9446.635386570233
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5808-56325-0383':
        gname = 'spec-5808-56325-0383'
        T_OIII = 13710.76825264594
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5344-55924-0017':
        gname = 'spec-5344-55924-0017'
        T_OIII = 17146.976141241066
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5741-55980-0039':
        gname = 'spec-5741-55980-0039'
        T_OIII = 12139.989053287203
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5767-56245-0126':
        gname = 'spec-5767-56245-0126'
        T_OIII = 10775.991876465538
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5865-56067-0822':
        gname = 'spec-5865-56067-0822'
        T_OIII = 11088.598210536786
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5785-56269-0387':
        gname = 'spec-5785-56269-0387'
        T_OIII = 12474.607714463225
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5963-56191-0105':
        gname = 'spec-5963-56191-0105'
        T_OIII = 9938.633223902725
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5852-56034-0153':
        gname = 'spec-5852-56034-0153'
        T_OIII = 10809.891938367335
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5438-56002-0426':
        gname = 'spec-5438-56002-0426'
        T_OIII = 10666.674179366131
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5888-56041-0705':
        gname = 'spec-5888-56041-0705'
        T_OIII = 11042.947242458555
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5796-56274-0105':
        gname = 'spec-5796-56274-0105'
        T_OIII = 9554.416964471655
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6140-56189-0595':
        gname = 'spec-6140-56189-0595'
        T_OIII = 11771.592309797865
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5743-56011-0352':
        gname = 'spec-5743-56011-0352'
        T_OIII = 14755.19713442155
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6005-56090-0873':
        gname = 'spec-6005-56090-0873'
        T_OIII = 8536.83517072629
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-5880-56042-0587':
        gname = 'spec-5880-56042-0587'
        T_OIII = 15753.754255768452
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6050-56089-0885':
        gname = 'spec-6050-56089-0885'
        T_OIII = 12373.265066298753
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6145-56266-0163':
        gname = 'spec-6145-56266-0163'
        T_OIII = 10705.415054489516
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6385-56356-0119':
        gname = 'spec-6385-56356-0119'
        T_OIII = 16165.98750224542
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6281-56295-0829':
        gname = 'spec-6281-56295-0829'
        T_OIII = 10942.115400143155
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6293-56561-0573':
        gname = 'spec-6293-56561-0573'
        T_OIII = 12200.656653687713
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6167-56189-0231':
        gname = 'spec-6167-56189-0231'
        T_OIII = 23004.82070247622
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6278-56266-0215':
        gname = 'spec-6278-56266-0215'
        T_OIII = 15601.007685041946
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6114-56209-0437':
        gname = 'spec-6114-56209-0437'
        T_OIII = 12395.714268927213
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6125-56273-0741':
        gname = 'spec-6125-56273-0741'
        T_OIII = 10394.99848836848
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6402-56334-0033':
        gname = 'spec-6402-56334-0033'
        T_OIII = 14470.454919090893
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6431-56311-0161':
        gname = 'spec-6431-56311-0161'
        T_OIII = 12211.71965486194
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6308-56215-0641':
        gname = 'spec-6308-56215-0641'
        T_OIII = 10427.699985253888
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6488-56364-0823':
        gname = 'spec-6488-56364-0823'
        T_OIII = 13156.23786996223
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6176-56264-0143':
        gname = 'spec-6176-56264-0143'
        T_OIII = 11316.515936324256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6380-56340-0182':
        gname = 'spec-6380-56340-0182'
        T_OIII = 11209.331901904156
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6409-56306-0367':
        gname = 'spec-6409-56306-0367'
        T_OIII = 15425.256691430432
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6439-56358-0552':
        gname = 'spec-6439-56358-0552'
        T_OIII = 15671.867498544503
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6461-56329-0529':
        gname = 'spec-6461-56329-0529'
        T_OIII = 11623.71334558001
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6475-56337-0385':
        gname = 'spec-6475-56337-0385'
        T_OIII = 10874.613131881746
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6519-56566-0861':
        gname = 'spec-6519-56566-0861'
        T_OIII = 14391.977260374884
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6622-56365-0320':
        gname = 'spec-6622-56365-0320'
        T_OIII = 9610.870398947525
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6458-56274-0655':
        gname = 'spec-6458-56274-0655'
        T_OIII = 16027.39159099714
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6495-56339-0478':
        gname = 'spec-6495-56339-0478'
        T_OIII = 13599.383140134987
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6664-56383-0845':
        gname = 'spec-6664-56383-0845'
        T_OIII = 10761.351701709767
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6688-56412-0259':
        gname = 'spec-6688-56412-0259'
        T_OIII = 11008.488339349717
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6730-56425-0898':
        gname = 'spec-6730-56425-0898'
        T_OIII = 16063.748579488312
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6620-56368-0337':
        gname = 'spec-6620-56368-0337'
        T_OIII = 11451.177854310561
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6813-56419-0498':
        gname = 'spec-6813-56419-0498'
        T_OIII = 10070.1086130086
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6456-56339-0275':
        gname = 'spec-6456-56339-0275'
        T_OIII = 12725.842355781295
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6671-56388-0473':
        gname = 'spec-6671-56388-0473'
        T_OIII = 15732.351333605193
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6706-56385-0307':
        gname = 'spec-6706-56385-0307'
        T_OIII = 16254.138617763192
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6727-56369-0267':
        gname = 'spec-6727-56369-0267'
        T_OIII = 11047.952729466202
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6976-56448-0597':
        gname = 'spec-6976-56448-0597'
        T_OIII = 11607.377298512543
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6821-56396-0363':
        gname = 'spec-6821-56396-0363'
        T_OIII = 12714.31359382933
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6716-56401-0835':
        gname = 'spec-6716-56401-0835'
        T_OIII = 13574.75403840583
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7049-56570-0467':
        gname = 'spec-7049-56570-0467'
        T_OIII = 10533.33397383276
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6649-56364-0285':
        gname = 'spec-6649-56364-0285'
        T_OIII = 18663.411020305255
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6728-56426-0277':
        gname = 'spec-6728-56426-0277'
        T_OIII = 14775.27069679468
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6739-56393-0901':
        gname = 'spec-6739-56393-0901'
        T_OIII = 10367.582572914509
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7132-56565-0369':
        gname = 'spec-7132-56565-0369'
        T_OIII = 19069.508398426922
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7128-56567-0023':
        gname = 'spec-7128-56567-0023'
        T_OIII = 16041.924501562851
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7261-56603-0267':
        gname = 'spec-7261-56603-0267'
        T_OIII = 12882.507679404664
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7132-56565-0949':
        gname = 'spec-7132-56565-0949'
        T_OIII = 23098.838536291518
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7236-56605-0764':
        gname = 'spec-7236-56605-0764'
        T_OIII = 13246.594050585314
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6790-56430-0085':
        gname = 'spec-6790-56430-0085'
        T_OIII = 13397.521980857759
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7295-57067-0731':
        gname = 'spec-7295-57067-0731'
        T_OIII = 11424.724868332823
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7315-56685-0165':
        gname = 'spec-7315-56685-0165'
        T_OIII = 24005.955475678864
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7153-56904-0984':
        gname = 'spec-7153-56904-0984'
        T_OIII = 12674.045113686818
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7564-56804-0220':
        gname = 'spec-7564-56804-0220'
        T_OIII = 18595.8716608284
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7604-56947-0944':
        gname = 'spec-7604-56947-0944'
        T_OIII = 17444.813265163477
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-8057-57190-0321':
        gname = 'spec-8057-57190-0321'
        T_OIII = 10686.027060658256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-8068-57185-0168':
        gname = 'spec-8068-57185-0168'
        T_OIII = 18731.195679661694
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-8199-57428-0845':
        gname = 'spec-8199-57428-0845'
        T_OIII = 12139.989053287203
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-8208-57430-0244':
        gname = 'spec-8208-57430-0244'
        T_OIII = 14688.482014654513
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7870-57016-0041':
        gname = 'spec-7870-57016-0041'
        T_OIII = 10295.898997914053
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-7705-57332-0036':
        gname = 'spec-7705-57332-0036'
        T_OIII = 10938.354305551531
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-3851-55302-0325':
        gname = 'spec-3851-55302-0325'
        T_OIII = 11378.223487652196
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4075-55352-0777':
        gname = 'spec-4075-55352-0777'
        T_OIII = 16671.908621018792
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-4871-55928-0299':
        gname = 'spec-4871-55928-0299'
        T_OIII = 13537.893993857477
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    elif galaxyname == 'spec-6290-56238-0843':
        gname = 'spec-6290-56238-0843'
        T_OIII = 11923.215666232256
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    else:
        print('Galaxy not known: {0:s}'.format(galaxyname))
        return None
    outdict['full_tbl'] = full_tbl
    outdict['T_OIII'] = T_OIII
    return outdict