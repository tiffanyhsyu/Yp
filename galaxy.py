import os
from astropy.table import Table


def load_AOS2015(galaxyname):
    outdict = dict()
    dir = '/test_data/optical+nir/'
    # if galaxyname == "IZw18SE1":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "SBS0335-052E1":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "SBS0335-052E3":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "J0519+0007":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    if galaxyname == "SBS0940+5442":
        gname = 'SBS0940+5442_nir'
        T_OIII = 18600.0
        full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "Tol65":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "SBS1415+437No13":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "SBS1415+437No2":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "CGCG007-025No2":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "Mrk209":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "SBS1030+583":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "Mrk71No1":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "SBS1152+579":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "Mrk59":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "SBS1135+581":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    # elif galaxyname == "Mrk450No1":
    #     gname =
    #     T_OIII =
    #     full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
    outdict["full_tbl"] = full_tbl
    outdict["T_OIII"] = T_OIII
    return outdict