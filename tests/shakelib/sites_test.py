#!/usr/bin/env python

# stdlib imports
import sys
import os.path

# third party imports
import numpy as np
import pytest

# local imports
from shakelib.sites import Sites
import shakelib.sites as sites


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)


def test_depthpars():
    vs30 = np.linspace(200, 700, 6)
    cy14 = sites.Sites._z1pt0_from_vs30_cy14_cal(vs30)
    cb14 = sites.Sites._z2pt5_from_vs30_cb14_cal(vs30)
    ask14 = sites.Sites._z1pt0_from_vs30_ask14_cal(vs30)
    cy08 = sites.Sites._z1pt0_from_vs30_cy08(vs30)
    cb07 = sites.Sites._z2pt5_from_z1pt0_cb07(cy08)

    cy14t = np.array(
        [509.34591289, 458.77871089, 355.71703571, 228.88509539,
         125.83403858, 63.32230756])
    cb14t = np.array(
        [2.79470046489, 1.75746573662, 1.26461099128, 0.97969726864,
         0.79525892497, 0.66668613683])
    ask14t = np.array(
        [494.57881929, 453.37522656, 365.19411327, 247.50118189,
         142.44704445, 73.48751385])
    cy08t = np.array(
        [336.55300716, 315.06784809, 215.89565373, 111.17436474,
         57.50382896, 32.18118661])
    cb07t = np.array(
        [1728.90806073, 1651.6689139, 1295.14487514, 918.67184124,
         725.72626513, 634.69136587])

    np.testing.assert_allclose(cy14, cy14t)
    np.testing.assert_allclose(cb14, cb14t)
    np.testing.assert_allclose(ask14, ask14t)
    np.testing.assert_allclose(cy08, cy08t)
    np.testing.assert_allclose(cb07, cb07t)


def test_sites(vs30file=None):
    vs30file = os.path.join(homedir, 'sites_data/Vs30_test.grd')
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 0.0083 * 5
    yspan = 0.0083 * 5
    mysite = Sites.fromCenter(cx, cy, xspan, yspan, dx, dy,
                              vs30File=vs30file, padding=True,
                              resample=False)
    grd = mysite.getVs30Grid().getData()
    grd_target = np.array([[426.64892578, 398.89712524, 428.88549805,
                            428.78335571, 428.58578491, 430.54354858,
                            433.59750366],
                           [426.20635986, 425.57946777, 428.21954346,
                            426.06726074, 421.86233521, 423.53192139,
                            426.25296021],
                           [428.14602661, 430.05944824, 429.3427124,
                            426.13626099, 409.76391602, 383.07299805,
                            372.39117432],
                           [432.64077759, 434.55209351, 432.21600342,
                            395.53771973, 419.31866455, 421.67749023,
                            426.23449707],
                           [345.14605713, 403.78097534, 385.49118042,
                            413.04779053, 428.22869873, 427.00268555,
                            426.8951416],
                           [336.48217773, 347.82220459, 425.96798706,
                            432.0640564, 429.40097046, 427.74179077,
                            427.00006104],
                           [330.57504272, 392.33255005, 430.33862305,
                            432.01391602, 429.43969727, 427.30435181,
                            425.96151733]])
    np.testing.assert_allclose(grd, grd_target)

    sc = mysite.getSitesContext()
    scr = mysite.getSitesContext(rock_vs30=760.0)

    grd = sc.backarc
    grdt = np.array([[False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False]],
                    dtype=bool)
    np.testing.assert_allclose(grd, grdt)

    grd = sc.lats
    grdt = np.array([34.125,
                     34.11666667,
                     34.10833333,
                     34.1,
                     34.09166667,
                     34.08333333,
                     34.075])
    np.testing.assert_allclose(grd, grdt)

    grd = sc.lons
    grdt = np.array([-118.225, -118.21666667, -118.20833333, -118.2,
                     -118.19166667, -118.18333333, -118.175])
    np.testing.assert_allclose(grd, grdt)

    grd = sc.vs30
    grdt = np.array(
        [[426.64892578, 398.89712524, 428.88549805, 428.78335571,
          428.58578491, 430.54354858, 433.59750366],
         [426.20635986, 425.57946777, 428.21954346, 426.06726074,
          421.86233521, 423.53192139, 426.25296021],
         [428.14602661, 430.05944824, 429.3427124, 426.13626099,
          409.76391602, 383.07299805, 372.39117432],
         [432.64077759, 434.55209351, 432.21600342, 395.53771973,
          419.31866455, 421.67749023, 426.23449707],
         [345.14605713, 403.78097534, 385.49118042, 413.04779053,
          428.22869873, 427.00268555, 426.8951416],
         [336.48217773, 347.82220459, 425.96798706, 432.0640564,
          429.40097046, 427.74179077, 427.00006104],
         [330.57504272, 392.33255005, 430.33862305, 432.01391602,
          429.43969727, 427.30435181, 425.96151733]])
    np.testing.assert_allclose(grd, grdt)

    grd = scr.vs30
    grdt = np.array([[760., 760., 760., 760., 760., 760., 760.],
                     [760., 760., 760., 760., 760., 760., 760.],
                     [760., 760., 760., 760., 760., 760., 760.],
                     [760., 760., 760., 760., 760., 760., 760.],
                     [760., 760., 760., 760., 760., 760., 760.],
                     [760., 760., 760., 760., 760., 760., 760.],
                     [760., 760., 760., 760., 760., 760., 760.]])
    np.testing.assert_allclose(grd, grdt)

    grd = sc.vs30measured
    grdt = np.array([[False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False],
                     [False, False, False, False, False, False, False]],
                    dtype=bool)
    np.testing.assert_allclose(grd, grdt)

    grd = sc.z1pt0_ask14_cal
    grdt = np.array(
        [[335.06579012, 366.39725582, 332.4593083, 332.57855835,
          332.80916279, 330.52077553, 326.93706615],
         [335.58036296, 336.30856215, 333.23643861, 335.74201097,
          340.60915547, 338.68122322, 335.52619949],
         [333.32217548, 331.08730312, 331.92526979, 335.66183031,
          354.37798647, 383.19730994, 394.00166561],
         [328.06152136, 325.81357093, 328.56025522, 370.03738338,
          343.53423548, 340.82221914, 335.54765968],
         [419.3225465, 361.03999163, 380.68893324, 350.67806177,
          333.2257608, 334.65418562, 334.77934134],
         [426.64586941, 416.9870982, 335.85735317, 328.73858131,
          331.85719411, 333.79341356, 334.65724021],
         [431.42758271, 373.4746626, 330.76064671, 328.79741743,
          331.81193751, 334.30299304, 335.86486938]])
    np.testing.assert_allclose(grd, grdt, atol=1e-2)

    grd = sc.z1pt0_cy08
    grdt = np.array(
        [[183.2043947, 217.27787758, 180.56690603, 180.68687904,
          180.91907101, 178.625999, 175.08452094],
         [183.72884833, 184.47314285, 181.34994816, 183.89385557,
          188.91901585, 186.91536614, 183.67358657],
         [181.43651087, 179.19139187, 180.03044953, 183.81199342,
          203.71493023, 237.03939223, 250.05192732],
         [176.18920323, 173.98674225, 176.68107716, 221.49002942,
          191.99154431, 189.14149784, 183.69548028],
         [280.25774719, 211.16371072, 234.04298229, 199.65723106,
          181.33916991, 182.78577761, 182.91298179],
         [288.5669674, 277.54780981, 184.01166928, 176.85723522,
          179.96216197, 181.91290358, 182.78888133],
         [293.80330679, 225.506479, 178.86520526, 176.91538893,
          179.91677656, 182.42922842, 184.01934871]])
    np.testing.assert_allclose(grd, grdt, atol=1e-2)

    grd = sc.z1pt0_cy14_cal
    grdt = np.array(
        [[322.09231215, 357.07647045, 319.22097485, 319.3522119,
          319.60603217, 317.08933532, 313.15733537],
         [322.65987959, 323.46347258, 320.07644714, 322.83822344,
          328.21884905, 326.08502532, 322.60012698],
         [320.17085962, 317.71195569, 318.63340828, 322.74975843,
          343.55342441, 376.19315199, 388.62006436],
         [314.38985785, 311.92697475, 314.93687894, 361.19727923,
          331.46256533, 328.4548676, 322.62380134],
         [418.15243574, 351.03313981, 373.32295288, 339.41631214,
          320.0646893, 321.63848529, 321.7764637],
         [426.80110483, 415.40447342, 322.96549285, 315.13252356,
          318.55852723, 320.68989724, 321.64185267],
         [432.47426567, 365.09924657, 317.35292206, 315.19707981,
          318.50874868, 321.25138528, 322.9737867]])
    np.testing.assert_allclose(grd, grdt, atol=1e-2)

    grd = sc.z2pt5_cb07
    grdt = np.array(
        [[1177.61979893, 1300.11396989, 1168.13802718, 1168.56933015,
          1169.4040603, 1161.16046639, 1148.4288528],
         [1179.50520975, 1182.18094855, 1170.95306363, 1180.09841078,
          1198.16386197, 1190.96074128, 1179.30654372],
         [1171.26425658, 1163.19305376, 1166.20946607, 1179.80411634,
          1251.35517419, 1371.15661506, 1417.9366787],
         [1152.40018562, 1144.48233838, 1154.16847239, 1315.25665576,
          1209.20960179, 1198.96368474, 1179.3852516],
         [1526.52660116, 1278.13354004, 1360.38452132, 1236.76774567,
          1170.91431583, 1176.11487052, 1176.57216955],
         [1556.3982478, 1516.78437627, 1180.52195107, 1154.8017606,
          1165.96397227, 1172.97688837, 1176.12602836],
         [1575.22288791, 1329.69579201, 1162.02041292, 1155.01082322,
          1165.80081172, 1174.83307617, 1180.5495586]])
    np.testing.assert_allclose(grd, grdt, atol=1e-2)

    grd = sc.z2pt5_cb14_cal
    grdt = np.array(
        [[1.17466154, 1.26861168, 1.1676564, 1.16797461, 1.16859058,
          1.16251358, 1.15315136],
         [1.17605704, 1.17803908, 1.16973403, 1.17649629, 1.18992131,
          1.18455663, 1.17590996],
         [1.16996381, 1.16401074, 1.166234, 1.17627836, 1.230198,
          1.32873844, 1.37243028],
         [1.15606906, 1.15025389, 1.15736893, 1.2809454, 1.19818264,
          1.19051805, 1.17596823],
         [1.49705669, 1.25107322, 1.31920733, 1.21901552, 1.16970542,
          1.1735483, 1.17388652],
         [1.54123541, 1.48388699, 1.17680997, 1.15783457, 1.16605299,
          1.17122878, 1.17355655],
         [1.57278236, 1.29292406, 1.1631469, 1.1579883, 1.16593269,
          1.17260055, 1.17683042]])
    np.testing.assert_allclose(grd, grdt, atol=1e-2)

    lldict = {'lats': np.array([34.1, 34.111]),
              'lons': np.array([-118.2, -118.222])}
    scsamp = mysite.getSitesContext(lldict)

    vs30 = scsamp.vs30
    vs30t = np.array([395.53771973, 428.14602661])
    np.testing.assert_allclose(vs30, vs30t)

    lldict = {'lats': np.array([34.1, 34.111]),
              'lons': np.array([-118.2, -118.222])}
    scrsamp = mysite.getSitesContext(lldict, rock_vs30=760)
    vs30 = scrsamp.vs30
    vs30t = np.array([760., 760.])
    np.testing.assert_allclose(vs30, vs30t)

    lats = np.array([34.1, 34.111])
    lons = np.array([-118.2, -118.222])
    lldict = {'lats': lats, 'lons': lons}
    scsamp = mysite.getSitesContext(lldict)
    grd = scsamp.vs30measured
    grdt = np.array([False, False], dtype=bool)
    np.testing.assert_allclose(grd, grdt)

    with pytest.raises(Exception) as e:  # noqa
        scsamp = mysite.getSitesContextFromLatLon(
            np.array([34.1, 34.111, 34.5]),
            np.array([-118.2, -118.222]))

    mysite = Sites.fromCenter(cx, cy, xspan, yspan, dx, dy,
                              vs30File=None, padding=True,
                              resample=False)
    grd = mysite.getVs30Grid().getData()
    grd_target = np.array(
        [[686., 686., 686., 686., 686., 686.],
         [686., 686., 686., 686., 686., 686.],
         [686., 686., 686., 686., 686., 686.],
         [686., 686., 686., 686., 686., 686.],
         [686., 686., 686., 686., 686., 686.],
         [686., 686., 686., 686., 686., 686.]])
    np.testing.assert_allclose(grd, grd_target)

    xmin = -118.2
    xmax = -118.12
    ymin = 34.05
    ymax = 34.1
    dx = 0.0083
    dy = 0.0083
    mysite = Sites.fromBounds(xmin, xmax, ymin, ymax, dx, dy,
                              vs30File=vs30file, padding=False,
                              resample=False)
    grd = mysite.getVs30Grid().getData()
    grd_target = np.array([[428.2287, 427.0027, 426.89514, 425.62024, 419.60953, 411.79617,
                            407.5351, 406.22122, 405.31622],
                           [429.40097, 427.7418, 427.00006, 419.50626, 411.1083, 407.489,
                            406.53305, 406.5966, 406.24887],
                           [429.4397, 427.30435, 425.96152, 426.15857, 427.56122, 397.67102,
                            399.21054, 404.54968, 407.18515],
                           [428.8627, 425.99606, 423.56927, 423.59836, 425.92758, 408.44885,
                            406.5581, 409.06946, 413.7521],
                           [427.91104, 424.53796, 419.47485, 418.17795, 424.14066, 428.57913,
                            432.953, 427.7773, 431.46524],
                           [423.15558, 424.48355, 419.27658, 418.6021, 423.86722, 428.06177,
                            432.4209, 438.54446, 448.37238]], dtype=np.float32)
    np.testing.assert_allclose(grd, grd_target, atol=1e-2)

    mysite = Sites.fromBounds(xmin, xmax, ymin, ymax, dx, dy,
                              vs30File=None, padding=False,
                              resample=False)
    grd = mysite.getVs30Grid().getData()
    grd_target = np.array(
        [[686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.],
         [686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.],
         [686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.],
         [686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.],
         [686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.],
         [686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.],
         [686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.],
         [686., 686., 686., 686., 686., 686., 686., 686., 686.,
          686., 686.]])
    np.testing.assert_allclose(grd, grd_target)

    # Test bounds
    nx, ny = mysite.getNxNy()
    assert nx == len(grd_target[0])
    assert ny == len(grd_target)

    # Check exception for invalid lat/lon shapes
    with pytest.raises(Exception):
        sample_dict = {'lats': np.array(['0', '0']),
                       'lons': np.array(['0'])}
        sc = mysite.getSitesContext(lldict=sample_dict, rock_vs30=None)

    # Check invalid file
    with pytest.raises(Exception):
        vs30file = os.path.join(homedir, 'sites_data/Wrong.grd')
        mysite = Sites.fromCenter(cx, cy, xspan, yspan, dx, dy,
                                  vs30File=vs30file, padding=True,
                                  resample=False)


if __name__ == '__main__':
    test_depthpars()
    test_sites()
