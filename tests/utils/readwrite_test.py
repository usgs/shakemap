#!/usr/bin/env python

import os.path

import numpy as np

from shakemap.utils.readwrite import read_nshmp_fault_xml
from shakemap.utils.readwrite import read_nshmp_grid_xml

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '../..'))
datdir = os.path.join(homedir, 'readwrite_data')

def test_read_nshmp_fault_xml():
    file = os.path.join(datdir, 'sub0_ch_mid.xml')
    sub = read_nshmp_fault_xml(file)

    target = np.array(
        [ 24.01755,  14.83582,  12.12714,  13.66332,  15.27711,  15.98822,
          19.36195,  17.45862,  20.20983,  22.92826,  22.9094 ,  23.94863,
          27.91269,  29.32002,  27.86888,  27.84144,  25.03354,  22.97471,
          23.69129])
    np.testing.assert_allclose(sub[0]['Geometry']['LowerTrace']['dep'], target)

    target = np.array(
      [ 40.347  ,  41.21782,  42.1117 ,  42.97969,  43.7    ,  43.86274,
        44.74235,  45.     ,  45.48932,  46.3    ,  46.36419,  46.7786 ,
        47.23766,  47.70893,  48.15165,  48.55325,  48.86098,  49.21444,
        49.73665]
    )
    np.testing.assert_allclose(sub[0]['Geometry']['LowerTrace']['lat'], target)

    target = np.array(
      [-123.829  , -124.17259, -124.38969, -124.51445, -124.50873,
       -124.49237, -124.35632, -124.4883 , -124.32961, -124.13677,
       -124.13849, -124.11635, -123.89909, -123.91646, -124.36144,
       -124.80386, -125.43884, -126.17539, -126.82846]
    )
    np.testing.assert_allclose(sub[0]['Geometry']['LowerTrace']['lon'], target)

    target = np.array(
      [ 6.74507,  5.     ,  5.     ,  5.     ,  5.     ,  5.     ,
        5.     ,  5.     ,  5.     ,  5.     ,  5.     ,  5.     ,
        5.     ,  5.     ,  5.     ,  5.     ,  5.     ,  5.     ,  5.     ]
    )
    np.testing.assert_allclose(sub[0]['Geometry']['Trace']['dep'], target)

    target = np.array(
      [ 40.35467,  41.2136 ,  42.09936,  42.95745,  43.7    ,  43.84118,
        44.71223,  45.     ,  45.45161,  46.3    ,  46.32987,  46.6294 ,
        46.94307,  47.27867,  47.66285,  47.99196,  48.27786,  48.68892,
        49.2529 ]
    )
    np.testing.assert_allclose(sub[0]['Geometry']['Trace']['lat'], target)

    target = np.array(
      [-125.09899, -125.08518, -125.04685, -125.28921, -125.41145,
       -125.42693, -125.49976, -125.52694, -125.57416, -125.72885,
       -125.73574, -125.81562, -125.92025, -126.04937, -126.23375,
       -126.439  , -126.6781 , -127.0801 , -127.64473]
    )
    np.testing.assert_allclose(sub[0]['Geometry']['Trace']['lon'], target)

    assert sub[0]['IncrementalMfd'][0]['floats'] == 'false'
    assert sub[0]['IncrementalMfd'][0]['m'] == '9.12'
    assert sub[0]['IncrementalMfd'][0]['rate'] == '0.0019'
    assert sub[0]['IncrementalMfd'][0]['type'] == 'SINGLE'
    assert sub[0]['IncrementalMfd'][0]['weight'] == '0.333'
    assert sub[0]['IncrementalMfd'][1]['floats'] == 'false'
    assert sub[0]['IncrementalMfd'][1]['m'] == '8.69'
    assert sub[0]['IncrementalMfd'][1]['rate'] == '0.0019'
    assert sub[0]['IncrementalMfd'][1]['type'] == 'SINGLE'
    assert sub[0]['IncrementalMfd'][1]['weight'] == '0.334'
    assert sub[0]['IncrementalMfd'][2]['floats'] == 'false'
    assert sub[0]['IncrementalMfd'][2]['m'] == '8.82'
    assert sub[0]['IncrementalMfd'][2]['rate'] == '0.0019'
    assert sub[0]['IncrementalMfd'][2]['type'] == 'SINGLE'
    assert sub[0]['IncrementalMfd'][2]['weight'] == '0.333'

def test_read_nshmp_grid_xml():
    file = os.path.join(datdir, 'USGS_grid.xml')
    sub = read_nshmp_grid_xml(file)

    # Test Settings
    mfds = sub['Settings']['DefaultMfds']
    srcp = sub['Settings']['SourceProperties']
    assert [m['a'] for m in mfds] == ['0.0', '0.0', '0.0', '0.0']
    assert [m['b'] for m in mfds] == ['1.0', '1.0', '1.0', '1.0']
    assert [m['dMag'] for m in mfds] == ['0.1', '0.1', '0.1', '0.1']
    assert [m['mMax'] for m in mfds] == ['6.45', '6.95', '7.45', '7.95']
    assert [m['mMin'] for m in mfds] == ['4.75', '4.75', '4.75', '4.75']
    assert [m['type'] for m in mfds] == ['GR', 'GR', 'GR', 'GR']
    assert [m['weight'] for m in mfds] == ['0.2', '0.5', '0.2', '0.1']
    assert srcp['focalMechMap'] == '[STRIKE_SLIP:1.0,NORMAL:0.0,REVERSE:0.0]'
    assert srcp['magDepthMap'] == '[10.0::[5.0:1.0]]'
    assert srcp['maxDepth'] == '22.0'
    assert srcp['ruptureScaling'] == 'NSHM_SOMERVILLE'
    assert srcp['strike'] == 'NaN'

    # Test nodes
    nodes = sub['Nodes']
    z = np.array([n['dep'] for n in nodes])
    ztarget = np.array(
      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
    np.testing.assert_allclose(z, ztarget)

    lat = np.array([n['lat'] for n in nodes])
    lattarget = np.array(
      [ 34.6,  34.6,  34.6,  34.6,  34.6,  34.6,  34.7,  34.7,  34.7,
        34.7,  34.7,  34.7,  34.7,  39.5])
    np.testing.assert_allclose(lat, lattarget)

    lon = np.array([n['lon'] for n in nodes])
    lontarget = np.array(
      [-109.8, -109.7, -109.6, -109.5, -109.4, -109.3, -110. , -109.9,
       -109.8, -109.7, -109.6, -109.5, -109.4, -110.7])
    np.testing.assert_allclose(lon, lontarget)

    a = np.array([float(n['a']) for n in nodes])
    atarget = np.array(
      [ 0.06516867,  0.06632903,  0.06741378,  0.06840201,  0.0692732 ,
        0.07000767,  0.06527287,  0.06659666,  0.06788012,  0.06910356,
        0.07024578,  0.07128526,  0.07220085,  0.27957186])
    np.testing.assert_allclose(lon, lontarget)

    t = np.array([n['type'] for n in nodes])
    ttarget = ['GR', 'GR', 'GR', 'GR', 'GR', 'GR', 'GR', 'GR', 'GR', 'GR', 'GR',
       'GR', 'GR', 'GR']
    assert all(t == ttarget)


def test_read_nshmp_rlme_xml():
    file = os.path.join(datdir, 'Charlevoix Seismic Zone.xml')
    sub = read_nshmp_grid_xml(file)

    # Test Settings
    mfds = sub['Settings']['DefaultMfds']
    srcp = sub['Settings']['SourceProperties']
    assert [m['floats'] for m in mfds] == ['false', 'false', 'false', 'false']
    assert [m['m'] for m in mfds] == ['6.75', '7.0', '7.25', '7.5']
    assert [m['rate'] for m in mfds] == ['0.0', '0.0', '0.0', '0.0']
    assert [m['type'] for m in mfds] == ['SINGLE', 'SINGLE', 'SINGLE', 'SINGLE']
    assert [m['weight'] for m in mfds] == ['0.2', '0.5', '0.2', '0.1']
    assert srcp['focalMechMap'] == '[STRIKE_SLIP:1.0,NORMAL:0.0,REVERSE:0.0]'
    assert srcp['magDepthMap'] == '[10.0::[5.0:1.0]]'
    assert srcp['maxDepth'] == '22.0'
    assert srcp['ruptureScaling'] == 'NSHM_POINT_WC94_LENGTH'
    assert srcp['strike'] == 'NaN'

    # Test Nodes
    nodes = sub['Nodes']
    z = np.array([n['dep'] for n in nodes])
    ztarget = np.array(
      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.])
    np.testing.assert_allclose(z, ztarget)

    lat = np.array([n['lat'] for n in nodes])
    lattarget = np.array(
      [ 47. ,  47.1,  47.1,  47.1,  47.2,  47.2,  47.2,  47.2,  47.2,
        47.3,  47.3,  47.3,  47.3,  47.3,  47.3,  47.3,  47.4,  47.4,
        47.4,  47.4,  47.4,  47.4,  47.4,  47.4,  47.4,  47.5,  47.5,
        47.5,  47.5,  47.5,  47.5,  47.5,  47.5,  47.5,  47.5,  47.6,
        47.6,  47.6,  47.6,  47.6,  47.6,  47.6,  47.6,  47.6,  47.6,
        47.7,  47.7,  47.7,  47.7,  47.7,  47.7,  47.7,  47.7,  47.7,
        47.7,  47.8,  47.8,  47.8,  47.8,  47.8,  47.8,  47.8,  47.8,
        47.8,  47.9,  47.9,  47.9])
    np.testing.assert_allclose(lat, lattarget)

    lon = np.array([n['lon'] for n in nodes])
    lontarget = np.array(
      [-70.4, -70.5, -70.4, -70.3, -70.6, -70.5, -70.4, -70.3, -70.2,
       -70.7, -70.6, -70.5, -70.4, -70.3, -70.2, -70.1, -70.7, -70.6,
       -70.5, -70.4, -70.3, -70.2, -70.1, -70. , -69.9, -70.7, -70.6,
       -70.5, -70.4, -70.3, -70.2, -70.1, -70. , -69.9, -69.8, -70.6,
       -70.5, -70.4, -70.3, -70.2, -70.1, -70. , -69.9, -69.8, -69.7,
       -70.5, -70.4, -70.3, -70.2, -70.1, -70. , -69.9, -69.8, -69.7,
       -69.6, -70.4, -70.3, -70.2, -70.1, -70. , -69.9, -69.8, -69.7,
       -69.6, -70.1, -70. , -69.9])
    np.testing.assert_allclose(lon, lontarget)

    rate = np.array([float(n['rate']) for n in nodes])
    ratetarget = np.array(
      [  2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05,   2.04776120e-05,   2.04776120e-05,
         2.04776120e-05])
    np.testing.assert_allclose(rate, ratetarget)

if __name__ == '__main__':

    test_read_nshmp_rlme_xml()
    test_read_nshmp_grid_xml()
    test_read_nshmp_fault_xml()
