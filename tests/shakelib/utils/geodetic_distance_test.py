#!/usr/bin/env python

import numpy as np

from shakelib.utils.distance import (geodetic_distance_fast, 
                                     geodetic_distance_haversine)


def test_distances():
    W = -102.0
    E = -100.0
    S = 62.0
    N = 64.0
    lonspan = E - W
    latspan = N - S
    lons1 = np.linspace(W, E, 4 * int(lonspan))
    lats1 = np.linspace(S, N, 4 * int(latspan))
    lons1_radians = np.radians(lons1).reshape((-1, 1))
    lats1_radians = np.radians(lats1).reshape((-1, 1))
    W = -105.0
    E = -104.0
    S = 66.0
    N = 67.0
    lonspan = E - W
    latspan = N - S
    lons2 = np.linspace(W, E, 4 * int(lonspan))
    lats2 = np.linspace(S, N, 4 * int(latspan))
    lons2_radians = np.radians(lons2).reshape((1, -1))
    lats2_radians = np.radians(lats2).reshape((1, -1))

    d1 = geodetic_distance_fast(lons1_radians, lats1_radians, 
                                lons2_radians, lats2_radians)

    d1_test = np.array(
      [[ 468.20225686,  498.86822304,  530.93902911,  564.15794707],
       [ 442.68140081,  472.03948363,  503.04222914,  535.38587971],
       [ 418.3695637 ,  446.19575113,  475.95071332,  507.27603453],
       [ 395.47553776,  421.50488633,  449.79749798,  479.9330591 ],
       [ 374.2446911 ,  398.16739986,  424.74280617,  453.4833895 ],
       [ 354.95998908,  376.42039946,  400.97881782,  428.07972608],
       [ 337.93853762,  356.53965837,  378.73409841,  403.90591865],
       [ 323.52096317,  338.83780046,  358.27658802,  381.18178394]])

    assert np.allclose(d1, d1_test)

    d2 = geodetic_distance_haversine(lons1_radians, lats1_radians, 
                                     lons2_radians, lats2_radians)

    d2_test = np.array(
      [[ 468.07174603,  498.75566141,  530.84565791,  564.08409095],
       [ 442.53667791,  471.9109949 ,  502.93199853,  535.29501382],
       [ 418.21264626,  446.05269576,  475.82431308,  507.16818137],
       [ 395.30897492,  421.34919191,  449.65619258,  479.80881082],
       [ 374.07149802,  398.00153662,  424.58843141,  453.34391388],
       [ 354.78352531,  376.24732442,  400.8137616 ,  427.92676557],
       [ 337.76230336,  356.36270821,  378.56125359,  403.74177579],
       [ 323.34828811,  338.66050755,  358.09925893,  381.00928211]])

    assert np.allclose(d2, d2_test)


if __name__ == '__main__':
    test_distances()
