#!/usr/bin/env python

# stdlib imports
import os.path
import sys

# third party imports
import numpy as np
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA
import pytest

# local imports
from shakelib.conversions.imc.boore_kishida_2017 import BooreKishida2017


do_test = True

amps_in = np.log(np.array([0.05, 0.1, 0.2, 0.4, 0.8, 1.6]))
sigmas_in = np.array([0.5, 0.55, 0.6, 0.65, 0.61, 0.7])
rrup_in = np.array([100.0, 50.0, 25.0, 12.0, 6.0, 1.0])
imc_in = [const.IMC.RotD50,
          const.IMC.RotD50,
          const.IMC.RotD100,
          const.IMC.RotD50,
          const.IMC.RotD100,
          const.IMC.GMRotI50,
          const.IMC.AVERAGE_HORIZONTAL]
imc_out = [const.IMC.GMRotI50,
           const.IMC.AVERAGE_HORIZONTAL,
           const.IMC.RotD50,
           const.IMC.GREATER_OF_TWO_HORIZONTAL,
           const.IMC.GREATER_OF_TWO_HORIZONTAL,
           const.IMC.GREATER_OF_TWO_HORIZONTAL,
           const.IMC.GREATER_OF_TWO_HORIZONTAL]
mags_in = np.array([5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0])
imt_in = [PGA(), PGV(), SA(0.3), SA(1.0), SA(3.05)]


amps_target = np.array(
    [[-2.99751949, -2.30634141, -1.61516334, -0.92410122, -0.23292315,
        0.45513398],
     [-3.00428338, -2.31469074, -1.6250981, -0.9357148, -0.24612216,
        0.43783666],
     [-3.01181065, -2.32300028, -1.6341899, -0.94563493, -0.25682456,
        0.42511215],
     [-3.01695413, -2.3278757, -1.63879727, -0.94995846, -0.26088004,
        0.42174957],
     [-3.01848387, -2.32967754, -1.64087122, -0.95232054, -0.26351421,
        0.41841202],
     [-3.00567931, -2.31596465, -1.62624999, -0.93673749, -0.24702283,
        0.43725142],
     [-3.01424793, -2.32494632, -1.63564472, -0.94656959, -0.25726798,
        0.42593853],
     [-3.02035418, -2.33256737, -1.64478057, -0.95730945, -0.26952264,
        0.40976817],
     [-3.02735109, -2.33901171, -1.65067234, -0.96261611, -0.27427674,
        0.40644244],
     [-3.03288451, -2.34510743, -1.65733036, -0.96986955, -0.28209248,
        0.39717318],
     [-3.17123427, -2.48581923, -1.8004042, -1.11544454, -0.4300295,
        0.24313038],
     [-3.1919086, -2.50757532, -1.82324205, -1.13942786, -0.45509458,
        0.21526899],
     [-3.19063745, -2.50581691, -1.82099636, -1.1366662, -0.45184565,
        0.21977749],
     [-3.20465738, -2.51876928, -1.83288119, -1.14742061, -0.46153252,
        0.21285019],
     [-3.21345032, -2.52681264, -1.84017497, -1.15392066, -0.46728298,
        0.20903738],
     [-2.89494534, -2.19513344, -1.49532154, -0.79511713, -0.09530523,
        0.61507],
     [-2.87776794, -2.17755962, -1.4773513, -0.77672712, -0.0765188,
        0.63488116],
     [-2.87967312, -2.17953214, -1.47939117, -0.7788383, -0.07869733,
        0.63252855],
     [-2.86729023, -2.16894327, -1.47059632, -0.77194313, -0.07359618,
        0.63299222],
     [-2.86010338, -2.16206455, -1.46402571, -0.76569879, -0.06765996,
        0.63813196],
     [-3.06695508, -2.37487533, -1.68279557, -0.99077868, -0.29869892,
        0.391689],
     [-3.07267199, -2.38127757, -1.68988315, -0.99859196, -0.30719755,
        0.3814188],
     [-3.06904939, -2.37723505, -1.68542071, -0.99368487, -0.30187053,
        0.38783131],
     [-3.07349526, -2.38240739, -1.69131952, -1.00035294, -0.30926507,
        0.37855886],
     [-3.0759971, -2.38446777, -1.69293844, -1.00150439, -0.30997506,
        0.37899003],
     [-2.90921793, -2.20743693, -1.50565593, -0.80336645, -0.10158544,
        0.61387984],
     [-2.86706637, -2.16330351, -1.45954064, -0.75515259, -0.05138972,
        0.6691986],
     [-2.88325616, -2.17877839, -1.47430061, -0.76915553, -0.06467775,
        0.6577586],
     [-2.85924461, -2.15682891, -1.4544132, -0.75145164, -0.04903593,
        0.66807004],
     [-2.839502, -2.13712231, -1.43474262, -0.7318192, -0.02943951,
        0.68757335],
     [-2.90521964, -2.20197522, -1.4987308, -0.79489171, -0.09164729,
        0.62760088],
     [-2.85536289, -2.151309, -1.44725511, -0.74255887, -0.03850498,
        0.68283565],
     [-2.88561201, -2.18011066, -1.47460931, -0.76838038, -0.06287903,
        0.66220321],
     [-2.85903848, -2.15588372, -1.45272896, -0.74898482, -0.04583006,
        0.67318634],
     [-2.83367318, -2.13026424, -1.4268553, -0.72284201, -0.01943307,
        0.70024037]])

sigs_target = np.array(
    [[0.5021999,  0.55200067,  0.60183448,  0.65169375,  0.6118045,
        0.70157305],
     [0.50366136,  0.55333061,  0.60305453,  0.65282062,  0.6130047,
        0.70261993],
     [0.50360749,  0.55328158,  0.60300954,  0.65277906,  0.61296044,
        0.70258132],
     [0.50376195,  0.55342217,  0.60313854,  0.65289823,  0.61308735,
        0.70269204],
     [0.50401437,  0.55365195,  0.60334938,  0.65309301,  0.61329477,
        0.70287302],
     [0.50692913,  0.5563067,  0.60578638,  0.65534505,  0.61569241,
        0.70496606],
     [0.50731408,  0.55665751,  0.60610855,  0.65564287,  0.6160094,
        0.70524292],
     [0.50833411,  0.55758728,  0.60696258,  0.65643246,  0.61684971,
        0.70597703],
     [0.50842484,  0.55766999,  0.60703856,  0.65650271,  0.61692448,
        0.70604236],
     [0.50734364,  0.55668445,  0.6061333,  0.65566575,  0.61603374,
        0.70526419],
     [0.50695541,  0.55633065,  0.60580838,  0.65536539,  0.61571405,
        0.70498496],
     [0.50726646,  0.55661411,  0.6060687,  0.65560603,  0.61597018,
        0.70520867],
     [0.50691702,  0.55629567,  0.60577625,  0.65533569,  0.61568244,
        0.70495735],
     [0.50696301,  0.55633757,  0.60581473,  0.65537126,  0.6157203,
        0.70499042],
     [0.50712286,  0.55648324,  0.60594851,  0.65549493,  0.61585193,
        0.70510538],
     [0.50873603,  0.55795372,  0.60729923,  0.65674375,  0.61718097,
        0.70626649],
     [0.50918036,  0.55835888,  0.60767149,  0.657088,  0.61754728,
        0.70658661],
     [0.50893265,  0.558133,  0.60746394,  0.65689607,  0.61734305,
        0.70640813],
     [0.50913428,  0.55831685,  0.60763287,  0.65705229,  0.61750928,
        0.7065534],
     [0.50939642,  0.55855592,  0.60785254,  0.65725544,  0.61772544,
        0.70674233],
     [0.50494378,  0.55449817,  0.604126,  0.65381054,  0.61405881,
        0.70353978],
     [0.50554264,  0.55504357,  0.60462663,  0.65427316,  0.61455135,
        0.70396972],
     [0.50540885,  0.55492171,  0.60451477,  0.65416978,  0.6144413,
        0.70387364],
     [0.50572465,  0.55520934,  0.60477882,  0.6544138,  0.61470108,
        0.70410043],
     [0.50592042,  0.55538768,  0.60494254,  0.6545651,  0.61486216,
        0.70424106],
     [0.51013851,  0.55923278,  0.60847457,  0.65783075,  0.61833753,
        0.70727739],
     [0.51241364,  0.56130895,  0.61038327,  0.65959665,  0.62021588,
        0.70892012],
     [0.51294738,  0.56179624,  0.61083142,  0.66001138,  0.62065693,
        0.70930601],
     [0.51335585,  0.56216922,  0.61117447,  0.66032888,  0.62099455,
        0.70960146],
     [0.51423859,  0.56297542,  0.61191611,  0.66101537,  0.62172448,
        0.71024033],
     [0.51337003,  0.56218217,  0.61118638,  0.6603399,  0.62100627,
        0.70961172],
     [0.51556409,  0.56418643,  0.61303045,  0.66204708,  0.62282127,
        0.71120063],
     [0.51854404,  0.56691086,  0.61553873,  0.66437032,  0.62529027,
        0.71336381],
     [0.51950106,  0.56778636,  0.61634515,  0.66511754,  0.62608414,
        0.71405977],
     [0.51892048,  0.5672552,  0.61585588,  0.66466418,  0.62560248,
        0.71363749]])


def test_bk17():
    amps_out = np.empty([0, 6])
    sigs_out = np.empty([0, 6])
    for i, imc in enumerate(imc_in):
        bk17 = BooreKishida2017(imc, imc_out[i])
        for imt in imt_in:
            tmp = bk17.convertAmps(imt, amps_in, rrup_in, mags_in[i])
            amps_out = np.vstack((amps_out, tmp))
            tmp = bk17.convertSigmas(imt, sigmas_in)
            sigs_out = np.vstack((sigs_out, tmp))

    if do_test is True:
        np.testing.assert_allclose(amps_out, amps_target, atol=1e-5)
        np.testing.assert_allclose(sigs_out, sigs_target, atol=1e-5)

        # Do a round trip
        bk17 = BooreKishida2017(imc_in[0], imc_out[0])
        tmp = bk17.convertAmps(imt_in[0], amps_in, rrup_in, mags_in[0])
        bk17 = BooreKishida2017(imc_out[0], imc_in[0])
        orig_amps = bk17.convertAmps(imt_in[0], tmp, rrup_in, mags_in[0])
        np.testing.assert_allclose(amps_in, orig_amps, atol=1e-5)
    else:
        print(repr(amps_out))
        print(repr(sigs_out))

    # Test that an invalid/unknown parameter is changed to AVERAGE_HORIZONTAL
    bk17 = BooreKishida2017('wrong', imc_out[0])
    assert bk17.imc_in == 'Average horizontal'
    assert bk17.imc_out == imc_out[0]
    bk17 = BooreKishida2017(imc_out[0], 'wrong')
    assert bk17.imc_in == imc_out[0]
    assert bk17.imc_out == 'Average horizontal'
    bk17 = BooreKishida2017('wrong', 'wrong')
    assert bk17.imc_in == 'Average horizontal'
    assert bk17.imc_out == 'Average horizontal'

    # Test that the correct input/output imc returns the right path
    bk17 = BooreKishida2017('Average horizontal',
            'Average Horizontal (RotD50)')
    assert len(bk17.path) == 2
    assert bk17.path[0] == 'Average horizontal'
    assert bk17.path[-1] == 'Average Horizontal (RotD50)'
    bk17 = BooreKishida2017('Average horizontal',
            'Horizontal Maximum Direction (RotD100)')
    assert len(bk17.path) == 3
    assert bk17.path[0] == 'Average horizontal'
    assert bk17.path[-1] == 'Horizontal Maximum Direction (RotD100)'
    if bk17.path[1] == 'Greater of two horizontal':
        correct = True
    elif bk17.path[1] == 'Average Horizontal (RotD50)':
        correct = True
    else:
        correct = False
    assert correct == True
    bk17 = BooreKishida2017('Horizontal',
            'Random horizontal')
    assert len(bk17.path) == 3
    assert bk17.path[0] == 'Horizontal'
    assert bk17.path[-1] == 'Random horizontal'
    if bk17.path[1] == 'Greater of two horizontal':
        correct = True
    elif bk17.path[1] == 'Average Horizontal (RotD50)':
        correct = True
    else:
        correct = False
    assert correct == True

    # Check output amps for chained conversion is the same as two seperate
    # conversions as suggested in the original class docstring
    bk17 = BooreKishida2017('Horizontal Maximum Direction (RotD100)',
            'Average Horizontal (RotD50)')
    mid = bk17.convertAmps(imt_in[0], amps_in, rrup_in, mags_in[0])
    bk17 = BooreKishida2017('Average Horizontal (RotD50)',
            'Average Horizontal (GMRotI50)')
    last = bk17.convertAmps(imt_in[0], mid, rrup_in, mags_in[0])
    bk17 = BooreKishida2017('Horizontal Maximum Direction (RotD100)',
            'Average Horizontal (GMRotI50)')
    full = bk17.convertAmps(imt_in[0],amps_in, rrup_in, mags_in[0])
    np.testing.assert_allclose(full, last, atol=1e-5)

    # Test exception for missing magnitude/rupture
    with pytest.raises(ValueError) as e:
        bk17.convertAmps(imt_in[0],amps_in, None, mags_in[0])
    with pytest.raises(ValueError) as e:
        bk17.convertAmps(imt_in[0],amps_in, rrup_in, None)
    with pytest.raises(ValueError) as e:
        bk17.convertAmps(imt_in[0],amps_in, None, None)
    # Test exception for unknown imt
    with pytest.raises(ValueError) as e:
        bk17.convertAmpsOnce('wrong', [10.0], rrup_in, mags_in[0])
    # Test amp conversion for unknown imc
    bk17 = BooreKishida2017('Average Horizontal',
            'Average Horizontal')
    a = [-2.99573227, -2.30258509, -1.60943791, -0.91629073,
            -0.22314355, 0.47000363]
    mag = 8.0
    # Average => Average
    target = bk17.convertAmpsOnce(PGA(), a, rrup_in, mag)
    # Average => wrong. The incorrect imc should use GMAR and be the same as
    # the above result
    bk17.imc_out = 'wrong'
    result = bk17.convertAmpsOnce(PGA(), a, rrup_in, mag)
    np.testing.assert_allclose(target, result, atol=1e-5)


if __name__ == '__main__':
    test_bk17()
