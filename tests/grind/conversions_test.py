
import shakemap.grind.conversions as cv


def test_newmarkhall1982():
    # Inputs
    PGVin = 10
    PSA10in = 0.1

    # NewmarkHall1982
    PSA10out = cv.NewmarkHall1982.pgv2psa10(PGVin)
    PGVout = cv.NewmarkHall1982.psa102pgv(PSA10in)
    vfact = cv.NewmarkHall1982().getVfact()
    
    assert abs(PSA10out - 0.1056348) < 0.0001
    assert abs(PGVout - 9.46658) < 0.001
    assert abs(vfact - 94.6658) < 0.001


def test_bommeralarcon2006():
    # Inputs
    PGVin = 10
    PSA05in = 0.1

    # BommerAlarcon2006
    PSA05out = cv.BommerAlarcon2006.pgv2psa05(PGVin)
    PGVout = cv.BommerAlarcon2006.psa052pgv(PSA05in)
    vfact = cv.BommerAlarcon2006().getVfact()

    assert abs(PSA05out - 0.2038735983690112) < 0.0001
    assert abs(PGVout - 4.905000) < 0.001
    assert abs(vfact - 49.050000) < 0.001
