
import shakemap.grind.conversions.imt.bommer_alarcon_2006 as ba06


def test_bommeralarcon2006():
    # Inputs
    PGVin = 10
    PSA05in = 0.1

    # BommerAlarcon2006
    PSA05out = ba06.BommerAlarcon2006.pgv2psa05(PGVin)
    PGVout = ba06.BommerAlarcon2006.psa052pgv(PSA05in)
    vfact = ba06.B∑ommerAlarcon2006().getVfact()

    assert abs(PSA05out - 0.2038735983690112) < 0.0001
    assert abs(PGVout - 4.905000) < 0.001
    assert abs(vfact - 49.050000) < 0.001
