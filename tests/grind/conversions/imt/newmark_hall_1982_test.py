
import shakemap.grind.conversions.imt.newmark_hall_1982 as nh82


def test_newmarkhall1982():
    # Inputs
    PGVin = 10
    PSA10in = 0.1

    # NewmarkHall1982
    PSA10out = nh82.NewmarkHall1982.pgv2psa10(PGVin)
    PGVout = nh82.NewmarkHall1982.psa102pgv(PSA10in)
    vfact = nh82.NewmarkHall1982.getVfact()

    assert abs(PSA10out - 0.1056348) < 0.0001
    assert abs(PGVout - 9.46658) < 0.001
    assert abs(vfact - 94.6658) < 0.001
