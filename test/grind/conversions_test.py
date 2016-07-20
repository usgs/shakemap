
import shakemap.grind.conversions as cv

# Inputs
PGVin = 10
PSA10in = 0.1

# NewmarkHall1982
PSA10out = cv.NewmarkHall1982.pgv2psa10(PGVin)
PGVout = cv.NewmarkHall1982.psa102pgv(PSA10in)

assert abs(PSA10out - 0.1056348) < 0.0001
assert abs(PGVout - 9.46658) < 0.001

# BommerAlarcon2006
PSA10out = cv.BommerAlarcon2006.pgv2psa10(PGVin)
PGVout = cv.BommerAlarcon2006.psa102pgv(PSA10in)

assert abs(PSA10out - 0.05096839) < 0.0001
assert abs(PGVout - 19.6200) < 0.001
