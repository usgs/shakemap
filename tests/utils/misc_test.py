
from shakemap.utils.misc import getCommandOutput
from shakemap.utils.timeutils import ShakeDateTime


def test_getCommandOutput():
    cmd = 'ls *.py'
    rc, so, se = getCommandOutput(cmd)
    assert rc == True

    cmd = 'ls asdf'
    rc, so, se = getCommandOutput(cmd)
    assert rc == False


def test_timeutils():
    time = 1470152355.244582
    sdt = ShakeDateTime.utcfromtimestamp(time)
    assert sdt == ShakeDateTime(2016, 8, 2, 15, 39, 15, 244582)
