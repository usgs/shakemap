
from shakemap.utils.misc import get_command_output
from shakemap.utils.timeutils import ShakeDateTime


def test_get_command_output():
    cmd = 'ls *.py'
    rc, so, se = get_command_output(cmd)
    assert rc == True

    cmd = 'ls asdf'
    rc, so, se = get_command_output(cmd)
    assert rc == False


def test_timeutils():
    time = 1470152355.244582
    sdt = ShakeDateTime.utcfromtimestamp(time)
    assert sdt == ShakeDateTime(2016, 8, 2, 15, 39, 15, 244582)
