#!/usr/bin/env python

import os
import shutil
import subprocess
import tempfile


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))

########################################################################
# Test sm_profile
########################################################################


def test_profile():

    with tempfile.TemporaryDirectory() as tmpdir:
        pfile = os.path.join(str(tmpdir), 'profile.conf')
        program = os.path.join(shakedir, 'bin', 'sm_profile')
        #
        # Do some stuff that won't work until we create a profile
        #
        cp = subprocess.run([program, '-f', pfile], shell=False)
        assert cp.returncode
        cp = subprocess.run([program, '-f', pfile, '-a'], shell=False)
        assert cp.returncode
        cp = subprocess.run([program, '-f', pfile, '-l'], shell=False)
        assert cp.returncode
        cp = subprocess.run([program, '-f', pfile, '-s', 'dummy'], shell=False)
        assert cp.returncode
        cp = subprocess.run([program, '-f', pfile, '-d', 'dummy'], shell=False)
        assert cp.returncode
        #
        # Now create a profile and try everything again
        # Presumably we should test the outputs and results of these actions,
        # but that's a lot of monkeying around for not much benefit
        #
        cp = subprocess.run([program, '-f', pfile, '-c',
                             'test_profile1', '-a', '-n'], shell=False)
        assert not cp.returncode
        # No args program run (should succeed)
        cp = subprocess.run([program, '-f', pfile], shell=False)
        assert not cp.returncode
        # -a by itself (should fail)
        cp = subprocess.run([program, '-f', pfile, '-a'], shell=False)
        assert cp.returncode
        # Create the same profile a second time (should fail
        cp = subprocess.run([program, '-f', pfile, '-c',
                             'test_profile1', '-a', '-n'], shell=False)
        assert cp.returncode
        # List should now work
        cp = subprocess.run([program, '-f', pfile, '-l'], shell=False)
        assert not cp.returncode
        # Make a couple more profiles
        cp = subprocess.run([program, '-f', pfile, '-c',
                             'test_profile2', '-a', '-n'], shell=False)
        assert not cp.returncode
        cp = subprocess.run([program, '-f', pfile, '-c',
                             'test_profile3', '-a', '-n'], shell=False)
        assert not cp.returncode
        # List should still work
        cp = subprocess.run([program, '-f', pfile, '-l'], shell=False)
        assert not cp.returncode
        # Delete a profile that doesn't exist (should fail)
        cp = subprocess.run([program, '-f', pfile, '-d',
                             'not_a_profile', '-a'], shell=False)
        assert cp.returncode
        # Delete a profile that we're not set to (should succeed)
        cp = subprocess.run([program, '-f', pfile, '-d',
                             'test_profile2', '-a'], shell=False)
        assert not cp.returncode
        # Switch to a profile that doesn't exist (should fail)
        cp = subprocess.run([program, '-f', pfile, '-s', 'dummy'], shell=False)
        assert cp.returncode
        # Switch to a profile that exists (should succeed)
        cp = subprocess.run([program, '-f', pfile, '-s',
                             'test_profile1'], shell=False)
        assert not cp.returncode
        # Delete profiles that we are set to (should succeed)
        cp = subprocess.run([program, '-f', pfile, '-d',
                             'test_profile1', '-a'], shell=False)
        assert not cp.returncode
        cp = subprocess.run([program, '-f', pfile, '-d',
                             'test_profile3', '-a'], shell=False)
        assert not cp.returncode
        # Create a profile and accept the defaults (should succeed)
        # Making the 'test' profile, so it will copy Northridge
        op = subprocess.Popen([program, '-f', pfile, '-c', 'test', '-n'],
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=False)
        op.communicate('\n\n'.encode('ascii'))
        assert not op.returncode
        # delete it, but change your mind (should succeed)
        op = subprocess.Popen([program, '-f', pfile, '-d', 'test', '-n'],
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=False)
        op.communicate('n\n'.encode('ascii'))
        assert not op.returncode
        # now delete it for real (should succeed)
        op = subprocess.Popen([program, '-f', pfile, '-d', 'test', '-n'],
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=False)
        op.communicate('\n'.encode('ascii'))
        assert not op.returncode
        # There should now be junk from 'test' and 'test_profileX' (where
        # X is 1, 2, 3) lying around. It will presumably get cleaned up
        # when we fall out of the context of the enclosing 'with' statement

        # Create a profile with directories we can't make (should fail)

        #
        # Note: some test systems (Jenkins, Azure) give you unlimited
        # powers, so the following tests (that should fail) actually succeed.
        # If we can make a directory at the root level, we have su powers,
        # so we just skip these tests.
        try:
            os.mkdir("/not_possible_xxyzggkfmtpz")
        except PermissionError:
            # Normal system privileges -- run the tests
            run_mkdir_tests = True
        except OSError:
            # Normal system privileges -- run the tests
            run_mkdir_tests = True
        else:
            # Godlike system privileges -- the tests will pass even if they
            # shouldn't, so don't run them since success is failure
            run_mkdir_tests = False
            # Should try to delete the directory we just made, but don't
            # want to bomb out if we can't for some bizarre reason.
            try:
                shutil.rmtree("/not_possible_xxyzggkfmtpz")
            except Exception:
                pass

        if run_mkdir_tests:
            op = subprocess.Popen([program, '-f', pfile, '-c',
                                   'junk_profile', '-n'],
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  shell=False)
            ipath = os.path.join(os.path.abspath(os.sep), 'zzx/yzzskx') + '\n'
            op.communicate((ipath + ipath + ipath).encode('ascii'))
            assert op.returncode

            op = subprocess.Popen([program, '-f', pfile, '-c',
                                   'junk_profile', '-n'],
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  shell=False)
            ipath = os.path.join(shakedir, 'junkprofile', 'install') + '\n'
            dpath = os.path.join(os.path.abspath(os.sep), 'zzxy/zzskx') + '\n'
            op.communicate((ipath + dpath + dpath + dpath).encode('ascii'))
            assert op.returncode


########################################################################
# main program
########################################################################
if __name__ == '__main__':

    test_profile()
