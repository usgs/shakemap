#!/usr/bin/env python

import os
import shutil
import subprocess

from shakemap.utils.config import get_config_paths

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))

########################################################################
# Test sm_profile
########################################################################
def test_profile():
    #
    # Rename existing profile(s)
    #
    smrc = os.path.join(os.path.expanduser('~'), '.shakemap')
    if os.path.isdir(smrc):
        os.rename(smrc, smrc + '_safe')
    #
    # Do some stuff that won't work until we create a profile
    #
    program = os.path.join(shakedir, 'sm_profile')
    cp = subprocess.run([program], shell=False)
    assert cp.returncode
    cp = subprocess.run([program, '-a'], shell=False)
    assert cp.returncode
    cp = subprocess.run([program, '-l'], shell=False)
    assert cp.returncode
    cp = subprocess.run([program, '-s', 'dummy'], shell=False)
    assert cp.returncode
    cp = subprocess.run([program, '-d', 'dummy'], shell=False)
    assert cp.returncode
    #
    # Now create a profile and try everything again
    # Presumably we should test the outputs and results of these actions,
    # but that's a lot of monkeying around for not much benefit
    #
    cp = subprocess.run([program, '-c', 'test_profile1', '-a'], shell=False)
    assert not cp.returncode
    # No args program run (should succeed)
    cp = subprocess.run([program], shell=False)
    assert not cp.returncode
    # -a by itself (should fail)
    cp = subprocess.run([program, '-a'], shell=False)
    assert cp.returncode
    # Create the same profile a second time
    cp = subprocess.run([program, '-c', 'test_profile1', '-a'], shell=False)
    assert cp.returncode
    # List should now work
    cp = subprocess.run([program, '-l'], shell=False)
    assert not cp.returncode
    # Make a couple more profiles
    cp = subprocess.run([program, '-c', 'test_profile2', '-a'], shell=False)
    assert not cp.returncode
    cp = subprocess.run([program, '-c', 'test_profile3', '-a'], shell=False)
    assert not cp.returncode
    # List should still work
    cp = subprocess.run([program, '-l'], shell=False)
    assert not cp.returncode
    # Delete a profile that doesn't exist (should fail)
    cp = subprocess.run([program, '-d', 'not_a_profile', '-a'], shell=False)
    assert cp.returncode
    # Delete a profile that we're not set to (should succeed)
    cp = subprocess.run([program, '-d', 'test_profile2', '-a'], shell=False)
    assert not cp.returncode
    # Switch to a profile that doesn't exist (should fail)
    cp = subprocess.run([program, '-s', 'dummy'], shell=False)
    assert cp.returncode
    # Switch to a profile that exists (should succeed)
    cp = subprocess.run([program, '-s', 'test_profile1'], shell=False)
    assert not cp.returncode
    # Delete profiles that we are set to (should succeed)
    cp = subprocess.run([program, '-d', 'test_profile1', '-a'], shell=False)
    assert not cp.returncode
    cp = subprocess.run([program, '-d', 'test_profile3', '-a'], shell=False)
    assert not cp.returncode
    # Create a profile and accept the defaults (should succeed)
    # Making the 'test' profile, so it will copy Northridge
    op = subprocess.Popen([program, '-c', 'test'], 
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=False)
    op.communicate('\n\n'.encode('ascii'))
    assert not op.returncode
    # delete it, but change your mind (should succeed)
    op = subprocess.Popen([program, '-d', 'test'], 
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=False)
    op.communicate('n\n'.encode('ascii'))
    assert not op.returncode
    # now delete it for real (should succeed)
    op = subprocess.Popen([program, '-d', 'test'], 
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=False)
    op.communicate('\n'.encode('ascii'))
    assert not op.returncode
    # Clean up test directories
    smpd = os.path.join(os.path.expanduser('~'), 'ShakeMap')
    shutil.rmtree(os.path.join(smpd, 'test_profile1'))
    shutil.rmtree(os.path.join(smpd, 'test_profile2'))
    shutil.rmtree(os.path.join(smpd, 'test_profile3'))
    shutil.rmtree(os.path.join(smpd, 'test'))
    # Create a profile with directories we can't make (should fail)
    op = subprocess.Popen([program, '-c', 'junk_profile'], 
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=False)
    ipath = os.path.join(os.path.abspath(os.sep), 'zzxyzzskx') + '\n'
    op.communicate((ipath + ipath + ipath).encode('ascii'))
    assert op.returncode

    op = subprocess.Popen([program, '-c', 'junk_profile'], 
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=False)
    ipath = os.path.join(shakedir, 'junkprofile', 'install') + '\n'
    dpath = os.path.join(os.path.abspath(os.sep), 'zzxyzzskx') + '\n'
    op.communicate((ipath + dpath + dpath + dpath).encode('ascii'))
    assert op.returncode

    #
    # Create a profile for testing other programs
    # We leave this in place for later tests
    #
    op = subprocess.Popen([program, '-c', 'test_profile'], 
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=False)
    ipath = os.path.join(shakedir, 'tests', 'data') + '\n'
    dpath = os.path.join(shakedir, 'tests', 'data', 'eventdata') + '\n'
    op.communicate((ipath + dpath).encode('ascii'))
    assert not op.returncode

    
########################################################################
# Test sm_assemble
########################################################################
def test_assemble():
    installpath, datapath = get_config_paths()
    program = os.path.join(shakedir, 'sm_assemble')

    # Process a non-existent event (should fail)
    cp = subprocess.run([program, 'not_an_event'], shell=False)
    assert cp.returncode
    # This would succeed, but we remove event.xml (should fail)
    event_file = os.path.join(datapath, 'wenchuan', 'current', 'event.xml')
    os.rename(event_file, event_file + '_safe')
    cp = subprocess.run([program, 'wenchuan'], shell=False)
    assert cp.returncode
    os.rename(event_file + '_safe', event_file)
    # Normal event (should succeed)
    cp = subprocess.run([program, 'wenchuan'], shell=False)
    assert not cp.returncode
    # Do an event with model.conf (not model_zc.conf) and no zoneinfo
    # (should succeed)
    cp = subprocess.run([program, 'nc72282711'], shell=False)
    assert not cp.returncode
    # Try not having an event-specific config (should succeed)
    model_file = os.path.join(datapath, 'nc72282711', 'current',
                              'model.conf')
    os.rename(model_file, model_file + '_safe')
    cp = subprocess.run([program, 'nc72282711'], shell=False)
    assert not cp.returncode
    os.rename(model_file + '_safe', model_file)
    # Remove the existing hdf file from a no-backup even (should succeed)
    hdf_file = os.path.join(datapath, 'nc72282711_dyfi', 'current',
                            'shake_data.hdf')
    if os.path.isfile(hdf_file):
        os.rename(hdf_file, hdf_file + '_safe')
    cp = subprocess.run([program, 'nc72282711_dyfi'], shell=False)
    assert not cp.returncode
    if os.path.isfile(hdf_file + '_safe'):
        os.rename(hdf_file + '_safe', hdf_file)
    #
    # Try some bad config files
    #
    # Should fail validation
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad0', model_file)
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert cp.returncode
    os.rename(model_file + '_safe', model_file)
    # Should fail vs30 filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad1', model_file)
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert cp.returncode
    os.rename(model_file + '_safe', model_file)
    # Should fail prediction locations filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad2', model_file)
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert cp.returncode
    os.rename(model_file + '_safe', model_file)
    #
    # Make sure the location file substitutions work (should succeed)
    #
    cp = subprocess.run([program, 'northridge_points'], shell=False)
    assert not cp.returncode

########################################################################
# Test sm_augment
########################################################################
def test_augment():
    installpath, datapath = get_config_paths()
    assemble_program = os.path.join(shakedir, 'sm_assemble')
    program = os.path.join(shakedir, 'sm_augment')

    # Process a non-existent event (should fail)
    cp = subprocess.run([program, 'not_an_event'], shell=False)
    assert cp.returncode
    # This would succeed, but we remove shake_data.hdf (should fail)
    hdf_file = os.path.join(datapath, 'wenchuan', 'current', 'shake_data.hdf')
    os.rename(hdf_file, hdf_file + '_safe')
    cp = subprocess.run([program, 'wenchuan'], shell=False)
    assert cp.returncode
    os.rename(hdf_file + '_safe', hdf_file)
    # Normal event (should succeed)
    cp = subprocess.run([program, 'wenchuan'], shell=False)
    assert not cp.returncode
    # Do an event with model.conf (not model_zc.conf) and no zoneinfo
    # (should succeed)
    cp = subprocess.run([assemble_program, 'nc72282711'], shell=False)
    assert not cp.returncode
    cp = subprocess.run([program, 'nc72282711'], shell=False)
    assert not cp.returncode
    #
    # Make sure the location file substitutions work (should succeed)
    #
    cp = subprocess.run([program, 'northridge_points'], shell=False)
    assert not cp.returncode
    #
    # Try some bad config files
    #
    cp = subprocess.run([assemble_program, 'nc72282711_nodata_nofault'], 
                        shell=False)
    assert not cp.returncode
    # Should fail validation
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad0', model_file)
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert cp.returncode
    os.rename(model_file + '_safe', model_file)
    # Should fail vs30 filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad1', model_file)
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert cp.returncode
    os.rename(model_file + '_safe', model_file)
    # Should fail prediction locations filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad2', model_file)
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert cp.returncode
    os.rename(model_file + '_safe', model_file)
    #
    # Switch originators (should succeed)
    #
    model_file = os.path.join(datapath, 'nc72282711', 
                              'current', 'model.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.cz', model_file)
    cp = subprocess.run([program, 'nc72282711'], shell=False)
    assert not cp.returncode
    os.rename(model_file + '_safe', model_file)
    

########################################################################
# Test sm_model
########################################################################
def test_model():
    installpath, datapath = get_config_paths()
    #
    # This is Northridge for a set of output points (not a grid)
    # Run in verbose mode to hit that code
    # Remove the products directory to hit that code
    #
    program = os.path.join(shakedir, 'sm_assemble')
    cp = subprocess.run([program, 'northridge_points'], shell=False)
    assert not cp.returncode

    products_dir = os.path.join(datapath, 'northridge_points', 'current', 
                            'products')
    if os.path.isdir(products_dir):
        shutil.rmtree(products_dir)
    program = os.path.join(shakedir, 'sm_model')
    cp = subprocess.run([program, 'northridge_points', '-v'], shell=False)
    assert not cp.returncode
    #
    # This is a small grid with station data only
    #
    program = os.path.join(shakedir, 'sm_assemble')
    cp = subprocess.run([program, 'nc72282711'], shell=False)
    assert not cp.returncode

    program = os.path.join(shakedir, 'sm_model')
    cp = subprocess.run([program, 'nc72282711', '-v'], shell=False)
    assert not cp.returncode
    #
    # This is a small grid with DYFI data only
    #
    program = os.path.join(shakedir, 'sm_assemble')
    cp = subprocess.run([program, 'nc72282711_dyfi'], shell=False)
    assert not cp.returncode

    program = os.path.join(shakedir, 'sm_model')
    cp = subprocess.run([program, 'nc72282711_dyfi', '-v'], shell=False)
    assert not cp.returncode
    #
    # Run with no data and no fault, and use the default extent.
    #
    program = os.path.join(shakedir, 'sm_assemble')
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert not cp.returncode

    program = os.path.join(shakedir, 'sm_model')
    cp = subprocess.run([program, 'nc72282711_nodata_nofault'], shell=False)
    assert not cp.returncode
    #
    # Set the bias and outlier magnitude limits low to test additional
    # code branches
    #
    program = os.path.join(shakedir, 'sm_assemble')
    cp = subprocess.run([program, 'nc72282711_nofault'], shell=False)
    assert not cp.returncode

    program = os.path.join(shakedir, 'sm_model')
    cp = subprocess.run([program, 'nc72282711_nofault', '-v'], shell=False)
    assert not cp.returncode
    #
    # This event exists, but we hide the hdf file (should fail)
    #
    hdf_file = os.path.join(datapath, 'nc72282711_dyfi', 'current', 
                            'shake_data.hdf')
    os.rename(hdf_file, hdf_file + '_safe')
    program = os.path.join(shakedir, 'sm_model')
    cp = subprocess.run([program, 'nc72282711_dyfi'], shell=False)
    assert cp.returncode
    os.rename(hdf_file + '_safe', hdf_file)
    #
    # This event doesn't exist (should fail)
    #
    program = os.path.join(shakedir, 'sm_model')
    cp = subprocess.run([program, 'not_an_event'], shell=False)
    assert cp.returncode


########################################################################
# Put the user's profile back the way it was when we started
########################################################################
def test_restore():
    smrc = os.path.join(os.path.expanduser('~'), '.shakemap')
    if os.path.isdir(smrc):
        shutil.rmtree(smrc)
    if os.path.isdir(smrc + '_safe'):
        os.rename(smrc + '_safe', smrc)
    
if __name__ == '__main__':
    test_profile()
    test_assemble()
    test_augment()
    test_model()
    test_restore()
