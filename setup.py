from distutils.core import setup
import os.path
import versioneer


setup(name='shakemap',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='USGS Near-Real-Time Ground Motion Mapping',
      author='Bruce Worden, Mike Hearne, Eric Thompson',
      author_email='cbworden@usgs.gov,mhearne@usgs.gov,emthompson@usgs.gov',
      url='http://github.com/usgs/shakemap',
      packages=['shakemap',
                'shakemap.mapping',
                'shakemap.utils',
                'shakemap.utils.extern.scp'],
      package_data={'shakemap': [os.path.join('tests', 'data', '*'),
                                 os.path.join('utils', 'configspec.ini'),
                                 ]},
      scripts=['sm_assemble', 'sm_augment', 'sm_model', 'sm_profile'],
      )
