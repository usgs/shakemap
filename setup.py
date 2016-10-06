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
                'shakemap.grind',
                'shakemap.grind.conversions',
                'shakemap.grind.conversions.imc',
                'shakemap.grind.conversions.imt',
                'shakemap.grind.correlation',
                'shakemap.grind.directivity',
                'shakemap.grind.gmice',
                'shakemap.genex',
                'shakemap.mapping',
                'shakemap.plotting',
                'shakemap.utils',
                'shakemap.utils.extern.scp',
                'shakemap.transfer'],
      package_data={'shakemap': [os.path.join('tests', 'data', '*'),
                                 os.path.join('utils', 'configspec.ini'),
                                 os.path.join('grind', 'data', 'ps2ff', '*.csv')]},
      scripts=['mkfault'],
      )
