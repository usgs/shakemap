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
      packages=['shakemap','shakemap.gmpe','shakemap.grind','shakemap.genex',
                'shakemap.mapping','shakemap.utils','shakemap.correlation',
                'shakemap.gmice','shakemap.directivity','shakemap.transfer',
                'shakemap.extern.scp'],
      package_data = {'shakemap':[os.path.join('utils', 'configspec.ini')],
                      'shakemap':[os.path.join('grind', 'data', 'ps2ff', '*.csv')]},
      scripts = ['mkfault', 'mkinputdir', 'mkscenariogrids'],
)
