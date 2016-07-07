from distutils.core import setup
import os.path

setup(name='shakemap',
      version='0.1dev',
      description='USGS Near-Real-Time Ground Motion Mapping',
      author='Bruce Worden, Mike Hearne',
      author_email='cbworden@usgs.gov,mhearne@usgs.gov',
      url='',
      packages=['shakemap','shakemap.gmpe','shakemap.grind','shakemap.genex',
                'shakemap.mapping','shakemap.utils','shakemap.correlation',
                'shakemap.gmice','shakemap.directivity'],
      package_data = {'shakemap':[os.path.join('utils', 'configspec.ini')],
                      'shakemap':[os.path.join('grind', 'data', 'ps2ff', '*.csv')]},
      scripts = ['mkfault', 'mkinputdir'],
)
