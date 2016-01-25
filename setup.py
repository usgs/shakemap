from distutils.core import setup

setup(name='shakemap',
      version='0.1dev',
      description='Python GMPE library',
      author='Bruce Worden, Mike Hearne',
      author_email='cbworden@usgs.gov,mhearne@usgs.gov',
      url='',
      package_data = {'shakemap':['shakelib/configspec.ini']},
      packages=['shakemap','shakemap.gmpe','shakemap.shakelib','shakemap.oldshakelib','shakemap.correlation','shakemap.gmice','shakemap.directivity'],
      scripts = ['mkfault.py'],
)
