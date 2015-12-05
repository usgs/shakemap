from distutils.core import setup

setup(name='pyshake',
      version='0.1dev',
      description='Python GMPE library',
      author='Bruce Worden, Mike Hearne',
      author_email='cbworden@usgs.gov,mhearne@usgs.gov',
      url='',
      package_data = {'pyshake':['shakelib/configspec.ini']},
      packages=['pyshake','pyshake.gmpe','pyshake.shakelib','pyshake.oldshakelib','pyshake.correlation','pyshake.gmice'],
      scripts = ['mkfault.py'],
)
