import os
from distutils.core import setup
import os.path
import versioneer
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import shutil

# This should be handled by conda when we install a platform-specific
# compiler, but apparently isn't on macs (yet?)
if shutil.which('clang') is None:
    os.environ['CC'] = 'gcc'
else:
    os.environ['CC'] = 'clang'

sourcefiles = ["shakemap/c/pcontour.pyx", "shakemap/c/contour.c"]

clib_source = ["shakemap/c/clib.pyx"]

ext_modules = [Extension("shakemap.c.pcontour",
                         sourcefiles,
                         libraries=["m"],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=[]),
               Extension("shakemap.c.clib",
                         clib_source,
                         libraries=['m', 'omp'],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=["-fopenmp"],
                         extra_link_args=["-fopenmp"]
                         )]

cmdclass = versioneer.get_cmdclass()
cmdclass['build_ext'] = build_ext

setup(name='shakemap',
      version=versioneer.get_version(),
      description='USGS Near-Real-Time Ground Motion Mapping',
      author='Bruce Worden, Mike Hearne, Eric Thompson',
      author_email='cbworden@usgs.gov,mhearne@usgs.gov,emthompson@usgs.gov',
      url='http://github.com/usgs/shakemap',
      packages=[
          'shakemap',
          'shakemap.utils',
          'shakelib',
          'shakelib.conversions',
          'shakelib.conversions.imc',
          'shakelib.conversions.imt',
          'shakelib.correlation',
          'shakelib.directivity',
          'shakelib.gmice',
          'shakelib.rupture',
          'shakelib.plotting',
          'shakelib.utils',
      ],
      package_data={
          'shakemap': [
              os.path.join('tests', 'shakemap', 'data', '*'),
              os.path.join('data', '*'),
          ],
          'shakelib': [
              os.path.join('test', 'data', '*'),
              os.path.join('utils', 'data', '*'),
              os.path.join('rupture', 'ps2ff', '*.csv')
          ]
      },
      scripts=[
          'bin/associate_amps',
          'bin/pgm2json',
          'bin/getdyfi',
          'bin/receive_amps',
          'bin/receive_origins',
          'bin/receive_origins_gsm',
          'bin/run_verification',
          'bin/shake',
          'bin/sm_check',
          'bin/sm_compare',
          'bin/sm_create',
          'bin/sm_migrate',
          'bin/sm_profile',
          'bin/sm_queue',
          'bin/sm_rupture',
          'bin/sm_batch',
          'bin/sm_sync'],
      cmdclass=cmdclass,
      ext_modules=cythonize(ext_modules)
      )
