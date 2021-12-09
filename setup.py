from setuptools import setup
from cmake_setuptools import *

setup(
    name='healpeak',
    version='',
    url='',
    license='',
    author='Pier Fiedorowicz',
    author_email='pierfied@email.arizona.edu',
    description='Fast peak/void counting for HEALPix maps.',
    ext_modules=[CMakeExtension('healpeak')],
    cmdclass={'build_ext': CMakeBuildExt}
)
