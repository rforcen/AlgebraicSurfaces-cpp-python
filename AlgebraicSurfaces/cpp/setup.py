# https://docs.python.org/3/distutils/apiref.html
from distutils.core import setup
from distutils.extension import Extension

AlgebraicSurfaces = Extension(
    'AlgebraicSurfaces',
    sources=['AlgebraicSurfaces.cpp'],
    libraries=['boost_python37-mt', 'boost_numpy37-mt'],
    extra_compile_args=['-std=c++17']  # lambda support required
)

setup(
    name='AlgebraicSurfaces',
    version='0.1',
    ext_modules=[AlgebraicSurfaces])

# call with: python3.7 setup.py build_ext --inplace
