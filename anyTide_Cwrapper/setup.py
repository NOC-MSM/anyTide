from distutils.core import setup, Extension
setup(name='polpredict', version='1.0',  \
      ext_modules=[Extension('polpredict', ['polpredict.cpp'])])
