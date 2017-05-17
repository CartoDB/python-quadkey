from distutils.core import setup, Extension

module1 = Extension('quadkey', sources = ['quadkey.c'])

setup (name = 'Quadkey',
        version = '1.0',
        description = 'quadkey',
        ext_modules = [module1])
