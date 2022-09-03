#!/usr/bin/env python3

from distutils.core import setup

def get_version_from_Cargo_toml():
    marker = 'version ='
    for line in open('Cargo.toml'):
        if line.startswith(marker):
            return eval(line[len(marker):])


script_names = '''viewraw
                  xenon_thickness_from_h5
'''.split()

setup(name        = 'rustpetalo-python-utils',
      version     = get_version_from_Cargo_toml(),
      scripts     = [f'src/{name}.py' for name in script_names],
      py_modules  = ['utils'],
      package_dir = {'': 'src'},
      )
