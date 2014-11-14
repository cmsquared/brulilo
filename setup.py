from setuptools import setup

setup(name='brulilo',
      version='0.1',
      description='Toy nuclear reaction network scaffolding.',
      url='http://bitbucket.org/ChrisMalone/brulilo',
      author='Chris Malone',
      author_email='chris.m.malone@gmail.com',
      packages=['brulilo'],
      install_requires=[
          'numpy',
          'matplotlib',
          'lxml',
          'scipy',
      ],
      license='TBD')
