from setuptools import setup

setup(name='spectrassembler',
      version='0.0.2',
      description='Tool (experimental) to compute layout from overlaps with spectral algorithm',
      url='https://github.com/antrec/spectrassembler',
      author='Antoine Recanati',
      author_email = 'antoine.recanati@inria.fr',
      license='MIT',
      scripts = ['spectrassembler.py', 'ioandplots.py', 'overlaps.py', 'spectral.py', 'consensus.py'],
      install_requires=[
        'numpy',
        'scipy',
        'biopython']
      )
