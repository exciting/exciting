from setuptools import setup, find_packages

setup(name='exciting_tools',
      version='0.0.1',
      description='Utilities for aiding in the construction of exciting inputs and the postprocessing exciting outputs.',
      author='SOL Group',
      author_email='abuccheri@physik.hu-berlin.de',  # Add a point of contact
      packages=find_packages(include=['exciting_tools', 'exciting_tools.*']),
      install_requires=[
          'numpy>=1.14.5',
          'matplotlib>=2.2.0']
      )
