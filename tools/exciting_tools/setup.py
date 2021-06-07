from setuptools import setup

setup(name='excitingtools',
      version='0.0.1',
      description='Utilities for aiding in the construction of exciting inputs and the postprocessing exciting outputs.',
      author='SOL Group',
      author_email='abuccheri@physik.hu-berlin.de',  # Add a point of contact
      packages=['excitingtools'],
      install_requires=[
          'numpy>=1.14.5',
          'matplotlib>=2.2.0']
      )
