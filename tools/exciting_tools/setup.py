from setuptools import setup
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='excitingtools',
      version='1.1.0',
      description='Utilities for aiding in the construction of exciting inputs and the postprocessing exciting '
                  'outputs.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='SOL Group',
      url="http://exciting.wikidot.com",
      author_email='peschelf@physik.hu-berlin.de',  # Add a point of contact
      include_package_data=True,
      python_requires=">=3.6",
      install_requires=['wheel>=0.35.0', 'numpy>=1.14.5', 'matplotlib>=2.2.0'],
      extras_require={"ase": ['ase>=3.20.0']}
      )
