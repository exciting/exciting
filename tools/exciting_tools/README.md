# exciting Tools
*exciting tools* is intended to be a collection of scripts to facilitate the generation of exciting inputs and 
the post-processing of exciting outputs. 

## Installation
If one wishes to import *exciting tools* in their own scripts, it can be installed from this project's root directory 
(`<EXCITINGROOT>/tools/exciting_tools`) with:

```angular2html
pip install -e .
```
## External Package Dependencies
If a new external dependency is introduced to the package, this also requires adding to `setup.py` such that pip is aware 
of the new dependency.

## Basic File Structure 
In general, modules should begin with a docstring giving an overview of the module's purpose. External python
libraries should then be imported, followed by a space, then local modules belonging to *exciting tools*. Local modules 
should be loaded with relative paths, rather than prepending the system path `sys.path.insert(0,'/path/to/module_directory')`:

```angular2html
"""
Functions that operate on lattice vectors 
"""
import numpy as np

from .maths.math_utils import triple_product
```
This may change in the future in favour of loading modules in `__init__.py`. 

## Code Formatting 
We're currently favouring [yapf](https://github.com/google/yapf) formatter, which by default applies PEP8 formatting to 
the code, however even formatter that applies the PEP8 standard is sufficient. 

## Documentation 

### Writing Documentation
All functions and classes should be documented. The favoured docstring is *reStructuredText*:

```angular2html
class SimpleEquation:
   def demo(self, a: int, b: int, c: int) -> list:
    """
    Function definition

    :param int a: quadratic coefficient
    :param int b: linear coefficient 
    :param c: free term
    :type c: int
    :return list y: Function values   
    """
```
where the type can be specified in the `param` description, or separately using the `type` tag. For more details on the
documentation syntax, please refer to this [link](https://devguide.python.org/documenting/). The [google style guide](
https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html) for *reStructuredText* docstrings is also 
acceptable to follow. 

### Generating Documentation 
- [ ] TODO(Alex) Issue 57 Set up generation of documentation from docstrings, and description here. 

## Testing 

Every function should have a test where possible, unless the function is correct by inspection. The naming convention 
for a module called `module.py` is to prepend the it with `test_`:
```angular2html
excitingtools/module.py       # Collection of functions
tests/test_module.py          # Collection of tests for functions in module.py
```
Tests are intended to be run using *pytest*, for which the documentation can be found [here](https://docs.pytest.org/en/stable/index.html). 
One is able to run `pytest` from the `exciting_tools` root with no arguments. By default, all test files, classes and functions defined in the specification,
`exciting_tools/pytest.ini`,  will get executed. 


## Parsers 

The parsers are used in the test suite. Therefore, they should only return dictionaries with a specific structure.
 
The tolerance comparison will only evaluate the values of lowest-nested keys. As such, one should consider how they structure the parsed data. 
For example, it makes more sense to structure data like:
```python3
{‘wannier1’: {‘localisation_vector’: np.array(shape=(3)),
              ‘Omega’: float
             }
}
```
such that the tolerances will be w.r.t. `localisation_vector`, and `Omega`, rather than using the structure:
```python3
{‘localisation_vector’: {‘wannier1’:  np.array(shape=(3))
                         ‘wannier2’:  np.array(shape=(3))
                        },
 ‘Omega’: {‘wannier1’:  float
           ‘wannier2’:  float
          }
}
```
which will results in tolerances defined w.r.t. `wannier1` and `wannier2`. One can see in the latter case, there is no distinction between `localisation_vector` and `Omega`. In general, we’re more likely to want to set different tolerances for different properties, rather than for different functions with the same set of properties.
One could also structure the data like:
```python3
{‘localisation_vector’: np.array(shape=(n_wannier, 3)),
 ‘Omega’: : np.array(shape=(n_wannier)
}
```
where the less serialised data removes the key nesting.

## Uploading to PyPi

excitingtools is available as a separate package on PyPi. In order to upload a new version:

```bash
# Ensure twine is installed
pip3 install twine
# Build the wheels
cd $EXCITINGROOT/tools/exciting_tools
python3 setup.py sdist bdist_wheel

# Test the distribution and uploading (one requires a test-PyPi account)
twine check dist/*
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

# Upload to PyPi
twine upload dist/*
```

Before doing so, please ensure the semantic versioning is appropriately updated in `setup.py`.


## Contributors
The following people (in alphabetic order by their family names) have contributed to excitingtools:

* Alexander Buccheri
* Hannah Kleine
* Martin Kuban
* Benedikt Maurer
* Fabian Peschel
* Daniel Speckhard
* Elisa Stephan
* Mara Voiculescu
