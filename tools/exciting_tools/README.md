# exciting Tools
*exciting tools* is intended to be a collection of scripts to facilitate the generation of exciting inputs and 
the post-processing of exciting outputs. 

## Installation
If one wishes to import *exciting tools* in their own scripts, it can be installed from this project's root directory 
(`exciting/tools/exciting_tools`) with:

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
The intention is to a) utilise NOMAD parser for exciting and b) move to structured output in exciting. In the meantime, 
it is reasonable to add parsers for exciting outputs as required, as long as these routines **only** perform parsing
and return dictionaries. This is to make the eventual switch to NOMAD parser much easier.
