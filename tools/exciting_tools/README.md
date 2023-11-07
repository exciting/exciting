# excitingtools
<span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span> is a collection of 
modules to facilitate the generation of <span style="font-family:american typewriter; font-size:1em;">**exciting**</span>
inputs and the post-processing of <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> outputs. 

<span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span> currently provides functionality for:

* Generation of the <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> input XML file 
  using Python classes:
  - Automatically supported for the whole input file through dynamic class construction
  - Currently tested for `groundstate`, `structure`, `BSE` and `bandstructure`


* Parsing of <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> outputs into Python dictionaries


* High-level class API for interacting with results:
  - Currently implemented for eigenvalues, band structure and DOS (without SO coupling)

making it is possible to define a calculation, run it, and parse the relevant outputs all from within Python. 

<span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span> is used by, or in conjunction with:
* <span style="font-family:american typewriter; font-size:1em;">**exciting's**</span> regression-testing framework
  * Parsing of output data
* <span style="font-family:american typewriter; font-size:1em;">**exciting's**</span> Jupyter notebook tutorials 
  * Data handling
* [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)
  * Input and output handling in ASE's <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> calculator
* [Jobflow](https://github.com/materialsproject/jobflow)
  * For the development of complex, automated <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> workflows  

## Installation
If one wishes to import <span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span> in their own scripts, it can be installed from this project's root directory 
(`$EXCITING_ROOT/tools/exciting_tools`).

Although not strictly necessary, it is strongly recommended to use a virtual environment to manage Python packages. 
There are several solutions available, including [venv](https://docs.python.org/3/library/venv.html), (mini)[conda](https://docs.conda.io/en/latest/miniconda.html)
and [pipx](https://pypa.github.io/pipx/). 

To set up a venv:

```bash
# Create a directory `excitingvenv` containing our venv
python3 -m venv excitingvenv
# Activate the environment
source excitingvenv/bin/activate
```

Before installing excitingtools, you should upgrade `pip` and `setuptools`:
```bash
python3 -m pip install --upgrade --force pip
pip install --upgrade setuptools
```

No matter which type of environment or none is used, it should always be verified that the Python version is compatible with <span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span> (>=3.7) and that `pip` is up-to-date (>=22.0) and also points to the same Python version.

Now you can proceed to install <span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span>
from <span style="font-family:american typewriter; font-size:1em;">**exciting's**</span>'s root:

```bash
python3 -m pip install tools/exciting_tools
```

Alternatively, you can download it directly from PyPI:

```bash
pip install excitingtools
```

## External Package Dependencies
If a new external dependency is introduced to the package, this also requires adding to `pyproject.toml` such that pip is aware 
of the new dependency.

## Basic File Structure 
In general, modules should begin with a docstring giving an overview of the module's purpose. External python
libraries should then be imported, followed by a space, then local modules belonging to <span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span>. Local modules 
should be loaded with absolute paths rather than relative paths or prepending the system path `sys.path.insert(0,'/path/to/module_directory')`:

```angular2html
"""
Functions that operate on lattice vectors 
"""
import numpy as np

from excitingtools.maths.math_utils import triple_product
```
Exposed modules, forming user API, should be defined in `__init__.py` where ever possible.

## Code Formatting 
We currently favour [yapf](https://github.com/google/yapf) formatter, which by default applies PEP8 formatting to 
the code.  

After installing yapf, if you are in the root directory of excitingtools, you can simply type:
```bash
yapf -i excitingtools/path/to/file.py
```
and it will do the formatting for you. Note: This will automatically use our custom `.style.yapf` style-file.

## Documentation 

### Writing Documentation
All functions and classes should be documented. The favoured docstring is *reStructuredText*:

```python3
class SimpleEquation:
   def demo(self, a: int, b: int, c: int) -> list:
    """Function definition.

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

Documentation can straightforwardly be generated using the [pdoc](https://docs.python.org/3/library/pydoc.html) package:

```bash
pip install pdoc
pdoc -o documentation -d restructuredtext --math excitingtools/
```

- [ ] TODO(Alex) Issue 57. Set up generation of documentation from docstrings, with Sphinx 

### Basic Usage

#### Input XML Generation

<span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span> maps the XML tags and attributes 
of `input.xml` onto Python classes, enabling the generation of XML-formatted inputs directly from Python. A simple 
ground state calculation could like this:

```python3
import ase
import numpy as np

from excitingtools import ExcitingStructure, ExcitingGroundStateInput, ExcitingInputXML

# Lattice and positions in angstrom, as expected by ASE
lattice = np.array([[3.168394160510246,   0.0,                0.0],
                    [-1.5841970805453853, 2.7439098312114987, 0.0],
                    [0.0,                 0.0,                39.58711265]])
positions = np.array([[0.00000000, 0.00000000, 16.68421565],
                      [1.58419708, 0.91463661, 18.25982194],
                      [1.58419708, 0.91463661, 15.10652203],
                      [1.58419708, 0.91463661, 22.90251866],
                      [0.00000000, 0.00000000, 24.46831689],
                      [0.00000000, 0.00000000, 21.33906353]])
symbols = ['W', 'S', 'S', 'Mo', 'S', 'S']
atoms = ase.atoms.Atoms(symbols=symbols, positions=positions, cell=lattice)

structure = ExcitingStructure(atoms, species_path='.')

ground_state = ExcitingGroundStateInput(
    rgkmax=8.0,
    do="fromscratch",
    ngridk=[6, 6, 6],
    xctype="GGA_PBE_SOL",
    vkloff=[0, 0, 0],
    tforce=True,
    nosource=False
    )

input_xml = ExcitingInputXML(structure=structure, 
                             groundstate=ground_state, 
                             title="My exciting Crystal")

input_xml.write("input.xml")
```
Here we defined the attributes required to perform a ground state calculation as seperate classes, and composed the 
final XML string with `ExcitingInputXML` class. If the user does not have access to ASE, they can instead use a 
`List[dict]` to define the container with atoms data:

```python3
atoms = [{'species': 'W', 'position': [0.00000000, 0.00000000, 16.68421565]},
         {'species': 'S', 'position': [1.58419708, 0.91463661, 18.25982194]},
         {'species': 'S', 'position': [1.58419708, 0.91463661, 15.10652203]},
         {'species': 'Mo','position': [1.58419708, 0.91463661, 22.90251866]},
         {'species': 'S', 'position': [0.00000000, 0.00000000, 24.46831689]},
         {'species': 'S', 'position': [0.00000000, 0.00000000, 21.33906353]}]

structure = ExcitingStructure(atoms, lattice, species_path='.')
```

Additional examples can be found in the test cases, `exciting_tools/tests/input`. We note that not all XML tags 
currently map onto Python classes. One can consult `exciting_tools/excitingtools/input` to see what is available. 
Development follows a continuous integration and deployment workflow, therefore if one wishes for additional features, 
please make a request on GitHub issues or open a merge request.

#### Binary Execution

Next we can define a runner and run our calculation:
```python3
from excitingtools.runner.runner import BinaryRunner

runner = BinaryRunner('exciting_smp', run_cmd='', omp_num_threads=4, time_out=500)
run_status = runner.run()
```

#### Parsing Outputs

After the successful completion of the calculation, we can parse the relevant output files as dictionaries, using
`parser_chooser`. These are the main files one would be interested in after performing a ground state calculation, for
example:

```python3
from excitingtools import parser_chooser

info_out: dict = parser_chooser("INFO.OUT")
eigval_info: dict = parser_chooser("eigval.xml")
atoms_info: dict = parser_chooser("atoms.xml")
```

A full list of parsers is provided in `excitingtools/exciting_dict_parsers/parser_factory.py`. If we wish to perform
analysis of the data, <span style="font-family:american typewriter; font-size:1em;">**excitingtools**</span> provides
classes with high-level API. To perform a band structure plot using the `BandData` class:

```python
""" Plot silicon band structure
"""
import matplotlib.pyplot as plt

from excitingtools.exciting_obj_parsers.ks_band_structure import parse_band_structure
from excitingtools.dataclasses.band_structure import BandData

band_data: BandData = parse_band_structure("bandstructure.xml")
vertices, labels = band_data.band_path()

ha_to_ev = 27.2114
fig, ax = plt.subplots(figsize=(6, 9))

ax.set_xticks(vertices)
ax.set_xticklabels(labels)
plt.ylabel('Energy (eV)')

# Font sizes
ax.yaxis.label.set_size(20)
ax.tick_params(axis='both', which='major', labelsize=20)

# Vertical lines at high symmetry points
for x in vertices:
    plt.axvline(x, linestyle='--', color='black')

# Fermi reference
e_fermi = 0.0
# Number of valence bands
n_valence = 4

# Colour valence and conduction bands differently
line_colour = {key:'blue' for key in range(0, n_valence)}
line_colour.update({key:'red' for key in range(n_valence, band_data.n_bands)})

for ib in range(0, band_data.n_bands):
    plt.plot(band_data.flattened_k_points, ha_to_ev * band_data.bands[:, ib], color=line_colour[ib])
```

Tests demonstrating further usage are present in `excitingtools/tests/dataclasses`. We note that the high-level objects
and their parsers are separated. In principle, the data classes should only define a sensible schema or API for 
accepting relevant data, rather than know anything about the parsing. Object parsers (defined in `obj_parsers`) by definition
should return to data classes, but the data classes dictate the format of the data, not vice versa. 

## Testing 

Every function should have a test where possible, unless the function is correct by inspection. The naming convention 
for a module called `module.py` is to prepend it with `test_`, which allows it to be automatically recognised and run
by *pytest*:

```bash
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
pip3 install --upgrade twine build
# Build the wheels
cd $EXCITINGROOT/tools/exciting_tools
python3 -m build

# Test the distribution and uploading (one requires a test-PyPi account)
twine check dist/*
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

# Upload to PyPi
twine upload dist/*
```

Before doing so, please ensure the semantic versioning is appropriately updated in `pyproject.toml`.


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
