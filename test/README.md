# exciting Application Test Suite

This is the regression-testing framework for the exciting application test suite. It depends on a number of python scripts that
automate running and evaluating system tests for the all-electron full-potential package exciting.

The test suite compares output files from an exciting calculation to results from a reference calculation, with certain 
tolerances.

Contact: [SOL group](https://sol.physik.hu-berlin.de/index.php?page=contact) at Humboldt-Universität zu Berlin, or
post queries on the [MatSci exciting board](https://matsci.org/c/exciting/47).

License: See README of exciting.

## Requirements:  
The test suite requires python 3.5 (or above) and the following packages:
 * argparse
 * collections
 * numpy
 * json
 * excitingtools (packaged with exciting)

and additionally, pytest, if one wishes to run the framework's unit tests (defined in `test_functions/`).

To install excitingtools, from exciting's root directory type:

    pip3 install -e tools/exciting_tools

## Quickstart

To run the test suite with the default exciting binary, in this directory type:

```bash
python3 runtest.py
```

## Running the Test Suite:

    python3 runtest.py (-a <action> -t <tests> -e <executable> -np <NP> -omp <OMP>)
        -h   --help          optional    show a help message and exit

        -a   --action        optional    Defines what action is done. That can be:
                                         'run'   - running tests 
                                         'ref'   - running references
                                                   References are always running with exciting_smp.
                                                   WARNING: This should only be done for NEW test cases or 
                                                   if all tests SUCCEEDED for the current version of the code.
                                         Default is 'run'.

        -t   --tests         optional    Test(s) to run. 
                                         This can be some partial string match,
                                         such as `PBE`, or a full path `groundstate/LDA_PW-
                                         PbTiO3`. The test suite will run any test which
                                         contain the substring in its name. 
                                         Omit this option to run all tests.

        -e   --executable    optional    exciting executable.
                                         'exciting_serial'  - for the serial binary;
                                         'exciting_smp'     - for the binary with SMP (default);
                        	             'exciting_mpi'     - for the binary with MPI parallisation, only
                                         'exciting_mpismp'  - for the binary with MPI and SMP parallisation
                                         Default is 'exciting_smp'.

        -np  --NP            optional    Number of cores for MPI run. Can only be used in
                                         combination with exciting_mpi or exciting_mpismp as
                                         executable. Default is 2 for MPI and MPI+OMP
                                         calculations, and 1 for serial or pure OMP.

        -omp --ompthreads    optional    Number of threads for open MP parallelisation. Default is 2. 

        -handle-errors       optional    Allow testing to continue if a test assertion fails. 

        -run-failing-tests   optional    Run tests tagged as failing in src/exciting_settings/failing_tests.py
                                         Failing tests are not run by default. 

        -make-test           optional    Run tests from Makefile. If this option is set, all
                                         other options will be ignored and the test suite will
                                         run all tests with default settings. 
                                         The executable will be chosen from the compiled binaries
                                         with the following hierarchy: 
                                         1. exciting_mpismp 
                                         2. exciting_smp
                                         3. exciting_mpi
                                         4. exciting_serial
                                         If excitingtools is not installed, the test suite will provide 
                                         instructions on how to install the package.


## Adding new Tests:

Use the script `newtestcase.py` to add a new test case.

    python3 newtestcase.py -d <test description> -r <reference dir> -t <target dir> -i <init_template.xml>

        -h, --help      Show a help message and exit.

        -for            Method of the test case, informing the choice of tolerance
                        file. If not specified, the framework will attempt to infer 
                        the method from the target directory.

        -r --reference  Path to an existing calculation that will be the reference.

        -t --target     Target location to copy the reference to.


For example:

```bash
python3 newtestcase.py -for groundstate -r test_reference -t test_farm/groundstate/new_test
```

will copy input files from `test_reference` to the directory `test_farm/groundstate/new_test/ref`, execute the input 
with the appropriate exciting binary to generate reference output data and will add the file 
`tolerance_ground_state.json` to the directory.  Reference calculations are always run with `exciting_smp`.
If `-for` is not specified, the script will try to infer the calculation type from the target path and generate the 
corresponding tolerance file. 

The `tolerance_*.json` file defines the test tolerances, and has the following structure:

```json
{"files_under_test":
    ["INFO.OUT", "evalcore.xml", "geometry.xml", "eigval.xml", "atoms.xml"],
  "output_file1": {
    "property1": {
      "tol": 1e-08,
      "unit": "unit1"
     }
   },
  "output_file2": {
    "property2": {
      "tol": 1e-06,
      "unit": "unit2"
    },
    "property3": {
      "tol": 1e-08,
      "unit": "unit3"
    },
    "property4": {
      "tol": 1e-08,
      "unit": "unit4"
    }
  }
}
```

Every output file for the calculation type is listed as a dictionary. These dictionaries contain all properties parsed 
from the respective output file (but not necessarily every quantity in the output file - this depends on how the file
parser is written). For every property a tolerance and unit is defined. The units are currently not compared, however 
the plan is to include this in the future. See issue #100 in exciting's Gitlab. 

Furthermore, developers are free to modify (within reason) tolerances written to JSON, following generation of the file. 
The tolerance files are generated according to default values, defined in `src/tolerance/templates/`. If a tolerance 
value is consistently not appropriate, a new property is added to an output file, or indeed a new method is added to 
exciting, the developer should modify the appropriate template. 

Tolerance files also specify which files should be regression-tested. In simple cases like the ground state, each 
tested file has a dictionary of tolerances. However, for a method like BSE, which has numerous output files, it often
does not make sense to specify tolerances for each variation of one file. 

In this instance, the default `files_under_test` will be output containing wildcards (as specified in 
`src/utilities/wildcard_processor`):

```json
 {
  "files_under_test": [
    "EPSILON_??.OUT",
    "EXCITON_??.OUT"
  ]
}
```

This tells the test framework to use the same tolerances for all files that match `EPSILON_??.OUT`, and for all for that
match `EXCITON_??.OUT`, respectively. If one file, for example `EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT.xml`,
requires a different set of tolerance values. This can be defined explicitly:

```json
 {"files_under_test": [
    "EPSILON_??.OUT",
    "EXCITON_??.OUT",
    "EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT.xml"    
  ],
  "EPSILON_??.OUT": epsilon_tolerances_dict,
  "EXCITON_??.OUT": epsilon_tolerances_dict,
  "EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT.xml": specific_epsilon_tol_dict
}
```

Files under test will be removed from the `tolerance_*.json` files in the near-future, and placed in a general 
configuration file. The documentation will be updated when this change is made. Finally, tolerance files are 
automatically generated when running `newtestcase.py`, however, one is also free to generate them independently:

```bash
python3 src/reference/generate_tolerance.py -for <method>  
```

which is appropriate if one manually prepares the inputs and outputs for a test case. 


### Tolerance Templates

The default tolerances for each output file are defined in `src/tolerance/templates`, as python dictionaries:

```python3
# Tolerances for output.file
output_file = {'property1': default.energy,
               'property2': default.float,
               'property3': default.string}
```

where the dictionary values are entries in a named tuple, defined at the start of the template:

```python
default = DefaultTolerances(integer=Tol(0),
                            float=Tol(1.e-8),
                            str=Tol(''),
                            total_energy=Tol(1.e-8, Unit.hartree),
                            energy=Tol(1.e-7, Unit.hartree)
                            )
```

Here, the developer is free to define as many named entries as appropriate for a given method. If a unit is not known or 
required, one can choose `Unit.null`. Finally, the end of the template should contain a dictionary that contains
tolerances for all tested output files of a given method. For example, `ground_state_tolerances` looks like:

```python
# Single dictionary for all tested ground state outputs, for dumping to JSON
ground_state_tolerances = {'files_under_test': ['INFO.OUT', 'eigval.xml', 'evalcore.xml', 'atoms.xml', 'geometry.xml'],
                           'INFO.OUT': info_out,
                           'eigval.xml': eigval,
                           'evalcore.xml': evalcore,
                           'atoms.xml': atoms,
                           'geometry.xml': geometry
                           }
```

This dictionary also defines the default output files to be regression-tested, per module. As stated above, the 
templates should only be modified if a new value is parsed from an existing output file, or a new file or method is 
added to exciting. 

## Parsed Data Structure Expected by the Test Framework

The test framework expects all data to be given to it (parsed) as python dictionaries. In that respect, using structured
output in files makes parsing trivial (JSON, XML, YAML). If one opts for non-structured output, then they should still 
ensure sufficient information is available in the file header, regarding the nature of the data structuring.

**If the developer writes a new output file, they are also responsible for writing its parser.**


### Modifying Existing Output File Keys

If one modifies an output file, they should be aware that in some cases (INFO.OUT, for example) they are also modifying
the keys of the associated parser. In general, no new parsers should auto-generate keys - they should be defined 
explicitly, such that the keys are decoupled from the formatting in the output file.

In all cases, if one changes a key in a parser, then they **must** change all associated tolerance files, plus the 
method template to be consistent. In general, the structure of parsed data should be considered carefully before 
committing to it.


## Guidelines for Writing a Parser

The tolerance comparison will only evaluate the values of lowest-nested keys (a deliberate design choice). As such, one 
should consider how they structure the parsed data. For example, it makes more sense to structure data like:

```JSON
{'wannier1': {'localisation_vector': np.array(shape=(3)),
              'Omega': float
             }
}
```

such that the tolerances can be assigned w.r.t. `localisation_vector`, and `Omega`, rather than using the structure:

```JSON
{'localisation_vector': {'wannier1':  np.array(shape=(3))
                         'wannier2':  np.array(shape=(3))
                        },
 'Omega': {'wannier1',:  float
           'wannier2':  float
          }
}
```

which will results in tolerances defined w.r.t. `wannier1` and `wannier2`. One can see in the latter case, there is 
no distinction between `localisation_vector` and `Omega`, one can only set a tolerance for each wannier function. In 
general, we want to set different tolerances for properties, rather than for different objects with the same set of 
properties. As a rule, it therefore makes sense to structure parsed data like:

```JSON
{'object1': {'attribute1': array,
             'attribute2': scalar
             },
'object2': {'attribute1': array,
            'attribute2': scalar
             }
}
```

One could also structure the data like:

```JSON
{ 'n_wannier': 5, 
  'localisation_vector': np.array(shape=(n_wannier, 3)),
  'Omega': : np.array(shape=(n_wannier))
}
```

where a reduction in serialisation avoids the problem.
