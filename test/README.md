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
`tolerance_ground_state.json` to the directory.  Reference calculations are always run with `exciting_serial`.
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

In simple cases like the ground state, each tested file has a dictionary of tolerances. However, for a method like BSE, 
which has numerous output files, it often does not make sense to specify tolerances for each variation of one file. 

In this instance, file name keys can be specified with wildcards (as defined in `src/utilities/wildcard_processor`). The 
example below tells the test framework to use the same tolerances for all files that match `EPSILON_??.OUT`, and for all 
for that match `EXCITON_??.OUT`, respectively. Furthermore, if one file,
`EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT.xml` for example, requires a different set of tolerance values it can 
be defined explicitly:

```json
{
  "EPSILON_??.OUT": epsilon_tolerances_dict,
  "EXCITON_??.OUT": epsilon_tolerances_dict,
  "EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT.xml": specific_epsilon_tol_dict
}
```

Finally, tolerance files are automatically generated when running `newtestcase.py`, however, one is also free to 
generate them independently:

```bash
python3 src/reference/generate_tolerance.py -for <method>  
```

which is appropriate if one manually prepares the inputs and outputs for a test case. 

## Manually Preparing Test Cases

The developer is also free to manually prepare test cases. Indeed, `newtestcase.py` will only copy `input.xml` and 
any species files present in the initial directory, so more elaborate tests must be done manually. To manually prepare 
a test case:

```bash
# Prepare inputs and tolerance in test_farm
python3 src/reference/generate_tolerance.py -for <method>
mkdir test_farm/method/test_case/ref
cp <INPUT FILES> test_farm/method/test_case/ref/*
mv tolerance_method.py test_farm/method/test_case/ref/.

# Run test
cd test_farm/method/test_case/
export OMP_NUM_THREADS= 2
./$EXCITINGROOT/bin/exciting_smp

# Remove outputs that will not be tested
rm <IRRELEVANT_OUPUTS> 
```

where `$EXCITINGROOT` is exciting's root directory, `<method>` corresponds to the specific exciting method,
`<INPUT FILES>` corresponds to all inputs required for the calculation, and `<IRRELEVANT_OUPUTS>` are all files not
under test.

### Rerunning Reference Data For a Test Case

The test framework does not provide a means of rerunning an existing test case to refresh reference data. Instead, this
is done straightforwardly from the terminal:

```bash
cd test_farm/method/test_case/ref
export OMP_NUM_THREADS= 2
./$EXCITINGROOT/bin/exciting_smp
```

The command:

```bash
python3 runtest.py -a ref 
```

should be **AVOIDED**, and is only intended for use when a bug is found that affects all test cases (in the ground state, 
for example). This feature has been flagged to be removed from the `runtest.py` script. See issue 116 in exciting's
Gitlab issue tracker. 

### Test Configuration File

Once reference and tolerance data are generated for a test, it should be added to the configuration file:

```bash
test/tests_config.yml
```

This file defines all tests that are part of the test suite (and should be consistent with those present in `test_farm/`).
Each test has several *properties*, which can either be explicitly defined, or omitted (causing the default to be used). 
Valid properties are:

```yaml
group:            Developer-defined test group
repeat:           True or False
files_under_test: Each output file to be compared to reference data
inputs:           All input files required for the calculation
failing_builds:   Any build types for which the test fails
comments:         Information on failing builds
```

At the time of writing, most tests have been specified with all of their properties omitted. This is **not** encouraged 
for new test cases. Be **explicit** with inputs and files under test. A working test case will look something like:

```yaml
method/test_name:
   files_under_test:
      - "some_output.OUT"
   inputs:
      - "input.xml"
      - "some_species_file.xml"
```

Test cases which fail for some builds will contain additional information:

```yaml
# YAML can include unparsed comments with the hash character
# TODO(Alex) Issue 101
groundstate/LDA_PW-collinear-Fe:
   group: NONE
   repeat: False
   files_under_test:
      - "INFO.OUT"
      - "evalcore.xml"
      - "geometry.xml"
      - "eigval.xml"
      - "atoms.xml"
   inputs:
      - "input.xml"
      - "Fe.xml"
   failing_builds:
      - intel_mpiandsmp
      - intel_serial
      - gcc_mpiandsmp
      - gcc_serial
   comments: 'Most energies differ to reference by ~1.e-7. \n
   scl%Sum of eigenvalues by ~ 1.e-6. \n
   DOS at Fermi differs by 5.7e-04. '
```

### Test Case Properties

Addressing each property:

#### group

`group` allows developer-defined grouping of tests. It could be that some tests are slow, or are only meant to be run on 
specific hardware. Grouping allows more flexibility with what the CI runs. All groups, and whether or not they run, is 
defined in `defaults_config.yml`. If the `group` property is not specified in a test case, `NONE` is assigned:

```yaml
group: NONE 
```

If a new group is added, the developer is also required to add the correspond entry to the enum class `TestGroup`, in
`src/runner/configure_tests.py`. 

#### repeat

`repeat` specifies if a test should be repeated if it fails. This is intended for flakely tests that **sometimes** fail, 
however in general, one should avoid committing flakely tests. The default is `False` (although it would make more sense 
to make it a number, and remove the `-repeat-tests N` command line arugment):

```yaml
repeat: False 
```

Note, tests will only be repeated if `repeat: True` is specified in the config file and `-repeat-tests N` is given as 
a command line argument to the test suite, where `N` is an integer `> 0`.

#### files_under_test

`files_under_test` specifies which exciting outputs should be regression-tested for each test case. Any output that has
a corresponding parser can be included. All defaults are specified in `test/defaults_config.yml`. The defaults are
used by omitting this property. Else, one specifies files like so:

```yaml
files_under_test:
   - "INFO.OUT"
   - "evalcore.xml"
   - "geometry.xml"
   - "eigval.xml"
   - "atoms.xml"
```

#### inputs

`inputs` defines the input files required to run a calculation (specifically, which files to copy to the `run/` directory
of each test case). This should include species files. Explictly specifying inputs is particularly useful if a test 
a) restarts or b) uses non-standard inputs. `inputs` are specified in the same way as `files_under_test`.

Default inputs are not defined in `test/defaults_config.yml`, rather the test framework will inspect the test `ref/` 
directory for `input.xml` and any species files, if the property is not specified. 

#### failing_builds and comments

In some cases, tests will fail with a given build stack. `failing_builds` can be used to specify if a test case fails 
only for a given build. This is required due to the test migration from the 2020 suite, and where possible should not be 
used for new tests. Valid failing builds choices are:

```yaml
failing_builds:
  - intel_serial
  - intel_smp
  - intel_purempi
  - intel_mpiandsmp
  - gcc_serial
  - gcc_smp
  - gcc_purempi
  - gcc_mpiandsmp
comments: "This test case fails for all builds - see issue i."
```

Additionally, `comments` can be added to give more information on what fails. This can be printed by the framework
at the end of execution. 


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


## Conclusion

This concludes the explanation of exciting's regression-testing framework. Happy testing. 
