# exciting Application Test Suite

This is the regression framework for the exciting application test suite. It depends on a number of python scripts that
automate running and evaluating system tests for the all-electron full-potential package exciting.

The test suite compares output files from an exciting calculation to results from a reference calculation, with certain 
tolerances.

Contact: [SOL group](https://sol.physik.hu-berlin.de/index.php?page=contact) at Humboldt-Universität zu Berlin. 

License: See README of exciting.

### Requirements:  
The test suite requires python 3.5 and the following packages:
 * argparse
 * collections
 * numpy
 * xml
 * unittest
 * excitingtools (packaged with exciting)

To install excitingtools, from exciting's root directory type:

    pip3 install -e tools/exciting_tools

### Using the Test Suite:

    python3 runtest.py (-a <action> -t <tests> -e <execuatble> -np <NP> -omp <OMP>)
        -h   --help          optional    show a help message and exit

        -a   --action        optional    Defines what action is done. That can be:
                                         'run'   - running tests 
                                         'ref'   - running references
                                                   References are always running with exciting_smp.
                                                   WARNING: This should only be done for NEW test cases or 
                                                   if all tests SUCCEEDED for the current version of the code.
                                         'clean' - cleaning test directories
                                         default is 'run'

        -t   --tests         optional    Test cases for the action. This can be either:
                                         <test case 1>, <test case 2>, ... - special test cases. This are the names of 
                                                                             the test case directories. Paths to the 
                                                                             test case directories work as well.
                                         'all'                             - all test cases
                                         default is 'all'

        -e   --executable    optional    exciting execuatbles.
	     		     		             'exciting_serial'  - for the serial binary;
					                     'exciting_smp'     - for the binary with SMP (default);
                        	             'exciting_mpi' - for the binary with MPI parallisation, only
					                     'exciting_mpismp'  - for the binary with MPI and SMP parallisation
					                     Default is 'exciting_smp'

        -np  --NP            optional    Number of cores for MPI parallelisation. Can only be used
                                         in combination with excitingmpi as executable.

        -omp --ompthreads    optional    Number of threads for open MP parallelisation. 
                                         Check envarOMP in runtest.py to make sure that the script 
                                         sets the right environment variable.
            				             Default is 2

        -handle-errors       optional    Allow assertion failures and skips to propagate to the
                                         end of the test suite. If the option is excluded, the
                                         default is to not allow error propagation

        -run-failing-tests   optional    Run tests tagged as failing.If the option is excluded,
                                         the default is not to run failing tests
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


### Adding new Tests:

Use the script `newtestcase.py` to add a new test case.

    python3 newtestcase.py -d <test description> -r <reference dir> -t <target dir> -i <init_template.xml>

        -h, --help        show a help message and exit

	    -d --description  Description of the new test case.

	    -r --reference    Path to an existing calculation that will be the reference.

	    -t --target       Target location to copy the reference to in the test farm.

	    -i --init_file    Init file for the new test case. The available init files can be found at '/xml/init_templates'. 
                          Default is init_groundstate.xml.

For example:

    python3 newtestcase.py -r test_reference -t test_farm/groundstate/new_test -i init_groundstate.xml

will copy (only) the input files from `test_reference` to the directory `test_farm/groundstate/new_test`, and subsequently
execute the input with the appropriate exciting binary to generate reference output data. Reference calculations are 
always running with exciting_smp.

An init.xml file defines the test tolerances. If `-i` is not specified, the script assumes a ground state calculation.
The default init.xml files can be found in `test/xml/init_templates`, but we provide a list of all possible pre-built init 
files here:

     init_groundstate.xml   - INFO.OUT, info.xml, atoms.xml, eigval.xml, evalcore.xml, geometry.xml
     init_dos.xml 	        - dos.xml
     init_bandstructure.xml - bandstructure.xml
     init_properties.xml    - RHO3D.xml, VCL3D.xml, VXC3D.xml, WF3D.xml, ELF3D.xml, EF3D.xml, LSJ.xml, EFG.xml, 
                              mossbauer.xml, expiqr.xml, effmass.xml, bandstructure.xml, dos.xml, KERR.OUT, 
                              EPSILON_11.OUT, EPSILON_12.OUT, EPSILON_33.OUT, CHI_111.OUT, ELNES.OUT
     init_GW.xml	        - EFERMI_GW.OUT, EVALQP.DAT, VXCNN.DAT, EPS00_GW.OUT
     init_BSE.xml	        - EPSILON/EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT, 
                              EPSILON/EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT, 
                              EPSILON/EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT, 
                              EXCITON/EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT, 
                              EXCITON/EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT, 
                              EPSILON/EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT, 
                              EXCITON/EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT,
                              EXCITON/EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT,
                              EXCITON/EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT
     init_wannier.xml	    - TDOS_WANNIER.OUT, WANNIER_INFO.OUT, POLARIZATION.OUT
     init_spintexture.xml   - spintext.xml

where the output files included in `init_groundstate.xml` are also included in all other init files. If the test needs 
special files and/or tolerances, it has to be added to the init.xml by hand. A brief description follows:

The root element is called <init> and has two attributes: 

| Attribute    | Description   |
|--------------|---------------|
|name          | Should be the same as directory name of the test case and is set automatically, if the script is used.  |
|description   | Includes a description of the calculation and the test. If nothing special is done, use the default explanation. |


For each exciting output file that is tested, there is a <test> element which has a number of attributes and elements. 
The most important attributes are:

| Attribute    | Type          |   Description   |
|--------------|---------------|-----------------|
|file          |  string       | Name of the exciting output file that is tested (e.g. atoms.xml). |
|tolFloat      | float/double  | Tolerance for float (double) values.|
|tolMSE        | float/double  | Tolerance for mean squred errors (e.g. energy bands). |
|condition     | two ints      | Allows for defining a number of values that are allowed to fail. First number for float failure, second for MSE failures.|
|ignore        | string        | Allows for defining values that will be ignored by the test (makes sense for machine dependent values like timings). The name of the values are the same as in the tested file. The names are seperated by a semicolon ";".|
