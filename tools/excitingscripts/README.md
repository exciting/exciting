# excitingscripts

`excitingscripts` is a collection of various Python scripts for executing various tasks with the [exciting](https://exciting-code.org) code from the command line.

## Installation

`excitingscripts` can be installed directly from PyPI:

```bash
pip install excitingscripts
```

Alternatively, it can be installed from the [exciting source][https://github.com/exciting/exciting]:

```bash
cd $EXCITINGROOT/tools/excitingscripts
pip install -e .
```

## Usage 

You can execute each script as a Python module, e.g., to run a single `exciting` calculation, use 

```bash
python3 -m excitingscripts.execute.single
```

Below, we list all scripts, and the functions that are used in each script.

NOTE: This README is automatically generated. To update the listing, change the docstrings on top of each script file, or function, respectively.

## Script Listing

### excitingscripts.checkfit


#### excitingscripts.checkfit.checkfit

Check-fit implementation.

##### arg_parser

Get the arg parser for checkfit scripts.  
  
:param quantity: what derivatives to check fit  
:return: the argparser

##### fit

Fit data to a polynomial.  
  
:param order: order of the fit  
:param x: x data  
:param y: y data  
:param order_of_derivative: determines which coefficient is taken from the fit  
:return: the chosen coefficient, if the fit is poorly conditioned "None"

##### get_unit_conversion_factor

Get the unit conversion factor.  
  
:param lattice_param_exp: exponent for the lattice parameter in the unit conversion.  
:return: unit conversion factor

##### print_info_to_stdout

Print some information to the terminal.  
  
:param max_displacement: maximum chosen displacement for the fit  
:param frequencies: the computed frequencies  
:param order_of_derivative: fit order of interest  
:param n_max: number of chosen displacement values

##### quantity_specific_checkfit

Specific checkfit implementation for either energy or force.  
  
:param quantity: name of the quantity key in the result file  
:param factor: factor for the unit conversion  
:param lattice_param_exp: exponential for the lattice parameter in the unit conversion  
return: checkfit function


#### excitingscripts.checkfit.energy_vs_displacement

Python script for extracting the derivatives at zero displacement of the energy-vs-displacement curves.

##### arg_parser

Get the arg parser for checkfit scripts.  
  
:param quantity: what derivatives to check fit  
:return: the argparser

##### quantity_specific_checkfit

Specific checkfit implementation for either energy or force.  
  
:param quantity: name of the quantity key in the result file  
:param factor: factor for the unit conversion  
:param lattice_param_exp: exponential for the lattice parameter in the unit conversion  
return: checkfit function


#### excitingscripts.checkfit.energy_vs_strain

Python script for extracting derivatives at zero strain of energy-vs-strain curves.

##### fit

Fit data to a polynomial.  
  
:param order: order of the fit  
:param x: x data  
:param y: y data  
:param order_of_derivative: determines which coefficient is taken from the fit  
:return: the chosen coefficient, if the fit is poorly conditioned "None"

##### parse_energy_vs_strain

Read the "energy_vs_strain" file  
  
:param directory: directory containing the file  
:param max_strain: value of maximum strain  
:return: strain and energy values

##### parse_info_elastic_constants

Gives the necessary info to be printed  
  
:param directory: path to the directory containing the file "INFO-elastic-constants"  
:return : dictionary containing the information

##### print_info_to_stdout

Print some information to the terminal.  
  
:param info:  
:param max_strain: maximum chosen strain for the fit  
:param derivatives: the computed derivatives  
:param order_of_derivative: fit order of interest  
:param n_max: number of chosen strain values

##### save_derivatives_to_json

Saving the data in fitted derivative data in a JSON format  
  
:param directory: path to save the file  
:param order_of_derivative: order of derivative for the fit  
:param results: results of the fits


#### excitingscripts.checkfit.energy_vs_volume

Python script for fitting energy-vs-volume curves

##### fit_energy_vs_volume

Fit the energy-vs-volume data using a polynomial fit.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:param order_of_fit: Order of the polynomial fit.  
:param isym: lattice symmetry code.  
  
:return:Polynomial fit, bulk modulus, minima, lattice constant, chi value

##### plot_energy_vs_volume

Plot the energy-vs-volume data and save the plot.  
  
:param volumes : List of volumes.  
:param energies : List of energies.  
:param curv : Polynomial fit curve.  
:param dmin : Minima of the fit.  
:param order_of_polynomial : Order of the polynomial fit.  
:param output_dir : Directory to save the plot.

##### print_info_to_stdout

Print information onto the screen  
  
:param isym: lattice symmetry code.  
:param dmin: minima of the fit.  
:param lattice_const: lattice constant for the given lattice symmetry code  
:param bulk_modulus: bulk modulus  
:param chi: chi-squared value indicating the goodness of fit


#### excitingscripts.checkfit.force_vs_displacement

Python script for extracting the derivatives at zero displacement of the force-vs-displacement curves.

##### arg_parser

Get the arg parser for checkfit scripts.  
  
:param quantity: what derivatives to check fit  
:return: the argparser

##### quantity_specific_checkfit

Specific checkfit implementation for either energy or force.  
  
:param quantity: name of the quantity key in the result file  
:param factor: factor for the unit conversion  
:param lattice_param_exp: exponential for the lattice parameter in the unit conversion  
return: checkfit function


#### excitingscripts.compare_transition_energies

Determine the transition energies for the transitions Γ→Γ and Γ→X for given directories in which exciting  
calculations have been performed.  
  
Located at `excitingscripts/compare_transition_energies.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.compare_transition_energies -r dir1 dir2 dir3  
```  
Where <code>dir1</code>, <code>dir2</code>, <code>dir3</code> take the place of the names of the  directories where exciting calculations have been performed, which need to be specified in order to calculate the transition energies. The script can be used for any number of directories.

##### determine_transition_energies

Determine the transition energies for the transitions Γ→Γ and Γ→X for a given directory in which an exciting  
calculation was performed.  
  
:param root_directory: Root directory.


### excitingscripts.convert


#### excitingscripts.convert.au2invcm

##### convert_q_points_from_atomic_units_to_inverse_cm

Display frequencies of the phonon modes from the exciting output file PHONON.OUT in inverse centimeters instead  
of atomic units.  
  
:param run_dir: Root directory.


#### excitingscripts.convert_xml2xsf

Convert **xml** files to **xsf**.  
  
Located at `excitingscripts/convert_xml2xsf.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.convert_xml2xsf -f file -d dimension  
```  
Where <code>file</code> is the **xml** file to be converted to **xsf** and <code>dimension</code> is the dimension of <code><span style="color:green">plot</span></code> sub-element in the <code><span style="color:green">properties</span></code> element for a given exciting calculation.

##### convert_xml2xsf

Convert a given XML file to XSF.  
  
:param file_to_convert: XML File to convert.  
:param dimension: Dimension of "plot" sub-element in "properties" element for the given exciting calculation.  
:param excitingroot: Environment variable string.


### excitingscripts.execute


#### excitingscripts.execute.convergence_test

Run a series of **exciting** calculations with different values of the main computational parameters.  
  
Located at `excitingscripts/execute/convergence_test.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.execute.convergence_test k_i k_f rgkmax_i rgkmax_f  
```  
Where <code>k_i</code> and <code>k_f</code> are the initial and final k-values for defining the <code><span style="color:green">groundstate</span></code> attribute <code><span style="color:MediumBlue">ngridk</span></code>, and <code>rgkmax_i</code> and <code>rgkmax_f</code> the initial and final values for the <code><span style="color:green">groundstate</span></code> attribute <code><span style="color:MediumBlue">rgkmax</span></code>.

##### execute_convergence_test

Execute a series of exciting calculations with varying values for the groundstate attributes ngridk and rgkmax  
 and return a list containing the total energy value for each set of parameters.  
  
:param k_initial: Initial k-value for defining the groundstate attribute ngridk.  
:param k_final: Final k-value for defining the groundstate attribute ngridk.  
:param rgkmax_initial: Initial value for the groundstate attribute rgkmax.  
:param rgkmax_final: Final value for the groundstate attribute rgkmax.  
:param root_directory: Root directory.  
:param excitingroot: Environment variable string.  
:returns: List containing total energy values for each set of parameters.

##### run_exciting

Execute an exciting calculation in a given running directory.  
  
:param root_directory: Root directory.  
:param excitingroot: Environment variable string.  
:param filename: Name of the exciting input file  
:param timeout: Maximum runtime in seconds


#### excitingscripts.execute.diamond_phonon

Execute phonon diamond calculations.

##### execute_diamond_phonon

Executes a series of exciting diamond calculations to get phonons.  
  
:param work_dir: Working directory containing the input files  
:param excitingroot: root directory of exciting  
:return: phonon results from the calculations

##### run_exciting

Execute an exciting calculation in a given running directory.  
  
:param root_directory: Root directory.  
:param excitingroot: Environment variable string.  
:param filename: Name of the exciting input file  
:param timeout: Maximum runtime in seconds


#### excitingscripts.execute.elastic_strain

Run a series of **exciting** calculations with different strain values.  
  
Located at `excitingscripts/execute/elastic_strain.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.execute.elastic_strain   
```

##### execute_elastic_strain

Execute a series of exciting calculations with different interlayer distances.  
  
:param root_directory: Root directory.  
:param dft_half: Boolean with "True" value for DFT-1/2 calculations.  
:param excitingroot: Environment variable string.  
:returns: Array with energy-strain data.

##### run_exciting

Execute an exciting calculation in a given running directory.  
  
:param root_directory: Root directory.  
:param excitingroot: Environment variable string.  
:param filename: Name of the exciting input file  
:param timeout: Maximum runtime in seconds


#### excitingscripts.execute.monitored


#### excitingscripts.execute.planar_average

Extract planar-averaged electrostatic potential in a given direction.  
  
Located at `excitingscripts/execute/planar_average.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.execute.planar_average direction  
```  
Where <code>direction</code> is the direction along which the plane-averaged potential will be visualized.

##### execute_planar_average

Extract planar-averaged electrostatic potential in a given direction.  
Script included in "execute" directory for consistency, due to classification of the old "tutorial scripts"  
  
:param potential_file: File containing electrostatic potential data.  
:param direction: Direction along which potential needs to be averaged.


#### excitingscripts.execute.single

Run a single **exciting** calculation.  
  
Located at `excitingscripts/execute/single.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.execute.single -r rundir  
```  
Where <code>rundir</code> is an optional parameter which specifies the running directory. If <code>rundir</code> is not specified, the calculation will run in the directory where the script is called.

##### run_exciting

Execute an exciting calculation in a given running directory.  
  
:param root_directory: Root directory.  
:param excitingroot: Environment variable string.  
:param filename: Name of the exciting input file  
:param timeout: Maximum runtime in seconds


#### excitingscripts.execute.volume_optimization

Run a series of **exciting** calculations for structures with different volumes.  
  
Located at `excitingscripts/execute/volume_optimization.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.execute.volume_optimization nr_vol  
```  
Where <code>nr_vol</code> is the number of volume values for which structures are generated by varying the lattice constant.

##### execute_volume_optimization

Execute a series of exciting calculations with different volumes obtained by varying the lattice constant.  
  
:param number_volume_values: Number of volume values for which structures are generated by varying the lattice constant.  
:param root_directory: Root directory.  
:param excitingroot: Environment variable string.  
:returns: NumPy array containing total energy values and corresponding volume values.

##### run_exciting

Execute an exciting calculation in a given running directory.  
  
:param root_directory: Root directory.  
:param excitingroot: Environment variable string.  
:param filename: Name of the exciting input file  
:param timeout: Maximum runtime in seconds


#### excitingscripts.free_from_dos

##### heat_capacity_distribution

Calculate distribution which is later multiplied by the phonon DOS to determine the heat capacity.  
  
:param omega: Frequency.  
:returns Distribution which is needed for determining the heat capacity.

##### normalized_frequency

Calculate normalized frequency for a given temperature and frequency.  
  
:param temp: Temperature.  
:param omega: Frequency.  
:returns Distribution which is needed for determining the heat capacity.

##### thermodynamic_properties

Calculate thermodynamic properties at a given temperature.  
  
:param temp: Temperature.  
:param omega_data: Frequency values.  
:param dos_data: DOS values.  
:param energy_unit: Energy unit.  
:returns: Vibrational free energy, vibrational internal energy, entropic contribution to the vibrational free  
energy, vibrational entropy, heat capacity, zero-point energy

##### vibrational_expression

Calculate distribution which is later multiplied by the phonon DOS to determine the vibrational free energy.  
  
:param omega: Frequency.  
:returns Distribution which is needed for determining the vibrational free energy.


### excitingscripts.lattice


#### excitingscripts.lattice.parameters

##### convert_exciting_to_sgroup

Convert an exciting input file  into an sgroup input file.  
  
:param input_file: Path to the exciting input.xml file.  
:param output_file: Path to the sgroup input file.

##### get_parameters

Extract lattice symmetry and parameters from the sgroup.out file.  
  
:param directory: Directory containing the sgroup.out file.  
:return: Dictionary with lattice symmetry and parameters.

##### print_parameters

Print the space group parameters in a formatted way.  
  
:param parameters: Dictionary with lattice symmetry and parameters.

##### run_sgroup

Run the sgroup program to extract the space group information from the input file.  
  
:param input_file: Path to the exciting input.xml file.  
:param output_dir: Path to the output directory.


### excitingscripts.optimize


#### excitingscripts.optimize.analyze

Python script for fitting the energy-vs-volume and energy-vs-strain curves.

##### fit_energy_data

Fit the energy data and write the optimized structure to a new .xml file.  
  
:param run_dir: Path to the directory containing the calculations.  
:param vol_flag: Flag to indicate if the calculation is a volume optimization.  
:param dir_name: Name of the directory indicating the type of optimization.  
:param fit_type: Type of fit.

##### fit_energy_vs_strain

Fit the energy data to a polynomial and return the optimized structure.  
  
:param strain: List of strains.  
:param energies: List of energies.  
:param out_file: Path to the output file.  
:param order_of_fit: Order of the polynomial fit.  
:return: Fitted parameters

##### fit_energy_vs_volume

Fit the energy data to an equation of state and return the optimized structure.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:param out_file: Path to the output file.  
:param eos: Equation of state.  
:return: Fitted parameters.

##### get_deformation_matrix

Get deformation matrix for a given strain value depending on the optimization parameter.  
  
:param eps: Strain value.  
:return: Dictionary containing deformation_matrix for each optimization parameter.

##### get_energy_data

Extracts and saves the energy-volume data from a series of exciting calculations.  
  
:param run_dir: Path to the directory containing the calculations.  
:param dir_name: Name of the directory indicating the type of optimization.  
:param vol_flag: Flag to indicate if the calculation is a volume optimization.

##### murnaghan_eos

Calculate energy using the Murnaghan equation of state.  
  
:param v: Volume.  
:param p: Parameters for the equation [v0, e0, b0, bp].  
:return: Calculated energy.

##### plot_energy

Plot the fitted energy data and save the plot  
  
:param energies: List of energies.  
:param x_data: List of volumes/strain.  
:param output_dir: Directory to save the plot.  
:param curv: Fitted curve.  
:param p: Fitted parameters.  
:param vol_flag: Flag to indicate if the calculation is a volume optimization.  
:param out_name:  
:param fit_type: Type of fit.

##### pressure_murnaghan_eos

Calculate pressure using the Murnaghan equation of state.  
  
:param v: Volume.  
:param p: Parameters for the equation [v0, e0, b0, bp].  
:return: Calculated pressure.

##### residuals

Calculate residuals for the least squares fit.  
  
:param p: Parameters for the Murnaghan EOS.  
:param e: Energies.  
:param v: Volumes.  
:param eos: Equation of state.  
:return: Residuals.


#### excitingscripts.optimize.setup

Python script for generating structures at different volume/strains.

##### check_monoclinic_compatibility

Checks if the given monoclinic structure is compatible with certain geometric criteria.  
  
:param base_vectors: The base vectors of the crystal structure.  
:param ref_scale: The reference scale for the unit cell.  
:param stretch: The stretch factors along each axis.  
:param threshold_angle: The threshold for determining if the angle is effectively 90 degrees.

##### get_crystal_system

Get crystal system based on the space group number.  
  
:param space_group_number: Space group number of the crystal structure.  
:return: Crystal system corresponding to the space group number.

##### get_deformation_matrix

Get deformation matrix for a given strain value depending on the optimization parameter.  
  
:param eps: Strain value.  
:return: Dictionary containing deformation_matrix for each optimization parameter.

##### setup_optimize_lattice

Set up the optimization of the lattice parameters.  
  
:param max_strain: The maximum physical strain.  
:param num_dist_str: The number of distorted structures.  
:param infile: Name of input file.  
:param opt_index: The index of the optimization parameter.


#### excitingscripts.optimize.submit

Python script for running a series of exciting calculations.

##### run_optimize_lattice

Run a series of exciting calculations with different volumes generated by the script  
'excitingscripts.optimize.setup'.  
  
:param run_dir: Directory where exciting runs.


### excitingscripts.plot


#### excitingscripts.plot.atomforce

This script allows for the visualization of the force-vs-displacement curve.

##### plot_atomforce

This script allows for the visualization of the force-vs-displacement curve.


#### excitingscripts.plot.band_structure

Plot single and multiple electronic and phonon band structures (BS)  
  
Require the following files:  
- for electronic BS: input.xml, BAND.OUT, BAND-QP.OUT, BAND_WANNIER.OUT, BANDLINES.OUT  
- for phonon BS: input.xml, PHDISP.OUT, PHLINES.OUT  
  
More details can be found **[here](https://www.exciting-code.org/home/the-python-script-plot.band_structure)**.  
  
Located at `excitingscripts/plot/band_structure.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.band_structure  
```

##### check_number_of_plots

print warning if the number_of_plots is larger than 4:  
un update of the colors list could be necessary.  
no_leg = input_options['no_legend']

##### extract_title_text

extract text in title from xmlfile

##### extract_xticks_labels

extract x-ticks labels for band structure from xmlfile

##### find_steps

find number of k values in the electronic energies of a single band

##### inquire_element

check if element is in xmlfile

##### inquire_file

inquire file existence

##### inquire_spin

check in input.xml if a spin-polarized calculation is performed

##### option_parser

Parse command line inputs   
  
Parse:   
    directory  
    eboundary  
    assign_type  
    phonon  
    no_legend  
    ezero  
    scale_box  
    kboundary  
    kpoint_boundary  
    eunit  
    funit  
    title  
    no_title  
    legend_position  
    scale_box  
    invert_colors  
    invert_plots  
    max_ticks_y  
    legend_label  
    show  
  
:return input_options: Dictionary of parsed command line arguments

##### plot_ezero

plot the energy zero

##### read_electronic_band

read electronic band-structure frome files infile

##### read_phonon_dispersion

read phonon-dispersion curves frome files infile

##### read_xticks_position

read position of ticks on the horizontal axis

##### set_energy_zero

set energy zero, assume in input bands aligned to the Fermi energy

##### set_legend_label

set legend label for the plot


#### excitingscripts.plot.bbirch

Python script for fitting energy-vs-volume curves

##### bulk_modulus_finite_difference

Calculate bulk modulus using the Finite Difference method  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:return: Calculated bulk modulus.

##### findindex

Finds the index of given value in the list upto specified tolerance  
  
:param x: value for which index is found  
:param y: list in which the value is searched for  
:param dymax: tolerance for which value can be found  
:return : index of the value in the list

##### fit_pressure_vs_volume

Fit the energy-vs-volume data using the Birch-Murnaghan Equation of State.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
  
:return:Birch fit, bulk modulus, bulk_modulus_pressure_deriv, minima, lattice constant, chi value

##### plot_bulk_modulus_vs_volume

Plot the bulk modulus vs volume data and save the plot.  
  
:param volumes : List of volumes.  
:param volumes_fd : List of volumes at which bulk modulus is calculated using finite difference.  
:param bulk_modulus_fd: List of bulk modulus calculated using finite difference.  
:param curv : Fitted curve.  
:param dmin : Minima of the fit.  
:param output_dir : Directory to save the plot.

##### print_info_to_stdout

Print information onto the screen  
  
:param dmin: minima of the fit.  
:param lattice_const: lattice constant for the given lattice symmetry code  
:param bulk_modulus: bulk modulus  
:param bulk_modulus_pressure_deriv: derivative of bulk modulus with respect to pressure  
:param chi: chi-squared value indicating the goodness of fit


#### excitingscripts.plot.birch

##### fit_energy_vs_volume

Fit the energy-vs-volume or energy-vs-strain data using the Birch-Murnaghan Equation of State.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:return:Birch fit, bulk modulus, bp, minima, lattice constant, chi value

##### parse_energy_file

Reads the "energy-vs-volume" or "energy-vs-strain" file  
  
:param data_file: path to the file  
:return: strain/volume and energy values

##### plot_energy

Plot the energy-vs-volume or energy-vs-strain data and save the plot.  
  
:param X : List of volumes.  
:param energies: List of energies  
:param curv : Fitted curve.  
:param dmin : Minima volume of the fit.  
:param output_dir : Directory to save the plot.  
:param dim: Dimension of the system  
:param v_eq: Volume at zero strain

##### print_info_to_stdout

Print information onto the screen  
  
:param dmin: minima of the fit.  
:param lattice_const: lattice constant for the given lattice symmetry code  
:param bulk_modulus: bulk modulus  
:param bp: derivative of bulk modulus with respect to pressure  
:param chi: chi-squared value indicating the goodness of fit  
:param dim:  
:param v_eq:

##### strain_to_volume

:param strain: Strain Value  
:param dim: Dimension of the system  
:param v0: Volume at zero strain  
:return: Volume value

##### volume_to_strain

:param volume: Volume value  
:param dim: Dimension of the system  
:param v0: Volume at zero strain  
:return: Strain Value


#### excitingscripts.plot.bondlength

Python visualization tool for following relative atomic coordinates of atoms during the relaxation process.

##### calculate_bond_lengths

Calculate bond lengths between specified atoms from position data.  
  
:param position_data: list of dictionaries containing atomic positions.  
:param atom1: first atom number.  
:param atom2: second atom number.  
:param lattice_matrix: unit cell containing the molecule  
:param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.  
:param threshold: threshold value for periodic boundary conditions  
  
:return: list of calculated bond lengths.

##### plot_bond_lengths

Plot the bond lengths and save the plot to a file.  
  
:param bond_lengths: list of bond lengths to plot.  
:param run_dir: directory where exciting runs.  
:param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.  
:param atom1:  first atom number  
:param atom2: second atom number  
:param show: whether to display the plot.  
:param dpi: resolution in DPI for the saved plot.


#### excitingscripts.plot.centerofmass

Python visualization tool for the cartesian components of the position of the center of mass during relaxation

##### plot_center_of_mass

Plot the center of mass components, then save the plot to a file.  
  
:param center_of_mass_data: A list of tuples containing the center of mass data.  
:param run_dir: directory where exciting runs  
:param show: Whether to display the plot.  
:param dpi: Resolution in DPI for the saved plot.


#### excitingscripts.plot.checkderiv

Plot derivatives.  
  
This is a very important tool that allows to represent the dependence of the calculated derivatives of the  
energy-vs-displacement and force-vs-displacement curves on  
  
 * the range of points included in the fitting procedure ("maximum displacement u"),  
 * the maximum degree of the polynomial used in the fitting procedure ("n").

##### find_first_none

Finds the first none in the list.  
  
:param inp: input list  
:return: the index of the first none in the inp list

##### plot_checkderiv

Plot the derivative of the energy-vs-displacement or force-vs-displacement curves.  
  
:param quantity: of interest, could be 'energy' or 'force' or 'strain'  
:param y_min_arg: lower limit of the y-axis for plotting  
:param y_max_arg: upper limit of the y-axis for plotting


#### excitingscripts.plot.compare_vdW

Visualize multiple energy-vs-distance curves.  
  
Located at `excitingscripts/plot/compare_vdW.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.compare_vdW -f file_name -r dir1 dir2 dir3  
```

##### plot_compare_vdw

Plot binding energy curve values for a given running directory.  
  
:param plot_file_path: Path to file containing data wanted for plot.  
:param color_index: Index needing for plotting curves with different colors for each calculation.


#### excitingscripts.plot.convergence

Visualize convergence results.  
  
Located at `excitingscripts/plot/convergence.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.convergence plot_mode  
```  
Where <code>plot_mode</code> is either <code>k</code> for plotting energy curves with varying values of the <code><span style="color:green">groundstate</span></code> attribute <code><span style="color:MediumBlue">ngridk</span></code>, <code>r</code> for varying values of the <code><span style="color:green">groundstate</span></code> attribute  <code><span style="color:MediumBlue">rgkmax</span></code> or <code>rk</code> for a 3D plot with varying values of both attributes.

##### plot_convergence_k

Plot energy curves for varying values of the groundstate attribute ngridk.

##### plot_convergence_r

Plot energy curves for varying values of the groundstate attribute rgkmax.

##### plot_convergence_rk

Plot energy curves for varying values of the groundstate attributes ngridk and rgkmax.


#### excitingscripts.plot.dos

Visualize desisty of states.  
  
More details can be found **[here](https://www.exciting-code.org/home/the-python-script-plot.dos)**.  
  
Located at `excitingscripts/plot/dos.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.dos  
```

##### check_number_of_plots

print warning if the number_of_plots is larger than 4:  
un update of the colors list could be necessary.  
no_leg = input_options['no_legend']

##### extract_title_text

extract text in title from xmlfile

##### find_steps



##### inquire_element

check if element is in xmlfile

##### inquire_file

inquire file existence

##### inquire_spin

check in input.xml if a spin-polarized calculation is performed

##### option_parser

Parse command line inputs   
  
Parse:   
    directory  
    eboundary  
    assign_type  
    phonon  
    no_legend  
    scale_box  
    dos_boundary  
    eunit  
    funit  
    title  
    no_title  
    legend_position  
    scale_box  
    reverse_colors  
    reverse_plots  
    no_fill  
    no_reverse_spin  
    max_ticks_x  
    max_ticks_y  
    legend_label  
    grid  
    show  
       
:return input_options: Dictionary of parsed command line arguments

##### set_legend_label

set legend label for the plot


#### excitingscripts.plot.energy

Visualize energy-vs-strain curves.  
  
Located at `excitingscripts/plot/energy.py`.   
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.energy  
```

##### sortstrain

Sort strain values and also sort energy values based on the index of the sorted strain list.  
  
:param strain: List containing strain values.  
:param strain: List containing energy values.  
:returns: Lists containing sorted strain and energy values.


#### excitingscripts.plot.exciton_weights

Visualize energy-vs-strain curves.  
  
Located at `excitingscripts/plot/exciton_weights.py`.   
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.exciton_weights structure_name file_name energy_min energy_max exciton_weights_size  
```  
Where <code>structure_name</code> is the name of the structure, <code>file_name</code> is the name of the file containing data needed for exciton visualization, <code>energy_min</code> and <code>energy_max</code> are the minimum and maximum energy values for setting plot axis limits, and <code>exciton_weights_size</code> is the size of excitonic weights.

##### plot_exciton_weights

Plot excitonic weights along a band structure path.  
  
Assumes presence of "bandstructure.dat" file in current running directory for plotting band structure.  
  
:param structure_name: Name of structure.  
:param exciton_file: File containing data needed for exciton visualization.  
:param energy_min: Minimum energy value for setting plot axis limit.  
:param energy_max: Maximum energy value for setting plot axis limit.  
:param exciton_weights_size: Size of excitonic weights.


#### excitingscripts.plot.files

Visualize data included in different files and different directories.  
  
More details can be found **[here](https://www.exciting-code.org/home/the-python-script-plot.files)**.  
  
Located at `excitingscripts/plot/files.py`.   
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.files  
```

##### check_number_of_plots

inquire file existence

##### inquire_file

inquire file existence

##### option_parser

Parse command line inputs  
  
Parse:  
    directory  
    files  
    legend_label  
    column_x  
    column_y  
    xboundary  
    yboundary  
    label_x  
    label_y  
    x_scale  
    y_scale  
    legend_position  
    title  
    no_title  
    max_ticks_x  
    max_ticks_y  
    no_legend  
    grid  
    scale_box  
    reverse_colors  
    reverse_plots  
    no_scientific  
    log_x  
    log_y  
    show  
    plot_name  
  
:return input_options: Dictionary of parsed command line arguments

##### read_data

read data from files infile


#### excitingscripts.plot.maxforce

Python visualization tool for the maximum amplitude of the force on the atoms during relaxation.

##### plot_forces

Plot the torque components and magnitude, then save the plot to a file.  
  
:param forces: list of forces and target.  
:param run_dir: directory where exciting runs.  
:param show: whether to display the plot.  
:param dpi: resolution in DPI for the saved plot.


#### excitingscripts.plot.multitask

Please, check https://www.exciting-code.org/the-python-script-plot.multitask  
to better understand how to use this script

##### Handle_complex

enum to treat handle_complex

##### Option_preprocess

enum to store which action to take as preprocessing

##### arrange_as_stack

Function to convert list of lists to a "stack" (only a list)  
:param list_of_lists: list with lists

##### convert_args_to_dict

Parsed arguments are converted to a dictionary  
:param args: command line arguments

##### fft

FFT - from time domain to (angular) frequency:  
:param t: time (array)  
:param f: function (array) f(t) to be fourier-transformed  
:param wcut: cut-off frequency (in Ha) for the low pass filter  
:return w, F: Tuple containing angular frequencies (in Ha) and the  
    fourier transform

##### format_plot



##### parse_input

Function to parse the arguments from the command line

##### preformat_plot

Set up some pyplot parameters, pre-formatting the plot

##### preprocess

Preprocess x and y, depending on the desired options  
:param x: list of arrays, data to be outputed but needed for preprocessing y  
:param y: list of arrays, data to preprocess  
:param str output: name of output file  
:param option_preproc: required kind of preprocessing  
:param wcut: smoothing parameter for the fourier transform  
:param diagonal_component: to obtain the dielectric tensor, we need to  
    know if the desired component belongs to the diagonal

##### read_file

FFT time to (angular) frequency:  
:param str file: name of the file to be read  
:param columns_to_plot: array of 2 integers with the columns to be read from file  
:param nlines_skip: number of lines to skip  
:param scale: array of two real numbers to scale x and y  
:param handle_complex: enum with the treatment of complex numbers  
:return x, y: Tuple x and y arrays as read from file with minimal processing

##### sanity_checks

Function to make some sanity checks of the command line arguments  
:param options: dictionary with the command line arguments

##### set_implicit_options

#-------------------------------------------------------------------------------


#### excitingscripts.plot.newbirch

Fit energy-vs-volume curves using the Birch-Murnaghan equation of state (**BM-EoS**) in polynomial form.  
  
Located at `excitingscripts/plot/newbirch.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.newbirch  
```

##### sortstrain




#### excitingscripts.plot.optimized_geometry

Python visualization tool for relaxed coordinates of atoms in the unit cell.

##### get_optimized_relative_coordinates

Get the optimized relative coordinates between two atoms.  
  
:param atom1: first atom number.  
:param atom2: second atom number.  
:param root_directory: root directory.  
:return: List of relative coordinates.

##### plot_optimized_geometry

Plot the optimized geometry.  
  
:param max_strains: List of strain values.  
:param optimized_geometries: List of relative coordinates.  
:param ymin: Minimum value for the y-axis.  
:param ymax: Maximum value for the y-axis.  
:param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.


#### excitingscripts.plot.pbirch

Python script for fitting energy-vs-volume curves

##### fit_pressure_vs_volume

Fit the energy-vs-volume data using the Birch-Murnaghan Equation of State.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
  
:return:Birch fit, bulk modulus, bulk_modulus_pressure_deriv, minima, lattice constant, chi value

##### plot_pressure_vs_volume

Plot the pressure-vs-volume data and save the plot.  
  
:param volumes : List of volumes.  
:param volumes_fd : List of volumes at which pressure is calculated using finite difference.  
:param pressure_fd: List of pressures calculated using finite difference.  
:param curv : Fitted curve.  
:param dmin : Minima of the fit.  
:param output_dir : Directory to save the plot.

##### pressure_finite_difference

Calculate pressure using the Finite Difference method  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:return: Calculated pressure.

##### print_info_to_stdout

Print information onto the screen  
  
:param dmin: minima of the fit.  
:param lattice_const: lattice constant for the given lattice symmetry code  
:param bulk_modulus: bulk modulus  
:param bulk_modulus_pressure_deriv: derivative of bulk modulus with respect to pressure  
:param chi: chi-squared value indicating the goodness of fit


#### excitingscripts.plot.phonon_anim

##### generate_visualization_files

Generate .axsf and .xyz files for visualizing phonon modes.  
  
:param supercell_dims: Dimensions of the supercell as (n1, n2, n3).  
:param scaling: Scaling factor for the atomic displacements.  
:param nsteps: Number of steps in the animation sequence.  
:param root_dir: Directory containing the files input.xml and PHONON.OUT.

##### parse_species_data

Parse data for each species to extract mass, atomic number, and count information.  
  
:param unique_species: List of unique species in the structure.  
:param excitingroot: exciting root directory.  
:param all_species: List of all species present in the structure.  
:return: Species data containing chemical symbol, mass, atomic number and number of atoms per species, as well as  
maximum atom count of any species in the structure.


#### excitingscripts.plot.poly

Python script for fitting energy-vs-volume curves

##### fit_energy_vs_volume

Fit the energy-vs-volume data using a polynomial fit.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:param order_of_fit: Order of the polynomial fit.  
  
:return:Polynomial fit, bulk modulus, bulk_modulus_pressure_deriv, minima, lattice constant, chi value

##### plot_energy_vs_volume

Plot the energy-vs-volume data and save the plot.  
  
:param volumes : List of volumes.  
:param energies : List of energies.  
:param curv : Polynomial fit curve.  
:param dmin : Minima of the fit.  
:param order_of_polynomial : Order of the polynomial fit.  
:param output_dir : Directory to save the plot.

##### print_info_to_stdout

Print information onto the screen  
  
:param dmin: minima of the fit.  
:param lattice_const: lattice constant for the given lattice symmetry code  
:param bulk_modulus: bulk modulus  
:param bulk_modulus_pressure_deriv: derivative of bulk modulus with respect to pressure  
:param chi: chi-squared value indicating the goodness of fit


#### excitingscripts.plot.relaxdistance

Python visualization tool for following relative atomic coordinates of atoms during the relaxation process.

##### calculate_relative_coordinates

Calculate bond lengths between specified atoms from position data.  
  
:param position_data: list of dictionaries containing atomic positions.  
:param atom1: first atom number.  
:param atom2: second atom number.  
:param lattice_matrix: unit cell containing the molecule  
:param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.  
:param threshold: threshold value for periodic boundary conditions  
  
:return: list of calculated bond lengths.

##### plot_relative_coordinates

Plot the bond lengths and save the plot to a file.  
  
:param delta_x: list of x-components.  
:param delta_y: list of y-components.  
:param delta_z: list of z-components.  
:param run_dir: directory where exciting runs.  
:param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.  
:param ymax:  
:param ymin:  
:param show: whether to display the plot.  
:param dpi: resolution in DPI for the saved plot.


#### excitingscripts.plot.spectra

##### plot_spectra

Plot imaginary part of the macroscopic dielectric function.  
  
:param plot_file_path: Path to file containing data wanted for plot.  
:param color_index: Index needing for plots with different colors for each calculation.


#### excitingscripts.plot.spintext

Produce plot of the spin texture.  
  
Located at `excitingscripts/plot/spintext.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.spintext -b ib -c context   
```  
Where <code>ib</code> defines the band index for the plot and <code>context</code> defines the context of the contour plot. Choises are <code>energy</code> and <code>spin_z</code>.

##### plane_transformation

Take reciprocal lattice vectors and ONS of a plane in rec. lat. coordinates where the first two vectors span the  
plane and the third is normal to them and calculate a matrix that transforms points in the plane to the xy plane  
in cartesian coordinates.  
  
:param rec_lat_vec: Reciprocal lattice vectors.  
:param plot_vec: ONS of the plotting plane.  
:return transformation_matrix: Matrix that transforms k and spin vectors to the plot plane.

##### plot_spintext

Plot the spin texture for a given band.  
  
:param root_directory: Directory of the exciting calculation.  
:param band:      Number of the band for the plot.  
:param contour:   Variable that will be plotted as contour. Can be either energy or spin_z.  
:param contour_threshold: Threshold for the contour plit. Can be either max or float.  
If max, the threshold is the absolute maximum value of the contour.

##### reciprocal_lattice_vectors

Get the reciprocal lattice vectors of real-space lattice vectors \{\mathbf{a}\}:  
  
\mathbf{b}_0 = 2 \pi \frac{\mathbf{a}_1 \wedge \mathbf{a}_2} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}  
\mathbf{b}_1 = 2 \pi \frac{\mathbf{a}_2 \wedge \mathbf{a}_3} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}  
\mathbf{b}_2 = 2 \pi \frac{\mathbf{a}_0 \wedge \mathbf{a}_1} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}  
  
:param lat_vec: Lattice vectors, stored column-wise  
:return: rec_lat_vec: Reciprocal lattice vectors, stored column-wise

##### triple_product

Vector triple product, defined as \mathbf{a} \cdot (\mathbf{b} \wedge \mathbf{c}).  
  
:param a: Vector a  
:param b: Vector b  
:param c: Vector c  
:return Triple product


#### excitingscripts.plot.status

Python visualization tool for the RMS deviations of the SCF potential as a function of the iteration  
number during the SCF loop.

##### plot_status

Python visualization tool for the RMS deviations of the SCF potential as a function  
 of the iteration number during the SCF loop.  
  
:param run_dir: directory where exciting runs


#### excitingscripts.plot.torque

Python visualization tool for the total torque during relaxation.

##### plot_torque

Plot the torque components and magnitude, then save the plot to a file.  
  
:param torque_data: list of torque data to plot.  
:param run_dir: directory where exciting runs.  
:param show: whether to display the plot.  
:param dpi: resolution in DPI for the saved plot.  
:param tol: determines the lowest value possible


#### excitingscripts.plot.vinet

Python script for fitting energy-vs-volume curves

##### fit_energy_vs_volume

Fit the energy-vs-volume data using the Vinet Equation of State.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
  
:return:Vinet fit, bulk modulus, bulk_modulus_pressure_deriv, minima, lattice constant, chi value

##### plot_energy_vs_volume

Plot the energy-vs-volume data and save the plot.  
  
:param volumes : List of volumes.  
:param energies : List of energies.  
:param curv : Polynomial fit curve.  
:param dmin : Minima of the fit.  
:param output_dir : Directory to save the plot.

##### print_info_to_stdout

Print information onto the screen  
  
:param dmin: minima of the fit.  
:param lattice_const: lattice constant for the given lattice symmetry code  
:param bulk_modulus: bulk modulus  
:param bulk_modulus_pressure_deriv: derivative of bulk modulus with respect to pressure  
:param chi: chi-squared value indicating the goodness of fit

##### vinet_eos

Vinet equation of state.  
  
:param v: Volume  
:param eq_vol: Equilibrium volume  
:param min_energy: Minimum energy  
:param bulk_modulus: Bulk modulus at equilibrium volume  
:param bulk_modulus_pressure_deriv: Pressure derivative of bulk modulus  
:return: Energy at volume v


#### excitingscripts.plot.volumecurves

Fit energy-vs-volume curves.  
  
Located at `excitingscripts/plot/volumecurves.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.plot.volumecurves -r dir1 dir2 dir3  
```  
Where <code>dir1</code>, <code>dir2</code>, <code>dir3</code> take the place of the names of the  directories where  
exciting calculations have been performed. The script can be used for any number of directories.

##### determine_functional

Determine the name of the XC functional in a given input file.  
  
:param input_file: Input file.  
:returns: String containing the name of the XC functional.

##### plot_volumecurves

Plot energy-vs-curve values for a given running directory.  
  
:param root_directory: Root directory.  
:param color_index: Index needing for plotting curves with different colors for each calculation.


### excitingscripts.setup


#### excitingscripts.setup.band_structure

Add band structure element to given input file by getting the band path from the input structure.

##### setup_band_structure

Add band structure element to given input file by getting the band path from the input structure.  
  
:param input_file: Input file.  
:param root_directory: Root directory.


#### excitingscripts.setup.convergence_test

Generate input files with different values of the main computational parameters.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.setup.convergence_test k_i k_f rgkmax_i rgkmax_f  
```  
Where <code>k_i</code> and <code>k_f</code> are the initial and final k-values for defining the <code><span style="color:green">groundstate</span></code> attribute  <code><span style="color:MediumBlue">ngridk</span></code>, and <code>rgkmax_i</code> and <code>rgkmax_f</code> the  initial and final values for the  <code><span style="color:green">groundstate</span></code> attribute  <code><span style="color:MediumBlue">rgkmax</span></code>.

##### setup_convergence_test

Create input files with varying values for the groundstate attributes ngridk and rgkmax and save them in  
corresponding directories.  
  
    :param input_file: Input file.  
    :param k_initial: Initial k-value for defining the groundstate attribute ngridk.  
    :param k_final: Final k-value for defining the groundstate attribute ngridk.  
    :param rgkmax_initial: Initial value for the groundstate attribute rgkmax.  
    :param rgkmax_final: Final value for the groundstate attribute rgkmax.  
    :param root_directory: Root directory.


#### excitingscripts.setup.dft_05

Generate a set of input files varying the attribute <code><span style="color:mediumblue">cut</span></code>.  
  
Located at `excitingscripts/setup/dft_05.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.setup.dft_05 r_cut_min r_cut_max number_r_cut_steps -s species -r root_dir  
```  
Where <code>r_cut_min</code> and <code>r_cut_max</code> are the minimum and maximum values for r_cut,  <code>number_r_cut_steps</code> is the number of r_cut values for which input files are generated, <code>species</code> is the species with regard to which r_cut is varied and <code>root_dir</code> is the root directory.

##### setup_dft_05

Generate a set of input files varying the attribute "cut".  
  
:param input_file: Input file.  
:param r_cut_min: Minimum r_cut value.  
:param r_cut_max: Maximum r_cut value.  
:param number_r_cut_steps: Number of r_cut values for which input files are generated by varying the attribute  
"cut".  
:param species: Species used for varying r_cut.  
:param root_directory: Root directory.


#### excitingscripts.setup.diamond_phonon

Diamond phonon setup script.

##### SetupFunc

Protocol for diamond phonon setup function, providing a flexible way of performing type checking.

##### arg_parser

Get the arg parser for phonon setup scripts.  
  
:param point: the point, gamma or x  
:return: the argparser

##### point_specific_setup

Setup function for diamond phonon calculations.  
  
:param get_new_positions: function to get the new positions  
:param set_supercell: how to set a supercell, only for X phonons  
:return: setup function


#### excitingscripts.setup.diamond_phonon_g

Script to set up a phonon calculation for diamond for a phonon at the gamma point.

##### arg_parser

Get the arg parser for phonon setup scripts.  
  
:param point: the point, gamma or x  
:return: the argparser

##### get_new_positions

Get the new displaced positions for a phonon with Gamma character.  
  
:param displacement: displacement of the atom(s)  
:param equilibrium_positions: the non-displaced positions  
:param phonon_mode: not used here  
:return: array with the new positions

##### point_specific_setup

Setup function for diamond phonon calculations.  
  
:param get_new_positions: function to get the new positions  
:param set_supercell: how to set a supercell, only for X phonons  
:return: setup function


#### excitingscripts.setup.diamond_phonon_x

Script to set up a phonon calculation for diamond at the X point.

##### arg_parser

Get the arg parser for phonon setup scripts.  
  
:param point: the point, gamma or x  
:return: the argparser

##### get_new_positions

Get the new displaced positions for a phonon with X character.  
  
:param d: displacement of the atoms  
:param equilibrium_positions: the non-displaced positions  
:param phonon_mode: Mode of the phonon  
:return: array with the new positions

##### point_specific_setup

Setup function for diamond phonon calculations.  
  
:param get_new_positions: function to get the new positions  
:param set_supercell: how to set a supercell, only for X phonons  
:return: setup function

##### set_supercell

Set a supercell according to an x phonon.  
  
:param input_obj: the input xml object


#### excitingscripts.setup.dos_band_structure

Add DOS and band structure element to given input file by getting the band path from the input structure.  
  
Located at `excitingscripts/setup/dos_band_structure.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.setup.dos_band_structure  
```

##### setup_dos_band_structure

Add DOS and band structure element to given input file by getting the band path from the input structure.  
  
:param input_file: Input file.  
:param root_directory: Root directory.


#### excitingscripts.setup.elastic_strain

Python script for generating strained structures.

##### get_deformation_mapping

Returns the deformation code string indicating the type of strains in Voigt notation  
  
:param deformation_code:  
:return :

##### get_deformation_matrix

Determines the deformation matrix from Langrangian strain matrix  
  
:param eta_matrix: langrangian strain matrix  
:param max_iter:  
:param tol:  
:return : 3x3 numpy array of deformation matrix

##### get_langrangian_strain_matrix

Determines the Langrangian strain matrix based on the deformation string provided  
  
:param eta: strain value  
:param deformation_str: string determining type deformation  
:return : 3x3 numpy array of langrangian strain matrix

##### setup_deformed_structures

Create input files with a series volume values for different strain values.  
  
:param input_file: Input file.  
:param maximum_strain:  
:param number_strain_values: Number of strain values for which structures are generated by varying the lattice constant.  
:param deformation_code:  
:param workdir: Working directory.


#### excitingscripts.setup.excitingroot

Replace placeholder "$EXCITINGROOT" in **input.xml** files by actual path.  
  
Located at `excitingscripts/setup/excitingroot.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.setup.excitingroot  
```

##### set_exciting_root

Replace all instances of the string '$EXCITINGROOT' in the file given by `input_file`  
and write to `output_file`.  
  
:param input_file: Input file.  
:param output_file: Input file, with '$EXCITINGROOT' replaced with excitingroot.  
:param excitingroot: Environment variable string.


#### excitingscripts.setup.graphene_along_c

Python script for generating strained structures.

##### setup_displaced_structures

:param input_file: Input file.  
:param eta_strain: Lagrangian strain.  
:param umax: Maximum displacement.  
:param number_of_displ: Number of displacements.  
:param workdir: Working directory.


#### excitingscripts.setup.interlayer_distance

Generate structures with different interlayer distances.  
  
Located at `excitingscripts/setup/interlayer_distance.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.setup.interlayer_distance dmin dmax nr_displ dinfty  
```  
Where <code>dmin</code> and <code>dmax</code> are the minimum and maximum values for the interlayer distance, <code>nr_displ</code> is the number of distances in the interval [<code>dmin</code>, <code>dmin</code>] and <code>dinfty</code> is the interlayer distance at infinity.

##### setup_interlayer_distance

Create input files for structures with different interlayer distances and save them in corresponding directories.  
  
:param input_file: Input file.  
:param dmin: Minimum interlayer distance in Bohr.  
:param dmax: Maximum interlayer distance in Bohr.  
:param displ_points: Number of distances in [dmin, dmax].  
:param dinfty: Interlayer distance at infinity in Bohr.  
:param root_directory: Root directory.


#### excitingscripts.setup.planar

Python script for setting calculations for two-dimensional materials.

##### setup_planar

Creates the planar files  
  
:param workdir: Working directory.


#### excitingscripts.setup.volume_optimization

Generate structures at different volumes.  
  
Located at `excitingscripts/setup/volume_optimization.py`.  
  
Call as:  
  
```bash  
python3 -m excitingscripts.setup.volume_optimization nr_vol  
```  
Where <code>nr_vol</code> is the number of volume values for which structures are generated by varying the lattice constant.

##### setup_volume_optimization

Create input files with a series volume values for structures generated at equally spaced intervals of the  
lattice constant with a variation of between -5% and +5% from the reference lattice constant and save them in  
corresponding directories.  
  
:param input_file: Input file.  
:param number_volume_values: Number of volume values for which structures are generated by varying the lattice  
constant.  
:param root_directory: Root directory.


### excitingscripts.utils


#### excitingscripts.utils.utils

General utils for exciting scripts.

##### T1



##### T2



##### birch_murnaghan_eos

Calculate energy using the Birch-Murnaghan equation of state.  
  
:param v: Volume.  
:param p: Parameters for the equation [eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv].  
:return: Calculated energy.

##### birch_murnaghan_fit

Perform the least squares fit using the Birch-Murnaghan EOS.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:return: Optimized parameters.

##### extract_values_from_line

Extract all numbers from a given line using regular expressions.  
  
:param line: input string from which to extract numbers.  
:return: list of values found in the input string.

##### get_decimal_decomposition

Decompose the number into mantissa and exponent.  
  
:param number: input number  
:return: tuple with shifted number (only one leading digit before the decimal point) and exponent

##### get_num_atoms

Extract the total number of atoms per unit cell from INFO.OUT.  
  
:param run_dir: directory where exciting runs.  
:return: number of atoms per unit cell.

##### get_prettified_scientific_notation

Decompose the number into mantissa and exponent and produce formatted string.  
  
:param number: input number  
:param unit: unit of the number  
:return: prettified string representation

##### get_structure_optimizations_properties

Read all lines from the INFO.OUT file, extract property for each optimization step.  
  
:param run_dir: directory where exciting runs.  
:param key: property name which is parsed for each optimization step. Available ones are:  
            "Maximum force",  
            "Center of mass",  
            "Total torque",  
            "Number of total scf iterations",  
            "Total atomic forces",  
            "Total energy",  
            "Atomic positions"  
  
:return: list of dictionaries containing properties.

##### initial_guess

Generate initial guess parameters for the Birch-Murnaghan EOS fit.  
  
:param volumes: List of volumes.  
:param energies: List of energies.  
:return: Initial guess parameters [eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv].

##### is_coordinate_cartesian

Check the coordinate type is cartesian from input.xml.  
  
:param run_dir: directory where exciting runs  
:return: coordinate type, either True for "cartesian" or False for "lattice" or other type.

##### parse_energy_vs_volume

Read the "energy_vs_volume" file  
  
:param directory: directory containing the file  
:return: volume and energy values

##### pressure_birch_murnaghan_eos

Calculate pressure using the Birch-Murnaghan equation of state.  
  
:param v: Volume.  
:param p: Parameters for the equation [eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv].  
:return: Calculated pressure.

##### residuals

Calculate residuals for the least squares fit.  
  
:param p: Parameters for the Birch-Murnaghan EOS.  
:param e: Energies.  
:param v: Volumes.  
:return: Residuals.

##### sort_lists_by_first_list

Sorts two lists, using the first list as reference  
  
:param first_list: first list to be sorted, used as reference  
:param second_list: second list to be sorted, uses first list as reference  
:return: sorted lists


