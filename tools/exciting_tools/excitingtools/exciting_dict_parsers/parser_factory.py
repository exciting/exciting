"""Function for selecting a parser, given a file name and parsing the data.

When adding a new file, it should be added to the dictionary ` _file_to_parser`

REQUIREMENTS. Parser function must:
 a) accept a file name (not a contents string)
 b) return a dictionary.
"""

import warnings
from fnmatch import fnmatch
from pathlib import Path
from typing import Callable, Union

from excitingtools.exciting_dict_parsers import (
    RT_TDDFT_parser,
    bse_parser,
    groundstate_parser,
    gw_eigenvalues_parser,
    gw_eps00_parser,
    gw_info_parser,
    gw_taskgroup_parser,
    gw_vxc_parser,
    hdf5_parser,
    input_parser,
    properties_parser,
    species_parser,
    state_parser,
)

# Map file name to parser function
# Note: more specific names should be higher, as the search will go through this map top-down
_file_to_parser = {
    "INFO.OUT": groundstate_parser.parse_info_out,
    "info.xml": groundstate_parser.parse_info_xml,
    "input.xml": input_parser.parse_input_xml,
    "species.xml": species_parser.parse_species_xml,
    "atoms.xml": groundstate_parser.parse_atoms,
    "evalcore.xml": groundstate_parser.parse_evalcore,
    "eigval.xml": groundstate_parser.parse_eigval,
    "geometry.xml": groundstate_parser.parse_geometry,
    "LINENGY.OUT": groundstate_parser.parse_linengy,
    "LO_RECOMMENDATION.OUT": groundstate_parser.parse_lo_recommendation,
    "*3D.xml": properties_parser.parse_plot_3d,
    "LSJ.xml": properties_parser.parse_lsj,
    "EFG.xml": properties_parser.parse_efg,
    "mossbauer.xml": properties_parser.parse_mossbauer,
    "expiqr.xml": properties_parser.parse_expiqr,
    "effmass.xml": properties_parser.parse_effmass,
    "bandstructure.xml": properties_parser.parse_bandstructure_depreciated,
    "dos.xml": properties_parser.parse_dos,
    "KERR.OUT": properties_parser.parse_kerr,
    "EPSILON_??.OUT": properties_parser.parse_epsilon,
    "CHI_111.OUT": properties_parser.parse_chi,
    "ELNES.OUT": properties_parser.parse_elnes,
    "SEEBECK_11.OUT": properties_parser.parse_seebeck,
    "ELECTCOND_11.OUT": properties_parser.parse_seebeck,
    "THERMALCOND_11.OUT": properties_parser.parse_seebeck,
    "Z_11.OUT": properties_parser.parse_seebeck,
    "ldos.out": properties_parser.parse_ldos,
    "band_edges.out": properties_parser.parse_band_edges,
    "spintext.xml": properties_parser.parse_spintext,
    "POLARIZATION.OUT": properties_parser.parse_polarization,
    "TDOS_WANNIER.OUT": properties_parser.parse_tdos_wannier,
    "WANNIER_INFO.OUT": properties_parser.parse_wannier_info,
    "coreoverlap.xml": properties_parser.parse_core_overlap,
    "wf1d-*.dat": properties_parser.parse_wf1d,
    "wf2d-*.xsf": properties_parser.parse_wf2d,
    "wf3d-*.xsf": properties_parser.parse_wf3d,
    "wf3d-*.cube": properties_parser.parse_cube,
    "INFOXS.OUT": bse_parser.parse_infoxs_out,
    "EPSILON_BSE*.OUT": bse_parser.parse_EPSILON_NAR,
    "EPSILON_NAR*.OUT": bse_parser.parse_EPSILON_NAR,
    "DICHROIC_*.OUT": bse_parser.parse_EPSILON_NAR,
    "OSCI_*.OUT": bse_parser.parse_EXCITON_NAR_BSE,
    "EXCITON_*.OUT": bse_parser.parse_EXCITON_NAR_BSE,
    "LOSS_*.OUT": bse_parser.parse_LOSS_NAR,
    "GW_INFO.OUT": gw_info_parser.parse_gw_info,
    "EFERMI_GW.OUT": gw_eigenvalues_parser.parse_efermi_gw,
    "EVALQP.DAT": gw_eigenvalues_parser.parse_evalqp,
    "VXCNN.DAT": gw_vxc_parser.parse_vxcnn,
    "EPS00_GW.OUT": gw_eps00_parser.parse_eps00_gw,
    "BARC_*": gw_taskgroup_parser.parse_barc,
    "SGI_*": gw_taskgroup_parser.parse_sgi,
    "EPSILON-GW_Q*": gw_taskgroup_parser.parse_epsilon,
    "EPSH.OUT": gw_taskgroup_parser.parse_epsilon,
    "EPSW1.OUT": gw_taskgroup_parser.parse_epsilon,
    "EPSW2.OUT": gw_taskgroup_parser.parse_epsilon,
    "INVERSE-EPS*": gw_taskgroup_parser.parse_inverse_epsilon,
    "JIND.OUT": RT_TDDFT_parser.parse_jind,
    "NEXC.OUT": RT_TDDFT_parser.parse_nexc,
    "ETOT_RTTDDFT.OUT": RT_TDDFT_parser.parse_etot,
    "EIGVAL_*": RT_TDDFT_parser.parse_eigval_screenshots,
    "PROJ_*": RT_TDDFT_parser.parse_proj_screenshots,
    "ATOM_*": RT_TDDFT_parser.parse_atom_position_velocity_force,
    "FCR_*": RT_TDDFT_parser.parse_force,
    "FEXT_*": RT_TDDFT_parser.parse_force,
    "FHF_*": RT_TDDFT_parser.parse_force,
    "FVAL_*": RT_TDDFT_parser.parse_force,
    "STATE.OUT": state_parser.parse_state_out,
    "bse_output.h5": hdf5_parser.parse_hdf5_file_as_dict,
    "fastBSE_output.h5": hdf5_parser.parse_hdf5_file_as_dict,
    "fastBSE_absorption_spectrum.out": bse_parser.parse_fastBSE_absorption_spectrum_out,
    "fastBSE_exciton_energies.out": bse_parser.parse_fastBSE_exciton_energies_out,
    "fastBSE_oscillator_strengths.out": bse_parser.parse_fastBSE_oscillator_strength_out,
}


def parse(full_file_name: str) -> dict:
    """Selects parser according to the name of the input file then returns the result of the parser.

    REQUIREMENTS. Parser function must:
     a) accept a file name (not a contents string)
     b) return a dictionary.

    :param full_file_name: file name prepended by full path
    :return: parsed data
    """

    full_file_path = Path(full_file_name.rstrip())
    if not full_file_path.exists():
        raise FileNotFoundError(f"File not found: {full_file_path}")

    file_name = full_file_path.name

    parser: Union[Callable[[str], dict], None] = None
    for pattern in _file_to_parser:
        if fnmatch(file_name, pattern):
            parser = _file_to_parser[pattern]
            break

    if not parser:
        raise KeyError(f"File does not have a parser: {file_name}")

    return parser(full_file_path.as_posix())


def parser_chooser(full_file_name: str) -> dict:
    """Old API. Selects parser according to the name of the input file then returns the result of the parser.

    :param full_file_name: file name prepended by full path
    :return: parsed data
    """
    warnings.warn(
        "Deprecated API. Use 'excitingtools.parse' instead. "
        "Support for this API will be removed in excitingtools 1.8.0",
        DeprecationWarning,
        stacklevel=2,
    )
    return parse(full_file_name)
