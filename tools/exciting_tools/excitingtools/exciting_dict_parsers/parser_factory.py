"""Function for selecting a parser, given a file name and parsing the data.

When adding a new file, it should be added to the dictionary ` _file_to_parser`

REQUIREMENTS. Parser function must:
 a) accept a file name (not a contents string)
 b) return a dictionary.
"""
import os

from excitingtools.utils.dict_utils import container_converter

from excitingtools.exciting_dict_parsers import \
    bse_parser, groundstate_parser, gw_eigenvalues_parser, gw_eps00_parser, gw_info_parser, gw_vxc_parser, \
    input_parser, properties_parser, RT_TDDFT_parser, species_parser, state_parser


# Map file name to parser function
_file_to_parser = {
    'INFO.OUT': groundstate_parser.parse_info_out,
    'info.xml': groundstate_parser.parse_info_xml,
    'input.xml': input_parser.parse_input_xml,
    'species.xml': species_parser.parse_species_xml,
    'atoms.xml': groundstate_parser.parse_atoms,
    'evalcore.xml': groundstate_parser.parse_evalcore,
    'eigval.xml': groundstate_parser.parse_eigval,
    'geometry.xml': groundstate_parser.parse_geometry,
    'LINENGY.OUT': groundstate_parser.parse_linengy,
    'LO_RECOMMENDATION.OUT': groundstate_parser.parse_lo_recommendation,
    'RHO3D.xml': properties_parser.parse_plot_3d,
    'VCL3D.xml': properties_parser.parse_plot_3d,
    'VXC3D.xml': properties_parser.parse_plot_3d,
    'WF3D.xml': properties_parser.parse_plot_3d,
    'ELF3D.xml': properties_parser.parse_plot_3d,
    'EF3D.xml': properties_parser.parse_plot_3d,
    'LSJ.xml': properties_parser.parse_lsj,
    'EFG.xml': properties_parser.parse_efg,
    'mossbauer.xml': properties_parser.parse_mossbauer,
    'expiqr.xml': properties_parser.parse_expiqr,
    'effmass.xml': properties_parser.parse_effmass,
    'bandstructure.xml': properties_parser.parse_bandstructure_depreciated,
    'dos.xml': properties_parser.parse_dos,
    'KERR.OUT': properties_parser.parse_kerr,
    'EPSILON_11.OUT': properties_parser.parse_epsilon,
    'EPSILON_12.OUT': properties_parser.parse_epsilon,
    'EPSILON_33.OUT': properties_parser.parse_epsilon,
    'CHI_111.OUT': properties_parser.parse_chi,
    'ELNES.OUT': properties_parser.parse_elnes,
    'SEEBECK_11.OUT': properties_parser.parse_seebeck,
    'ELECTCOND_11.OUT': properties_parser.parse_seebeck,
    'THERMALCOND_11.OUT': properties_parser.parse_seebeck,
    'Z_11.OUT': properties_parser.parse_seebeck,
    'ldos.out': properties_parser.parse_ldos,
    'band_edges.out': properties_parser.parse_band_edges,
    'spintext.xml': properties_parser.parse_spintext,
    'POLARIZATION.OUT': properties_parser.parse_polarization,
    'TDOS_WANNIER.OUT': properties_parser.parse_tdos_wannier,
    'WANNIER_INFO.OUT': properties_parser.parse_wannier_info,
    'coreoverlap.xml': properties_parser.parse_core_overlap,
    'INFOXS.OUT': bse_parser.parse_infoxs_out,
    'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_FXCMB1_OC11_QMT001.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_FXCMB1_OC22_QMT001.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_FXCMB1_OC33_QMT001.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_NLF_FXCMB1_OC11_QMT001.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_NLF_FXCMB1_OC22_QMT001.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_NAR_NLF_FXCMB1_OC33_QMT001.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC12.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_BSE-singlet-TDA-BAR_SCR-full_OC12.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC12.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EPSILON_NAR,
    'DICHROIC_K_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC12.OUT': bse_parser.parse_EPSILON_NAR,
    'OSCI_K_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'OSCI_K_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'OSCI_K_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'OSCI_K_BSE-singlet-TDA-BAR_SCR-full_OC12.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'OSCI_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'OSCI_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'OSCI_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'OSCI_KBAR_BSE-singlet-TDA-BAR_SCR-full_OC12.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EPSILON_BSE-IP_SCR-full_OC11.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_BSE-IP_SCR-full_OC22.OUT': bse_parser.parse_EPSILON_NAR,
    'EPSILON_BSE-IP_SCR-full_OC33.OUT': bse_parser.parse_EPSILON_NAR,
    'LOSS_NAR_FXCMB1_OC11_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_FXCMB1_OC22_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_FXCMB1_OC33_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC11_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC22_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC33_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-IP_SCR-full_OC11.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-IP_SCR-full_OC22.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-IP_SCR-full_OC33.OUT': bse_parser.parse_LOSS_NAR,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-IP_SCR-full_OC11.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-IP_SCR-full_OC22.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-IP_SCR-full_OC33.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'GW_INFO.OUT': gw_info_parser.parse_gw_info,
    'EFERMI_GW.OUT': gw_eigenvalues_parser.parse_efermi_gw,
    'EVALQP.DAT': gw_eigenvalues_parser.parse_evalqp,
    'VXCNN.DAT': gw_vxc_parser.parse_vxcnn,
    'EPS00_GW.OUT': gw_eps00_parser.parse_eps00_gw,
    'JIND.OUT': RT_TDDFT_parser.parse_jind,
    'NEXC.OUT': RT_TDDFT_parser.parse_nexc,
    'ETOT_RTTDDFT.OUT': RT_TDDFT_parser.parse_etot,
    'EIGVAL_': RT_TDDFT_parser.parse_eigval_screenshots,
    'PROJ_': RT_TDDFT_parser.parse_proj_screenshots,
    'ATOM_': RT_TDDFT_parser.parse_atom_position_velocity_force,
    'FCR_': RT_TDDFT_parser.parse_force,
    'FEXT_': RT_TDDFT_parser.parse_force,
    'FHF_': RT_TDDFT_parser.parse_force,
    'FVAL_': RT_TDDFT_parser.parse_force,
    'wf1d-0001-0001.dat': properties_parser.parse_wf1d,
    'wf1d-0003-0001.dat': properties_parser.parse_wf1d,
    'wf2d-0001-0001.xsf': properties_parser.parse_wf2d,
    'wf2d-0003-0001.xsf': properties_parser.parse_wf2d,
    'wf3d-0001-0001.xsf': properties_parser.parse_wf3d,
    'wf3d-0003-0001.xsf': properties_parser.parse_wf3d,
    'wf3d-0001-0001.cube': properties_parser.parse_cube,
    'wf3d-0003-0001.cube': properties_parser.parse_cube,
    'STATE.OUT': state_parser.parse_state_out,
}


def truncate_fnames_with_exts(file_name: str) -> str:
    """ Truncate file names that have open-ended extensions.

    For example:
      EIGVAL_00.dat -> EIGVAL_
      EIGVAL_01.dat -> EIGVAL_

    :param file_name: File name containing fixed prefix and
    an extension beginning with '_'.
    :return file_name: File name prefix, else input file name.
    """
    prefixes = [
        'EIGVAL_',
        'PROJ_',
        'ATOM_',
        'FCR_',
        'FEXT_',
        'FHF_',
        'FVAL_'
    ]
    if any( [ prefix in file_name for prefix in prefixes ] ) :
        file_name_prefix = file_name.split('_')[0] + '_'
        return file_name_prefix

    return file_name


def parser_chooser(full_file_name: str) -> dict:
    """ Selects parser according to the name of the input file then returns the result of the parser.

    REQUIREMENTS. Parser function must:
     a) accept a file name (not a contents string)
     b) return a dictionary.

    param: str, full_file_name: file name prepended by full path
    return: parsed data
    """
 
    full_file_name = full_file_name.rstrip()
    if not os.path.exists(full_file_name):
        raise FileNotFoundError(f'File not found: {full_file_name}')

    file_name = os.path.split(full_file_name)[1]
    file_name = truncate_fnames_with_exts(file_name)

    files_with_parsers = [name for name in _file_to_parser.keys()]
    if file_name not in files_with_parsers:
        raise KeyError(f"File does not have a parser: {file_name}")

    parser = _file_to_parser[file_name]
    data = parser(full_file_name)

    #  TODO(Alex) Issue 135 Ensure all parsers return appropriate values, not strings.
    #   container_converter should therefore be used as a decorator on parsers with values that are strings.
    #   That will massively speed up parsing
    return container_converter(data)

