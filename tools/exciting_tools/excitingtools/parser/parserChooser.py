"""
API for selecting a parser, given a file name and parsing the data

When adding a new file, it should be added to the dictionary ` _file_to_parser`

If the parser takes the read-in file string rather than the file name,
it should be added to the function `parser_expects_file_str`.
"""
import os
import sys
from xml.etree.ElementTree import ParseError
from typing import Callable

from excitingtools.dict_utils import container_converter

from . import groundstate_parser
from . import properties_parser
from . import bse_parser
from . import gw_parser
from . import RT_TDDFT_parser
from .erroneous_file_error import ErroneousFileError
from .parser_utils import generic_parser


# Map file name to parser function
_file_to_parser = {
    'INFO.OUT': groundstate_parser.parse_info_out,
    'info.xml': groundstate_parser.parse_info_xml,
    'atoms.xml': groundstate_parser.parse_atoms,
    'evalcore.xml': groundstate_parser.parse_evalcore,
    'eigval.xml': groundstate_parser.parse_eigval,
    'geometry.xml': groundstate_parser.parse_geometry,
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
    'bandstructure.xml': properties_parser.parse_bandstructure,
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
    'LOSS_NAR_FXCMB1_OC11_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_FXCMB1_OC22_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_FXCMB1_OC33_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC11_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC22_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_NAR_NLF_FXCMB1_OC33_QMT001.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_LOSS_NAR,
    'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC33.out': bse_parser.parse_LOSS_NAR,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT': bse_parser.parse_EXCITON_NAR_BSE,
    'GW_INFO.OUT': gw_parser.parse_gw_info,
    'EFERMI_GW.OUT': gw_parser.parse_efermi_gw,
    'EVALQP.DAT': gw_parser.parse_evalqp,
    'VXCNN.DAT': gw_parser.parse_vxcnn,
    'EPS00_GW.OUT': gw_parser.parse_eps00_gw,
    'JIND.OUT': RT_TDDFT_parser.parse_jind,
    'NEXC.OUT': RT_TDDFT_parser.parse_nexc,
    'ETOT_RTTDDFT.OUT': RT_TDDFT_parser.parse_etot,
    'EIGVAL_': RT_TDDFT_parser.parse_eigval_screenshots,
    'PROJ_': RT_TDDFT_parser.parse_proj_screenshots
}


def parser_expects_file_str(file_name: str) -> bool:
    """
    Distinguish between parsers that expect the filename and those which expect
    the read-in file as a string.

    Could be implemented as a dictionary {file_name: bool} but whilst the
    number of functions that accept file strings instead of file names is small,
    this results in less code.

    :param str file_name: Name of exciting file.
    :return bool: If the parser function expects parsed file string as the argument.
    """
    # Files with parsers expecting the read-in file string as input
    parsers = ['GW_INFO.OUT', 'EPS00_GW.OUT']
    if file_name in parsers:
        return True
    else:
        return False


def parser_chooser(full_file_name: str) -> dict:
    """
    Selects parser according to the name of the input file
    then returns the result of the parser.

    Could probably implement with decorators.

    param: str, full_file_name: file name prepended by full path
    return: parsed data
    """
    full_file_name = full_file_name.rstrip()
    file_name = os.path.split(full_file_name)[1]

    if ('EIGVAL_' in file_name) or ('PROJ_' in file_name):
        file_name = file_name.split('_')[0] + '_'

    files_with_parsers = [name for name in _file_to_parser.keys()]
    if not (file_name in files_with_parsers):
        sys.exit(f"File does not have a parser: {file_name}")

    parser = _file_to_parser[file_name]

    if parser_expects_file_str(file_name):
        return generic_parser(full_file_name, parser)

    else:
        if not os.path.exists(full_file_name):
            raise OSError(f'File path not valid: {full_file_name}')
        data = parser(full_file_name)
        return container_converter(data)
