"""
Test bse parsers.
The infoxs reference is boiled down to the necessary parts for a test.
Execute tests from exciting_tools directory:
pytest --capture=tee-sys
"""
import pytest
from excitingtools.exciting_dict_parsers.bse_parser import parse_infoxs_out

infoxs_file_str_success = """================================================================================
| EXCITING NITROGEN-14 started for task xsgeneigvec (301)                      =
| version hash id: 1775bff4453c84689fb848894a9224f155377cfc                    =
|                                                                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
================================================================================
One-shot GS runs for BSE calculations
================================================================================
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
One-shot GS runs for k+qmt/2 grids
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Info(xsgeneigvec): Generating eigenvectors for Q-point      1
--------------------------------------------------------------------------------
Info(xsgeneigvec): Generation of eigenvectors finished
 
Info(xsfinit): task Nr.     301 stopped gracefully
 
 Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    301                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task writepmatxs (320)                      =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
Info(writepmatxs): Momentum matrix elements finished
Info(xsfinit): task Nr.     320 stopped gracefully
 
 Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:24
     CPU time               : 4.46 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 04 s )
     wall time              : 0.84 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 01 s )
     CPU load               : 533.19 %
     CPU time  (cumulative) : 116.04 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 56 s )
     wall time (cumulative) : 2.96 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 03 s )
     CPU load  (cumulative) : 533.19 %
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    320                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task scrgeneigvec (401)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
 
Info(xsinit): mapping screening-specific parameters
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    401                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task scrwritepmat (420)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
 
Info(xsinit): mapping screening-specific parameters
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    420                                 =
================================================================================
 
================================================================================
| EXCITING NITROGEN-14 started for task bse (445)                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
 
Info(xsinit): mapping BSE-specific parameters

================================================================================
= EXCITING NITROGEN-14 stopped for task    445                                 =
================================================================================
"""

infoxs_file_str_fail = """================================================================================
| EXCITING NITROGEN-14 started for task xsgeneigvec (301)                      =
| version hash id: 1775bff4453c84689fb848894a9224f155377cfc                    =
|                                                                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
================================================================================
One-shot GS runs for BSE calculations
================================================================================
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
One-shot GS runs for k+qmt/2 grids
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Info(xsgeneigvec): Generating eigenvectors for Q-point      1
--------------------------------------------------------------------------------
Info(xsgeneigvec): Generation of eigenvectors finished

Info(xsfinit): task Nr.     301 stopped gracefully

 Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %

================================================================================
= EXCITING NITROGEN-14 stopped for task    301                                 =
================================================================================

================================================================================
| EXCITING NITROGEN-14 started for task writepmatxs (320)                      =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
Info(writepmatxs): Momentum matrix elements finished
Info(xsfinit): task Nr.     320 stopped gracefully

================================================================================
| EXCITING NITROGEN-14 started for task xsgeneigvec (301)                      =
| version hash id: 1775bff4453c84689fb848894a9224f155377cfc                    =
|                                                                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
================================================================================
One-shot GS runs for BSE calculations
================================================================================
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
One-shot GS runs for k+qmt/2 grids
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Info(xsgeneigvec): Generating eigenvectors for Q-point      1
--------------------------------------------------------------------------------
Info(xsgeneigvec): Generation of eigenvectors finished

Info(xsfinit): task Nr.     301 stopped gracefully

 Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %

================================================================================
= EXCITING NITROGEN-14 stopped for task    301                                 =
================================================================================


================================================================================
| EXCITING NITROGEN-14 started for task writepmatxs (320)                      =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
Info(writepmatxs): Momentum matrix elements finished
Info(xsfinit): task Nr.     320 stopped gracefully
 
 Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:24
     CPU time               : 4.46 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 04 s )
     wall time              : 0.84 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 01 s )
     CPU load               : 533.19 %
     CPU time  (cumulative) : 116.04 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 56 s )
     wall time (cumulative) : 2.96 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 03 s )
     CPU load  (cumulative) : 533.19 %
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    320                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task scrgeneigvec (401)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
 
Info(xsinit): mapping screening-specific parameters
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    401                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task scrwritepmat (420)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
"""

reference_parsed_infoxs_file_success = {
    'tasks': [{'name': 'xsgeneigvec', 'number': 301,
               'finished': True}, {'name': 'writepmatxs', 'number': 320,
                                   'finished': True},
              {'name': 'scrgeneigvec', 'number': 401,
               'finished': True}, {'name': 'scrwritepmat', 'number': 420,
                                   'finished': True},
              {'name': 'bse', 'number': 445, 'finished': True}],
    'success': True,
    'last_finished_task': 'bse'
    }

reference_parsed_infoxs_file_fail = {
    'tasks': [{'name': 'xsgeneigvec', 'number': 301,
               'finished': True}, {'name': 'writepmatxs', 'number': 320,
                                   'finished': False},
              {'name': 'xsgeneigvec', 'number': 301,
               'finished': True}, {'name': 'writepmatxs', 'number': 320,
                                   'finished': True},
              {'name': 'scrgeneigvec', 'number': 401,
               'finished': True}, {'name': 'scrwritepmat', 'number': 420,
                                   'finished': False}],
    'success': False,
    'last_finished_task': 'scrgeneigvec'
    }


@pytest.mark.parametrize(["infoxs_file_str", "reference_parsed_dict"],
                         [(infoxs_file_str_success,
                           reference_parsed_infoxs_file_success),
                          (infoxs_file_str_fail,
                           reference_parsed_infoxs_file_fail)])
def test_parse_info_xs_out(infoxs_file_str, reference_parsed_dict, tmp_path):
    infoxs_file_path = tmp_path / "INFOXS.OUT"
    infoxs_file_path.write_text(infoxs_file_str)
    info_xs_out = parse_infoxs_out(infoxs_file_path.as_posix())
    assert info_xs_out == reference_parsed_dict
