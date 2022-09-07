import pytest
from excitingtools.utils.test_utils import MockFile
import json

from .utilities import re_input


@pytest.fixture
def mock_notebook(tmp_path):
    # Notebook in dict form
    nb = {'cells': [{'cell_type': 'markdown', 'metadata': {'collapsed': True, 'pycharm': {'name': '#%% md\n'}},
                     'source': '# How to Start an exciting Calculation\n**by <span style="color:darkgreen">Andris Gulans</span>, '
                               '<span style="color:darkgreen">Jürgen Spitaler</span>, & <span style="color:darkgreen">'
                               'Pasquale Pavone</span>; for [<span style="color:darkgoldenrod">exciting *oxygen*</span>]'
                               '(http://exciting.wikidot.com/oxygen)**\n<hr style="border:2px solid #DDD"> </hr>'},
                    {'cell_type': 'markdown',
                     'source': '<div style="text-align: justify">\n\nFor the very first **exciting** run, '
                               'you will use an already prepared example of an input file that sets up a total-energy '
                               'calculation of diamond. Input files for **exciting** are written in the **XML** format '
                               'and are typically named **input.xml**. The **XML** format allows your data to be written '
                               'in a structured way. Figuratively speaking, an **exciting** input is pretty much like an '
                               'article with its sections and subsections. In case of **XML** data, sections and subsections '
                               'are called <code><span style="color:green">elements</span></code>.\n\n\n```xml\n<input>\n\n   '
                               '<title>Diamond</title>\n\n   <structure speciespath="/home/tutorial/exciting/species">\n\n      '
                               '<crystal scale="6.7274">\n         '
                               '<basevect>0.0   0.5   0.5</basevect>\n         '
                               '<basevect>0.5   0.0   0.5</basevect>\n         '
                               '<basevect>0.5   0.5   0.0</basevect>\n      '
                               '</crystal>\n\n      <species speciesfile="C.xml">\n         '
                               '<atom coord="0.00 0.00 0.00"/>\n         '
                               '<atom coord="0.25 0.25 0.25"/>\n      </species>\n\n   </structure>\n\n   '
                               '<groundstate\n      ngridk="4 4 4"\n      outputlevel="normal"\n      '
                               'xctype="GGA_PBE_SOL">\n   </groundstate>\n\n</input>\n\n```',
                     'metadata': {'collapsed': False, 'pycharm': {'name': '#%% md\n'}}}],
          'metadata': {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                       'language_info': {'codemirror_mode': {'name': 'ipython', 'version': 2}, 'file_extension': '.py',
                                         'mimetype': 'text/x-python', 'name': 'python', 'nbconvert_exporter': 'python',
                                         'pygments_lexer': 'ipython2', 'version': '2.7.6'}}, 'nbformat': 4,
          'nbformat_minor': 0}
    file = tmp_path / "mock_nb.ipynb"
    with file.open("w", encoding="UTF-8") as fid:
        json.dump(nb, fid)
    return MockFile(file, 'dummy string')


# Expected return string
ref_string = """<input>

   <title>Diamond</title>

   <structure speciespath="/home/tutorial/exciting/species">

      <crystal scale="6.7274">
         <basevect>0.0   0.5   0.5</basevect>
         <basevect>0.5   0.0   0.5</basevect>
         <basevect>0.5   0.5   0.0</basevect>
      </crystal>

      <species speciesfile="C.xml">
         <atom coord="0.00 0.00 0.00"/>
         <atom coord="0.25 0.25 0.25"/>
      </species>

   </structure>

   <groundstate
      ngridk="4 4 4"
      outputlevel="normal"
      xctype="GGA_PBE_SOL">
   </groundstate>

</input>
"""


def string_assert(output: str, reference: str):
    """ Compare two Strings.
    String comparison is very brittle.
    Compare line by line, ignoring leading and trailing
    white space.
    """
    ref_list = reference.split()
    output_list = output.split()
    assert len(ref_list) == len(output_list)

    for i in range(len(ref_list)):
        assert output_list[i].strip() == ref_list[i].strip()


def test_re_input(mock_notebook):
    input_title = "Diamond"
    input_str = re_input(mock_notebook.full_path, input_title)
    string_assert(input_str, ref_string)
