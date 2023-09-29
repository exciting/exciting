import os

import pygments
from nbconvert.preprocessors import TagRemovePreprocessor

NB_ROOT = os.path.join(os.environ["EXCITINGROOT"], "tools", "excitingjupyter", "excitingjupyter")
OUTDIR = os.path.join(os.environ["EXCITINGROOT"], "tools", "excitingjupyter", "static_html")

os.makedirs(OUTDIR, exist_ok=True)

c = get_config()  # type:ignore # pylint: disable=E0602

c.NbConvertApp.notebooks = [
    # Convert all tutorials in 01_getting_started and 02_ground_state (may be slow if combined with the --execute flag)
    os.path.join(NB_ROOT, "01_getting_started", "*_tutorial_*.ipynb"),
    os.path.join(NB_ROOT, "02_ground_state", "*_tutorial_*.ipynb")
    # for fast testing use, e.g., (because it has images etc.):
    # os.path.join(NB_ROOT, "01_getting_started", "03_tutorial_electronic_band_structure_and_density_of_states.ipynb")
]
c.NbConvertApp.export_format = "html"

c.TemplateExporter.extra_template_basedirs.append(
    os.path.join(os.environ["EXCITINGROOT"], "tools", "excitingjupyter", "templates"))
c.TemplateExporter.template_name = "static_web"

c.CSSHTMLHeaderPreprocessor.style = pygments.styles.get_style_by_name("default")

c.FilesWriter.build_directory = OUTDIR


class AfterExecutionTagRemovePreprocessor(TagRemovePreprocessor):
    """Custom preprocessors are called _after_ the default preprocessors.
    Thus, if we define a copy of the default class, we can just use it in the later loop.
    """

    remove_cell_tags = {"remove_cell"}
    remove_all_outputs_tags = {"remove_output"}


c.Exporter.preprocessors = [AfterExecutionTagRemovePreprocessor]

## you probably want to check:: $EXCITINGROOT/tools/excitingjupyter/venv/excitingvenv/share/jupyter
## customize: https://nbconvert.readthedocs.io/en/latest/customizing.html ,  write client: https://nbconvert.readthedocs.io/en/latest/nbconvert_library.html
