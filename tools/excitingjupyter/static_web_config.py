import os

import pygments
from nbconvert.preprocessors import TagRemovePreprocessor

NB_ROOT = os.path.join(os.environ["EXCITINGROOT"], "tools", "excitingjupyter", "excitingjupyter")
OUTDIR = os.path.join(os.environ["EXCITINGROOT"], "tools", "excitingjupyter", "static_html")

os.makedirs(OUTDIR, exist_ok=True)

c = get_config()  # type:ignore # pylint: disable=E0602  # noqa: F821

c.NbConvertApp.notebooks = [
    # Convert all tutorials (may be slow if combined with the --execute flag)
    os.path.join(NB_ROOT, "01_getting_started", "*.ipynb"),
    os.path.join(NB_ROOT, "02_ground_state/01_Methods", "*.ipynb"),
    os.path.join(NB_ROOT, "02_ground_state/02_Electronic_Properties", "*.ipynb"),
    os.path.join(NB_ROOT, "02_ground_state/03_Lattice_Optimization", "*.ipynb"),
    os.path.join(NB_ROOT, "02_ground_state/04_Molecules", "*.ipynb"),
    os.path.join(NB_ROOT, "02_ground_state/05_Lattice_Dynamics", "*.ipynb"),
    os.path.join(NB_ROOT, "02_ground_state/06_Elastic_Properties", "*.ipynb"),
    os.path.join(NB_ROOT, "03_excited_states/01_GW", "*.ipynb"),
    os.path.join(NB_ROOT, "03_excited_states/02_BSE", "*.ipynb"),
    os.path.join(NB_ROOT, "03_excited_states/03_TDDFT", "*.ipynb"),
    os.path.join(NB_ROOT, "03_excited_states/04_EPH", "*.ipynb"),
    os.path.join(NB_ROOT, "03_excited_states/05_Others", "*.ipynb"),
    os.path.join(NB_ROOT, "04_additional_features", "*.ipynb"),
    os.path.join(NB_ROOT, "05_tools_and_packages", "*.ipynb"),
    # for fast testing use, e.g., (because it has images etc.):
    # os.path.join(NB_ROOT, "01_getting_started", "electronic_band_structure_and_density_of_states.ipynb")
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
