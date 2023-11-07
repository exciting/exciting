# ``excitingjupyter`` - Jupyter notebook based tutorials for the **`exciting`** code

This package implements interactive tutorials for the **``exciting``** code.

One design goal of this package is to allow for double functionality: The tutorials should be executable via Jupyter and also provide a stand-alone version for the [exciting-code](exciting-code.org) website.

This requires that all functionalities that need custom Python code have to be implemented as scripts that are executable from the command line. 
Examples and instructions can be found in `$EXCITINGROOT/tools/excitingscripts`. 
Before implementing something new there, check if a more general function already exists and can be adapted.

## Installation

It is highly recommended to use a vitual environment for executing the **``exciting``** tutorials, or building the website source-code. If your virtual environment does not exist yet, you can create it with:

```bash
cd $EXCITINGROOT/tools/excitingjupyter
source create_env.sh
```

If it already exists, you can activate it with:

```bash
source venv/excitingvenv/bin/activate
```

More explanations can be found in [`excitingjupyter/00_before_starting.md`](excitingjupyter/00_before_starting.md).

## Run the interactive version

(LOAD ANY REQUIRED MODULES BEFORE STARTING JUPYTER)

After creating the virtual environment you can run Jupyter notebook with:

```bash
jupyter-notebook
```

This should automatically open a new tab in your browser, showing a GUI with the contents of the current folder. If not, the output in your terminal should give you all instructions to find a solution.

In the GUI, navigate to: `excitingjupyter`, and consequently the sub folders with the topic of your interest (_e.g._ `01_getting_started`). In there, you find several files with the `.ipynb` extension. Click on them, to start the interactive tutorial. In the running version, you can do (and save) changes.

## Create static `html` contents 

To generate the `html` pages, make sure you are in the `excitingjupyter` root-directory and run:

```bash
jupyter nbconvert --execute --config static_web_config.py
```

NOTE: This may take some time because all tutorials will be executed!

(for development purposes just omit the `--execute` flag.)

### Changing the template

Take a look at [the documentation](https://nbconvert.readthedocs.io/en/latest/customizing.html) of `nbconvert`. The templates are located at `$EXCITINGROOT/tools/excitingjupyter/templates`. The default template for `excitingjupyter`s `html` output (as defined in `./static_web_config.py`) is `static_web` (because it produces the source code for a static website). To edit the `.css` files (which should be about the only required changes to the template), please refer to `./templates/static_web/static/excitingjupyter.css`.
