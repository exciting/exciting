# <span style="color:#4056A1">0. Installation of the excitingJupyter Tutorial Environment</span>
**<span style="color:firebrick">Read the following paragraphs before starting with the tutorials!</span>**

<div style="text-align: justify">

## Compilation of exciting

Before starting, be sure that **`exciting`** is already compiled according to the procedure reported in 
**[<span style="color:#D79922">Download and compile exciting</span>](http://exciting.wikidot.com/neon-download-and-compile-exciting)**. 
This is also documented in exciting's `INSTALL` file in the repository root.

## Running the Wikidot Tutorial Scripts

In order to run the scripts used in the **`exciting`** tutorials it is important that the relevant environment variables 
are already defined in your **~/.bashrc** file as specified in **[<span style="color:#D79922">How to set environment 
variables for tutorials scripts</span>](http://exciting.wikidot.com/neon-tutorial-scripts-and-environment-variables)**. 
However, this is not necessary for running the tutorials inside Jupyter Notebooks.

## Installing the excitingJupyter Package

All Jupyter tutorials require Python 3.7 or above to run.  As a first step in running the Jupyter tutorials, it is 
useful to create a virtual environment in which you can install and run the notebooks. Each venv has its own Python binary 
(which matches the version of the binary that was used to create this environment) and can have its own independent set 
of installed Python packages in its site directories. To create a venv, move to:

```bash
cd $EXCITINGROOT/tools/excitingjupyter
```

You can create the venv using an executable bash script:

```bash
source create_env.sh
```

This script will automatically activate the environment, so you can start right ahead. 

Alternatively, you can do it by hand.

### Manually installing excitingjupyter

```bash
# Create python venv for running excitingjupyter
mkdir venv && cd venv
python3 -m venv excitingvenv
source excitingvenv/bin/activate
cd ..

python3 -m pip install --upgrade --force pip
pip3 install --upgrade setuptools
# Install excitingtools 
pip3 install -e ../exciting_tools
# Install excitingscripts
pip3 install -e ../excitingscripts
# Install excitingjupyter
pip3 install .
# Install local kernal for jupyter
python3 -m ipykernel install --user --name=excitingjupyter
```

If you see the error `invalid command: bdist_wheel` when installing excitingtools, you need to run `pip3 install wheel`, 
then try again.

One can leave the venv at any time by typing `deactivate`. If you are repeating the installation procedure for the venv, 
remember to exit with `deactivate` before trying to regenerate it. Creating a venv for the exciting notebooks to run in 
is only required once, however you should ensure it is activated **every** time you wish to run a notebook. From 
exciting's root:

```bash
source $EXCITINGROOT/tools/excitingjupyter/venv/excitingvenv/bin/activate
```
### Adding Custom CSS Style
In order to add the layout style designed for the Jupyter tutorials, type the following commands starting from 
`tools/excitingjupyter`:

```bash
# Find path for custom CSS file:
path=$(python -c "import notebook; print(notebook.__file__)")
notebookpath=${path::-11}
csspath="${notebookpath}static/custom/."
# Add custom CSS style:
cp excitingjupyter/custom.css "$csspath"
cp ../../docs/logo/logotransp.png "$csspath"
```
## Runtime Libraries

Please note that exciting requires certain libraries at runtime (for the SMP version, openBLAS or MKL), and setting them 
in a terminal shell is not sufficient as each Jupyter cell creates a new shell instance. The easiest way to ensure they 
are present is to add them to your `.bashrc`. For example, the SOL group uses the TCL module system so one would add:

```bash
module load intel/2019
```

to the `.bashrc` (which loads everything required). **Please take an equivalent approach on your platform.**

It is required to do this _in the same terminal where Jupyter is started_, before starting with the tutorials.

## Starting Jupyter

To start working with the Jupyter notebooks, move to:

```bash
cd $EXCITINGROOT/tools/excitingjupyter/excitingjupyter
```

and execute:

```bash
jupyter-notebook
```

This will open your browser, where you can select the tutorial you want to work on.
To start, _e.g._, with the first tutorials, click on the folder **01_getting_started**,
and select the notebook: **tutorial_how_to_start_an_exciting_calculation.ipynb**.

This should launch an executable version of the notebook in a new tab of your browser.

</div>
