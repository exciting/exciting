# <span style="color:#4056A1">0. Installation of the excitingJupyter Tutorial Environment</span>
**<span style="color:firebrick">Read the following paragraphs before starting with the tutorials!</span>**

<div style="text-align: justify">

## Getting exciting

This guide explains how to download and unpack the latest version of the **exciting** code. If you have already done it, please skip to the next section.

### 1. Choose a Directory

Choose a directory where you want to install or unpack the current release of *exciting sodium*. You can use any directory where you have write permissions. You can obtain *exciting sodium* using one of the two following methods. **Choose either Step 2 or Step 3** based on your preference.

### 2. Download and Unpack the Release Archive

After downloading the release archive file (i.e., `exciting-sodium.tar.gz`), extract it using:

```bash
tar -xvf exciting-sodium.tar.gz
```

This will unpack the contents into a directory named `exciting-sodium`.

### 3. Cloning from GitHub

Alternatively, you can clone the latest version of **exciting sodium** from the official GitHub repository:

```bash
git clone https://github.com/exciting/exciting -b sodium
```

This will download the contents into a directory named `exciting` containing the `sodium` release.

## Compilation of exciting

Before starting, be sure that **`exciting`** is already compiled according to the procedure reported in 
**[<span style="color:#D79922">Download and compile exciting</span>](https://exciting-code.org/uploads/exciting/tutorial_notebooks/00_tutorial_download_and_compile_exciting.html)**. 
This is also documented in exciting's `INSTALL` file in the repository root.

## Setting Environment Variables 

Before running the tutorials, it is important to set the necessary environment variables.
To do this, move to the exciting root directory and run the following commands:
```bash
cd tools/excitingjupyter
source set_env_vars.sh
```

## Installing the excitingJupyter Package

All Jupyter tutorials require Python 3.10 or above to run. As a first step in running the Jupyter tutorials, it is
useful to create a virtual environment (venv) in which you can install and run the notebooks. Each venv has its own Python binary 
(which matches the version of the binary that was used to create this environment) and can have its own independent set 
of installed Python packages in its site directories. To create a venv, move to:

```bash
cd $EXCITINGROOT/tools/excitingjupyter
```

You can create the venv using an executable bash script:

```bash
source create_env.sh
```

This script will automatically activate the environment for you, so you can start right ahead. 

Alternatively, you can do it by hand by following the instruction given below.

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
remember to exit with `deactivate` before trying to regenerate it. Creating a venv for the exciting notebooks 
is required only once, however you must ensure it is activated **every** time you wish to run a notebook. From 
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
notebookpath=${path::-20}
csspath="${notebookpath}nbclassic/static/custom/."
# Add custom CSS style:
rm -f $csspath/custom.css
cp excitingjupyter/custom.css "$csspath"
cp ../../docs/logo/logotransp.png "$csspath"
```
## Runtime Libraries

Please note that exciting requires certain libraries at runtime (for the SMP version, openBLAS or MKL), and setting them 
in a terminal shell is not sufficient as each Jupyter cell creates a new shell instance. The easiest way to ensure they 
are present is to add them to your `.bashrc`. For example, the SOL group uses the TCL module system so one would add:

```bash
module load intel-oneapi/2025.1
```

to the `.bashrc` (which loads everything required). **Please take an equivalent approach on your platform, and note that
one needs to load the same modules used for the compilation step.**

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
and select the notebook: **how_to_start_an_exciting_calculation.ipynb**.

This should launch an executable version of the notebook in a new tab of your browser.

</div>
