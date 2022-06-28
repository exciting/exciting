# Create python venv for running excitingjupyter
mkdir venv && cd venv
python3 -m venv excitingvenv
source excitingvenv/bin/activate
cd ..
python3 -m pip install --upgrade --force pip
pip3 install --upgrade setuptools
# Install excitingtools
pip3 install -e ../exciting_tools
# Install excitingjupyter
pip3 install .
# Install local kernal for jupyter
python3 -m ipykernel install --user --name=excitingjupyter
# Find path for custom CSS file:
path=$(python -c "import notebook; print(notebook.__file__)")
notebookpath=${path::-11}
csspath="${notebookpath}static/custom/."
# Add custom CSS style:
cp excitingjupyter/custom.css "$csspath"
cp ../../docs/logo/logotransp.png "$csspath"
