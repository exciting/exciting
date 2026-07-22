# Create python venv for running excitingjupyter
mkdir venv && cd venv || exit 1
python3 -m venv excitingvenv || exit 1
source excitingvenv/bin/activate || exit 1
cd .. || exit 1
python3 -m pip install --upgrade --force pip || exit 1
pip3 install --upgrade setuptools || exit 1
# Install excitingtools 
pip3 install -e ../exciting_tools || exit 1
# Install excitingscripts
pip3 install -e ../excitingscripts || exit 1
# Install excitingjupyter
pip3 install -e . || exit 1
# Install FSvisual
pip3 install --upgrade fsvisual || exit 1
# Install local kernal for jupyter
python3 -m ipykernel install --user --name=excitingjupyter || exit 1
# Find path for custom CSS file:
path=$(python -c "import notebook; print(notebook.__file__)")
notebookpath=${path::-20}
csspath="${notebookpath}nbclassic/static/custom/."
# Add custom CSS style:
rm -f $csspath/custom.css || exit 1
cp excitingjupyter/custom.css "$csspath" || exit 1
cp ../../docs/logo/logotransp.png "$csspath" || exit 1
