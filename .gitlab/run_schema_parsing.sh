#!/usr/bin/env bash

# Run the schema parsing script from excitingtools,
# to check whether the schema has changed and excitingtools stays up-to-date
# Note: Should be run from the repository root.

echo "Checking excitingtools is consistent with exciting's schema ..."

# Need to copy the old file out for reference
cp tools/exciting_tools/excitingtools/utils/valid_attributes.py .

if ! python3 -m excitingtools.utils.schema_parsing; then
  rm -f valid_attributes.py
  exit 1
fi

difference=$(diff valid_attributes.py tools/exciting_tools/excitingtools/utils/valid_attributes.py)
# clean up
mv valid_attributes.py tools/exciting_tools/excitingtools/utils/

if [ -n "$difference" ]; then
  echo "The schema has changed:"
  echo "$difference"
  echo ""
  echo "A discrepancy in the parsed schema was found. That usually means that you have touched the schema of exciting, "
  echo "and not updated the schema's representation in excitingtools."
  echo ""
  echo "In order to do so, go in the exciting root directory, upgrade pip and xmlschema, install excitingtools in "
  echo "editable mode and run the schema_parsing with the following commands:"
  echo "cd <your_development_exciting_root_dir>"
  echo "python3 -m pip install --upgrade pip xmlschema"
  echo "python3 -m pip install -e tools/exciting_tools"
  echo "python3 -m excitingtools.utils.schema_parsing"
  echo ""
  echo "NOTE: If you are using python3.7, you need an older version of xmlschema via:"
  echo "python3 -m pip install xmlschema==2.5.1"
  echo ""
  exit 1
fi

echo "Passed: excitingtools is up-to-date with exciting's schema."
exit 0