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
  echo "In order to do so, install excitingtools and xmlschema with the following commands:"
  echo "python3 -m pip install xmlschema"
  echo "python3 -m pip install \$EXCITINGROOT/tools/exciting_tools"
  echo "python3 -m excitingtools.utils.schema_parsing"
  echo ""
  exit 1
fi

echo "Passed: excitingtools is up-to-date with exciting's schema."
exit 0