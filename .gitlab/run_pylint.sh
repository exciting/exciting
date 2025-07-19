#!/usr/bin/env bash

# Run pylint in Gitlab continuous integration env, and return an error
# if pylint returns text.
#
# Notes
# ----------
#  Should be run from the repository root.
#  Accepts current branch as a script argument.
#  It might make sense to have pylint output the report in JSON and have a 
#  python script evaluate it, replacing this script (could also unit test it).

reference_branch="development"
current_branch=$1
changed_files=$(git diff --name-only origin/${reference_branch} origin/"${current_branch}" -- '*.py')

echo "Running pylint on python files with diff w.r.t. ${reference_branch}:"

# This works locally, but not in the bash shell of the CI
# therefore procecss the file names on the fly (see below)
# # Split string of changed files w.r.t. white space
# py_files=(${changed_files// / })
# # Apply pylint to python files
# for file in "${py_files[@]}"
# ...

# Apply pylint to python files
for file in $(echo "$changed_files" | tr " " "\n")
do
   # Skip files in any 'external' folder
   # This is to avoid issues with packages containing python
   if [[ "$file" == external/* || "$file" == */external/* ]]; then
    echo "Skipping external file: $file"
    continue
  fi

  echo "Checking file: $file"
  # F0001 disables import errors. This is a workaround for when files
  # have been deleted between commits (and there's nothing 
  # check in that instance)
  # E0402 checks for relative imports beyond top-level package
  # E0401 checks more imports
  # E0001 are syntax errors, needed for the old python2 scripts
  error_msg=$(pylint -E --disable=E0611 --disable=E1120 --disable=F0001 --disable=E0402 --disable=E0401 --disable=E0001 "$file")
  # pylint returns nothing if there are no errors
  if [ -n "$error_msg" -a "$error_msg" != " " ]; then
      echo "pylint has experienced an error:"
      echo "$error_msg"
      exit 1
  fi
done

echo "pylint has found no errors in files."
exit 0
