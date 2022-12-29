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
  echo "Checking file: $file"
  error_msg=$(pylint -E "$file")
  # pylint returns nothing if there are no errors
  if [ -n "$error_msg" -a "$error_msg" != " " ]; then
      echo "pylint has experienced an error:"
      echo "$error_msg"
      exit 1
  fi
done

echo "pylint has found no errors in files."
exit 0
