#!/usr/bin/env bash

# Create sim link for exciting binary.
# Expects a single argument of the form exciting$(SUFFIX)
# implying that it's executed from the bin directory.
# A Buccheri March 2021   

input_args=$1
input_array=(${input_args//_/ })
input_array_length=${#input_array[@]}

exciting_binary=${input_array[0]}
build_type=${input_array[${#input_array[@]}-1]}

if [ "$build_type" = "mpismp" ] && [ $input_array_length -lt 3 ]; then
    ln -s exciting_mpismp exciting
    echo "Creating symbolic link: exciting -> exciting_mpismp"
fi

if [ "$build_type" = "mpismp" ] && [ $input_array_length == 3 ]; then
    ln -s exciting_debug_mpismp exciting_debug
    echo "Creating symbolic link exciting_debug -> exciting_debug_mpismp"
fi
