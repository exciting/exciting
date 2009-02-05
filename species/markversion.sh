#! /bin/sh
# (C) 2009 S. Sagmeister and C. Ambrosch-Draxl
VERSION=`awk -F"/" '/data version/ {print $2}' ../src/mod_misc.F90 | sed 's/ //g;s/,/./g'`
stamp="# Exciting code version : $VERSION"

for fil in *.in; do
	signed=`grep 'Exciting code version' $fil | wc -l`
	if [ $signed -eq 0 ]; then
		echo >> $fil
		echo $stamp >> $fil
	else
		echo "file $fil already marked with Exciting version "`grep 'Exciting code version' $fil | awk -F ':' '{print $2}'`
	fi
done

