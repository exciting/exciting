#! /bin/sh

#
# purpose of this script: remove exciting binary output files
#

if [ $# -ne 1 ]; then
	echo
	echo "Syntax: `basename $0` list_of_filename_patterns"
	echo
	exit 1
fi

# with #-commented list of filenames
flist="$1"

outlist="exciting-binary-outputfiles.txt"
rmf="/bin/rm -f"

# remove comments (#) and blank lines
flist=`cat $flist | grep -v "^ *#" | awk 'NF>0'`

# remove old list
if [ -f $outlist ]; then
	$rmf $outlist
fi

echo
echo "Checking for the following binary output files:"
echo

# generate list
touch $outlist
for files in $flist; do
	echo "$files"
	# we also remove files like "YYYY.OUT*" - note the asterisk at the end!
	find . -type f -name "$files*" >> $outlist
done


echo
echo "Exciting binary output files written to: $outlist"
echo
