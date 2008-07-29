ref=`cat  ../../.git/HEAD |sed s/ref://`
ref=../../.git/$ref
ref=`echo $ref | sed s/\\\s//`
hash=`cat $ref | cut -c 1-20`
echo \#define GITHASH \"$hash\" >../../src/version.inc