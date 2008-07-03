ref=`cat  ../../.git/HEAD |sed s/ref://`
hash=`grep $ref  ../../.git/info/refs | cut -c 1-20`
echo \#define GITHASH \"$hash\" >../../src/version.inc