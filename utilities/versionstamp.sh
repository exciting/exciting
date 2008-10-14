#! /bin/sh
hash=`git log | head -1 | awk '{print $2}' | cut -c 1-20`
hash2=`git log | head -1 | awk '{print $2}' | cut -c 21-40`
echo "#define GITHASH \"$hash\"" >../../src/version.inc
echo "#define GITHASH2 \"$hash2\"" >>../../src/version.inc
if [ `git-diff | wc -l` -ne 0 ]; then
	echo "#define LOCALCHG" >> ../../src/version.inc
fi
