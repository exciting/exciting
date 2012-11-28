#!/bin/csh

echo char \*helpstring=\"\"
sed 's/^/\"/g' | sed 's/$/\\n\"/g'
echo \;
