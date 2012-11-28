#!/bin/sh
if [ ! -e /bin/csh ]; then
  echo /bin/csh cannot be found. Please install it or create a symbolic link from /bin/csh to /usr/bin/tcsh , if installed.
  exit 1
fi
./foolproof $*
