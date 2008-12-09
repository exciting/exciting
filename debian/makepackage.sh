#! /bin/sh


mkdir -p ./debian/usr/share/doc/exciting/
mkdir -p ./debian/usr/bin/

cp ../docs/exciting/exciting.pdf 	\
../docs/spacegroup/spacegroup.pdf \
../docs/Brillouin/* \
./debian/usr/share/doc/exciting/

cp ../bin/* ./debian/usr/bin/
cp ../COPYING ./debian/usr/share/doc/exciting/copyright

chmod a-s -R debian/
dpkg-deb --build debian exciting.deb
