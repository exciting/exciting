#! /bin/sh

rm -r ./debian/usr
mkdir -p ./debian/usr
mkdir -p ./debian/usr/share
mkdir -p ./debian/usr/share/exciting/doc
mkdir -p ./debian/usr/bin
cat version control.part2 >./debian/DEBIAN/control

cp ../docs/exciting/excitinginput.pdf 	\
../docs/spacegroup/spacegroup.pdf \
../docs/exciting/excitingsubroutines.pdf \
../docs/Brillouin/* \
./debian/usr/share/exciting/doc

cp -r ../examples ./debian/usr/share/exciting/examples
cp -r ../species ./debian/usr/share/exciting/species

cp ../bin/* ./debian/usr/bin/
ln -s ./debian/usr/bin/excitingser ./debian/usr/bin/exciting
cp ../COPYING ./debian/usr/share/exciting/copyright

chmod -R  a-s ./debian/
dpkg-deb --build debian exciting.deb
