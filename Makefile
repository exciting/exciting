
.NOTPARALLEL:

default: build/make.inc all

# Include platform specific variable settings
-include build/make.inc

info:
	@echo ""
	@echo "___ Starting compilation _________________________________________________"
	@echo ""
	@sleep 3

# exciting binary  

all: info mpiandsmp smp

serial:
	cd build/serial; $(MAKE)

smp:
	cd build/smp; $(MAKE)

mpi:
	cd build/mpi; $(MAKE)

mpiandsmp:
	cd build/mpiandsmp; $(MAKE)

debug:
	cd build/debug; $(MAKE)

debugmpiandsmp:
	cd build/debugmpiandsmp; $(MAKE)

# Test suite 

test ::
	cd test/; $(MAKE)

test_mpi ::
	cd test/; $(MAKE) mpi

test_mpiandsmp ::
	cd test/; $(MAKE) mpiandsmp

test_reference ::
	cd test/; $(MAKE) reference 

cleantests ::
	cd test/; $(MAKE) cleantests

# HTML Documentation, limited to a subsection of the source 

# Ford documentation. Note, the document build directory (1st argument) is relative
# to the location of ford_settings.md. Documentation is therefore built in docs/exciting_ford. 
ford :: 
	ford -o exciting_ford/ docs/ford_settings.md

cleanford ::
	cd docs/exciting_ford/;\
	rm -r blockdata css fonts interface js lists module proc program sourcefile src type;\
	rm favicon.png index.html search.html

# Source !BOP documentation and all XML inputs, parsed into latex

doc:  spacegroupdoc stateconvertdoc stateinfodoc inputdoc excitingfuncdoc split_inputdoc speciesdoc

excitingfuncdoc::
	$(MAKE) -f build/Make.common doc

expandedschema::
	xsltproc xml/schema/schemaexpand.xsl xml/schema/input.xsd >xml/excitinginput.xsd

inputdoc::expandedschema
	cd docs/exciting/;\
	xsltproc --stringparam importancelevels "essential expert" ../../xml/schematolatex.xsl ../../xml/excitinginput.xsd >excitinginput.tex;\
	xsltproc --stringparam importancelevels "essential expert" ../../xml/schematowikidot.xsl ../../xml/excitinginput.xsd >inputref.wikidot;\
	pdflatex excitinginput.tex;\
	pdflatex excitinginput.tex

split_inputdoc::
	cd xml/schema && $(MAKE)

inputdocwiki:xml/schema/*.xsd
	cd xml/schema; $(MAKE)

# ----------------------------------

# Utility Programs and their docs
# TODO(Alex) Issue 45. Decouple standalone programs from exciting's src and build system 
# Consider replacing with a single command, `make auxilliary` 

spacegroupdoc::
	cd src/spacegroup; $(MAKE) doc;\
	mv spacegroup.pdf ../../docs/spacegroup
	xsltproc xml/schema/schemaexpand.xsl xml/schema/symmetries.xsd>  xml/sgroupinput.xsd
	cd docs/spacegroup; \
	xsltproc --stringparam importancelevels "spacegroup" ../../xml/schematolatex.xsl ../../xml/sgroupinput.xsd > spacegroupinput.tex; \
	xsltproc --stringparam importancelevels "spacegroup" ../../xml/schematowikidot.xsl ../../xml/sgroupinput.xsd > ../../xml/schema/wiki/spacegroup ;\
	pdflatex spacegroupinput.tex; \
	pdflatex spacegroupinput.tex; \

speciesdoc::
	cd docs/species;\
	xsltproc --stringparam importancelevels "spacegroup" ../../xml/schematolatex.xsl ../../xml/species.xsd > species.tex;\
	pdflatex species.tex;pdflatex species.tex;

stateconvertdoc::
	cd src/stateconvert; $(MAKE) doc;\
	mv stateconvert.pdf ../../docs/stateconvert

stateinfodoc::
	cd src/stateinfo; $(MAKE) doc;\
	mv stateinfo.pdf ../../docs/stateinfo

spacegroup::
	cd src/spacegroup; $(MAKE)

stateinfo::
	cd src/stateinfo; $(MAKE)

stateconvert::
	cd src/stateconvert; $(MAKE)

species::
	cd src/species; $(MAKE)

# ---------------------------

libs:
	cd build/serial; $(MAKE) libs

# TODO(Alex) Issue 45. Decouple standalone programs from exciting's src and build system 
clean:
	cd build/serial; $(MAKE) clean cleanlibs
	cd build/mpi; $(MAKE) clean cleanlibs
	cd build/smp; $(MAKE) clean cleanlibs
	cd build/debug; $(MAKE) clean cleanlibs
	cd build/debugmpiandsmp; $(MAKE) clean cleanlibs
	cd build/mpiandsmp; $(MAKE) clean cleanlibs
	cd src/spacegroup; $(MAKE) clean
	cd src/species; $(MAKE) clean
	cd src/stateinfo; $(MAKE) clean
	cd src/stateconvert; $(MAKE) clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe exdg.*
	rm -f docs/exciting/*
	rm -f docs/spacegroup/*
	cd test; $(MAKE) cleantests

libxcclean:
	cd src/libXC && make clean

tgz::doc
	tar --exclude-from=".gitignore"  --transform 's,^,exciting/,' -c -v -f ./exciting.tar *
	gzip  -f --best ./exciting.tar
	du -h ./exciting.tar.gz
