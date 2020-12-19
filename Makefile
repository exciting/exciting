
.NOTPARALLEL:

default: build/make.inc all

# Include platform specific variable settings (will be passed on to other make calls)
-include build/make.inc

info:
	@echo ""
	@echo "___ Starting compilation _________________________________________________"
	@echo ""
	@sleep 3

all: info serial mpiandsmp 

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

test::
	cd test/; $(MAKE) summary

doc: inputdoc split_inputdoc speciesdoc
#TODO(Sebastian) Issue #15. Fix documentation compilation for subroutines 
#doc:  spacegroupdoc stateconvertdoc stateinfodoc inputdoc excitingfuncdoc split_inputdoc speciesdoc

excitingfuncdoc::
	$(MAKE) -f build/Make.common doc

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

expandedschema::
	xsltproc xml/schema/schemaexpand.xsl xml/schema/input.xsd >xml/excitinginput.xsd ;\

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

stateconvertdoc::
	cd src/stateconvert; $(MAKE) doc;\
	mv stateconvert.pdf ../../docs/stateconvert

stateinfodoc::
	cd src/stateinfo; $(MAKE) doc;\
	mv stateinfo.pdf ../../docs/stateinfo

eos::
	cd src/eos; $(MAKE)

spacegroup::
	cd src/spacegroup; $(MAKE)

stateinfo::
	cd src/stateinfo; $(MAKE)

stateconvert::
	cd src/stateconvert; $(MAKE)

species::
	cd src/species; $(MAKE)

libs:
	cd build/serial; $(MAKE) libs

clean:
	cd build/serial; $(MAKE) clean cleanlibs
	cd build/mpi; $(MAKE) clean cleanlibs
	cd build/smp; $(MAKE) clean cleanlibs
	cd build/debug; $(MAKE) clean cleanlibs
	cd build/debugmpiandsmp; $(MAKE) clean cleanlibs
	cd build/mpiandsmp; $(MAKE) clean cleanlibs
	cd src/eos; $(MAKE) clean
	cd src/spacegroup; $(MAKE) clean
	cd src/species; $(MAKE) clean
	cd src/src_vdwdf; $(MAKE) clean
	cd src/stateinfo; $(MAKE) clean
	cd src/stateconvert; $(MAKE) clean
	cd src/leblaiklib; $(MAKE) clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe exdg.*
	rm -f interfaces/*
	rm -f docs/exciting/*
	rm -f docs/spacegroup/*
	cd test; $(MAKE) cleantests

libxcclean:
	cd src/libXC && make clean

tgz::doc #libxcclean
	tar --exclude-from=".gitignore"  --transform 's,^,exciting/,' -c -v -f ./exciting.tar *
	gzip  -f --best ./exciting.tar
	du -h ./exciting.tar.gz

vdwdf:
	cd src/src_vdwdf
	$(MAKE)
