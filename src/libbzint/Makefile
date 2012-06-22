include ../../make.inc


DESTDIR  = ./
LIBNAME = libbzint.a
ARCH     = ar
ARCHFLAGS= -rc
###ARCHFLAGS= crv
RANLIB   = ar s
DIST_FILE =$(DESTDIR)/src_libbzint.tgz


FFLAGS = $(FGEN) $(FOPT)
LIB   = $(DESTDIR)$(LIBNAME)

#
#  Files defining modules
#
MODFILES = order.f90 kgen_internals.f90 tetra_internal.f90 polyhedron.f90 

#
#  Files for generation of k-point grids
#
KGENFILES = kgen.f90 kgenq.f90 kqgen.f90 kqgen_exciting.f90

#
# Subroutines for use of kgenfiles only
# 

KGENSUBS =  coorskp.f90 divisi.f90 jget.f90 idkp.f90 \
            jset.f90  sym2int.f90 tetinit.f90 tgen.f90 \
            tgenq.f90

#
#  Files for BZ integrations and convolutions
#

INTFILES = tetcw.f90 tetcw_1k.f90 tetcw_1kbc.f90 tetcorecw.f90 tetiw.f90 tetiwsurf.f90 tetraint.f90 \
           tetraqint.f90 calcdosve.f90 calcidosve.f90 dostet.f90 idos.f90 fermitet.f90 

#
# Subroutine for use of intfiles only
#

INTSUBS = bloechlcor.f90 convw1t.f90 convw.f90 convcorew.f90 dos1t.f90 \
	convw_1k.f90 convw_1kbc.f90 \
          generictetra.f90 genericfunf.f90 genericprism.f90 intdos1t.f90 intweight1t.f90 \
          intw.f90 intwsurf.f90 ksurf.f90 relnodes.f90 setnodes.f90 sortnodes.f90 \
          sortsurf.f90 surfnodes.f90 tlinap.f90 unrepnodes.f90 sorteq.f90 redifwt.f90 \
          redifwtx.f90 redifwty.f90 redifwtz.f90 edifwt.f90 edifwtaylor.f90 edifwtx.f90 edifwty.f90 \
          edifwtz.f90 redifwtaylor.f90 reduk.f90 redusym.f90

#           
#  Files that can be called independently (and are also used by the others programs)
#
FREEFILES = cartezian.f90 gbass.f90 intern.f90 rbass.f90 factorize.f90

#
# All files
#

FILES = $(MODFILES) $(KGENFILES) $(INTFILES) $(FREEFILES) $(KGENSUBS) $(INTSUBS)
#
OBJS = $(FILES:.f90=.o)

#
.SUFFIXES:	.f90 .mod

#..............................................................................
#
#  Build executable
#
##$(LIB):	$(OBJS)
##	$(ARCH) $(ARCHFLAGS) $(LIB) $(OBJS)
##	$(RANLIB) $(LIB)

$(LIB):	$(OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB) $(OBJS)

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#

$(OBJS): $(MODFILES)

#..............................................................................
#
#  remove object files, preprocessed source files and files needed for
#  static semantics checking
#
clean:  cleanmod
	rm  -f $(OBJS) *.a  

cleanmod:
	rm  -f *.mod

#..............................................................................
#
#  generate latex documentation
doc:
	cp tetra.tex $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\section{Routine/Function Prologues}' >> $(LIB:.a=.tex)
	protex -b -n $(MODFILES) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines for the generation of k-point grids}' >> $(LIB:.a=.tex)
	protex -b -n $(KGENFILES) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines used by the k-point generation subroutines}' >> $(LIB:.a=.tex)
	protex -b -n $(KGENSUBS) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines for the Brillouin Zone integration}' >> $(LIB:.a=.tex)
	protex -b -n $(INTFILES) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Subroutines used by the integration subroutines}' >> $(LIB:.a=.tex)
	protex -b -n $(INTSUBS) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\subsection{Other useful subroutines that can also be called separatedly}' >> $(LIB:.a=.tex)
	protex -b -n $(FREEFILES) >> $(LIB:.a=.tex)
	@echo '\newpage' >> $(LIB:.a=.tex)
	@echo '\input{tests}' >> $(LIB:.a=.tex)
	@echo '\end{document}' >> $(LIB:.a=.tex)

#..............................................................................
#
#  generate distribution file
#
dist:   
	tar -cPzvf $(DIST_FILE) Makefile $(FILES) tetra.tex tests.tex 
ifeq ($(BUILDSMBUILDSMPP),true)
F90=$(F90MT)
endif

.f90.o: Makefile
	$(F90) $(F90_OPTS) -c $<
