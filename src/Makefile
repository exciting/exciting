include ../make.inc

MAKE = make

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .F90
.F90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC_modules = modmain.F90 modmpi.F90  modxcifc.F90 

SRC_main = main.F90

SRC_routines = \
 autoradmt.F90 findprim.F90 findsym.F90 findsymcrys.F90 genppts.F90 rfinp.F90 \
 rvfcross.F90 seceqn.F90 symrf.F90 wavefmt.F90

SRC_routines_flux = \
 elfplot.F90 rhonorm.F90 energy.F90 spinchar.F90 \
 zpotclmt.F90 writegeom.F90 nfftifc.F90 zfftifc.F90 allatoms.F90 gridsize.F90 \
 poteff.F90 genrmesh.F90 readfermi.F90 potcoul.F90 gensfacgp.F90 checkmt.F90 \
 zfinp.F90 match.F90 force.F90 forcek.F90 writeefg.F90 packeff.F90 \
 bandchar.F90 findsymlat.F90 \
 genlofr.F90 atom.F90 writefermi.F90 writekpts.F90 fsmfield.F90 mossbauer.F90 \
 occupy.F90 writelinen.F90 writeinfo.F90 readinput.F90 charge.F90 moment.F90 \
 writesym.F90 genidxlo.F90 gencore.F90 addrhocr.F90 gengvec.F90 \
 genshtmat.F90 plot1d.F90 plot2d.F90 plot3d.F90 updatpos.F90 writeiad.F90 \
 symvect.F90 vecplot.F90 genylmg.F90 linengy.F90 init0.F90 init1.F90 \
 gengpvec.F90 genpmat.F90 ggamt.F90 ggair.F90 genveffig.F90 gencfun.F90 \
 genapwfr.F90 seceqnfv.F90 seceqnsv.F90 seceqnss.F90 getngkmax.F90 \
 dos.F90 rhoinit.F90 potplot.F90 writestate.F90 \
 potxc.F90 zpotcoul.F90 gndstate.F90 gndstate_check_for_convergence.F90\
 gndstate_force_relax.F90 gndstate_gencore_wf_density.F90 \
 gndstate_output_state.F90 \
 gndstate_output_timing_information.F90 gndstate_solvepotential.F90 \
 gndstate_write_latticevectors.F90 rhovalk.F90 readstate.F90 bandstr.F90 \
 writeeval.F90 rhoplot.F90 gwf2val.F90 gwf2cr.F90 \
 rfarray.F90 fermisurf.F90 \
 symrvf.F90 symrfmt.F90 wfplot.F90 hmlrad.F90 olprad.F90 olpistl.F90 olpaa.F90 \
 olpalo.F90 olplolo.F90 hmlistl.F90 hmlaa.F90 hmlalo.F90 hmllolo.F90 \
 hamiltonandoverlapsetup.F90 init2.F90 \
 writelsj.F90 wavefcr.F90 gdft.F90 gradrf.F90 dbxcplot.F90 effmass.F90 \
 projsbf.F90 seceqnhf.F90 writelat.F90 energynn.F90 findigp0.F90 \
 writeengy.F90 writechg.F90 writeforce.F90 genjlgpr.F90 findkpt.F90 \
 putevalfv.F90 getevalfv.F90 putevalsv.F90 getevalsv.F90 putoccsv.F90 \
 getoccsv.F90 putevecfv.F90 getevecfv.F90 putevecsv.F90 getevecsv.F90 \
 delevec.F90 outfilenamestring.F90 genvmatk.F90 genwfsv.F90 \
 findsymsite.F90 geomplot.F90

SRC_phonon = \
 phonon.F90 dyntask.F90 phcell.F90 phdisp.F90 readdyn.F90 dynqtor.F90 \
 dynrtoq.F90 dynsymapp.F90 dynsym.F90 dyndiag.F90 phdos.F90 sumrule.F90 \
 writephn.F90

SRC_optics = \
 writepmat.F90 linopt.F90 linoptk.F90 moke.F90

SRC_advanced = \
 vnlrho.F90 vnlrhomt.F90 genwiq2.F90 \
 exxengy.F90 exxengyk.F90 oepmain.F90 oepresk.F90 oepvnl.F90 \
 oepvnlk.F90 oepmagmt.F90 oepmag.F90 zrhoqint.F90
SRC_MPI=\
 mpiresumeevecfiles.F90 mpisumrhoandmag.F90 mpisync_evalsv_spnchr.F90
# rdmft.F90 rdminit.F90 genvmat.F90 genkinmat.F90 rdmvnl.F90 rdmvnlk.F90 \
# rdmengy.F90 rdmengyx.F90 rdmvaryn.F90 rdmdedn.F90 rdmminn.F90 rdmknrot.F90 \
# rdmminwf.F90 rdmvarywf.F90 rdmdedc.F90 rdmdexcdc.F90 rdmevecsv.F90 \
# rdmvnlwf.F90 rdmdexcdn.F90


SRC_lib = \
 euler.F90 wigner3j.F90 gaunt.F90 gauntyry.F90 r3mm.F90 r3mtm.F90 r3mmt.F90 \
 r3mv.F90 r3mtv.F90 r3cross.F90 r3dist.F90 r3taxi.F90 r3dot.F90 r3minv.F90 \
 r3mdet.F90 r3frac.F90 i3mdet.F90 factnm.F90 factr.F90 hermite.F90 \
 brzint.F90 sphcrd.F90 sphcover.F90 erf.F90 clebgor.F90 \
 sbessel.F90 sbesseldm.F90 genylm.F90 genrlm.F90 radmesh.F90 zmatinp.F90 \
 lopzflm.F90 sortidx.F90 gcd.F90 zfmtinp.F90 rfmtinp.F90 findband.F90 \
 gradzfmt.F90 gradrfmt.F90 ztorflm.F90 rtozflm.F90 zflmconj.F90 rotzflm.F90 \
 polynom.F90 sdelta.F90 stheta.F90 sdelta_mp.F90 stheta_mp.F90 sdelta_fd.F90 \
 stheta_fd.F90 rdiracint.F90 rdiracdme.F90 rdirac.F90 rschrodint.F90 \
 rschroddme.F90 rschrodapp.F90 reciplat.F90 \
 connect.F90 \
 flushifc.F90 spline.F90 writewiq2.F90 rfinterp.F90 rfmtctof.F90 fderiv.F90 \
 mixer.F90 fsmooth.F90 rotaxang.F90 i3minv.F90

SRC_lib_flux =

SRC_xc = \
  xc_pzca.F90 xc_pwca.F90 xc_pbe.F90 xc_am05.F90 xc_xalpha.F90 xc_wc06.F90 \
  x_wc06.F90 x_pbe.F90 c_pbe.F90 c_pbe_gcor.F90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines) $(SRC_routines_flux) $(SRC_MPI)\
 $(SRC_lib) $(SRC_lib_flux) $(SRC_advanced) $(SRC_phonon) $(SRC_optics) \
 $(SRC_xc)

OBJ = $(SRC:.F90=.o)
EXE = exciting
ITERLIB=-L./src_iterative_solver -literative_solver
exciting:	$(OBJ) ./src_iterative_solver/libiterative_solver.a
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(ITERLIB) $(LIB_LPK) $(LIB_FFT) $(LIB_SYS)

blas:
	cd BLAS; $(MAKE) # ; cp blas.a ..

lapack:
	cd LAPACK; $(MAKE) # ; cp lapack.a ..

fft:
	cd fftlib; $(MAKE) # ; cp fftlib.a ..

xc:
	cd xclib; $(MAKE) # ; cp xclib.a ..

all:	blas lapack fft exciting

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out *.aux *.dvi *.log *.pdf \
         *.tex *.toc $(EXE)

cleanall:
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	$(MAKE) clean

doc:
	./protex -s $(SRC_main) $(SRC_modules) $(SRC_routines) \
         $(SRC_routines_flux) $(SRC_lib) $(SRC_lib_flux) $(SRC_xc) \
         $(SRC_advanced) $(SRC_phonon) > exciting.tex
	 latex exciting;bibtex exciting;latex exciting;latex exciting
	 dvipdf exciting.dvi

doctemp:
	./protex -s $(SRC_main) $(SRC_routines) $(SRC_advanced)  > exciting.tex
	latex exciting

backup:
	tar -czf ex.tgz $(SRC) fftlib eos spacegroup species junk BLAS LAPACK \
         Makefile notes.txt docs docedit examples protex COPYING README \
         exciting.bib ../release

VERSION = $$(awk -F"/" '/data version/ {print $$2}' main.F90 | sed 's/ //g;s/,/./g')$
RELEASE = ../release/exciting

release:
	rm -rf $(RELEASE)
	mkdir $(RELEASE)
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	cd eos; $(MAKE) clean
	cd spacegroup; $(MAKE) clean
	cd species; $(MAKE) clean
	mkdir $(RELEASE)/src
	cp -r BLAS $(RELEASE)/src
	cp -r LAPACK $(RELEASE)/src
	cp -r fftlib $(RELEASE)/src
	cp -r eos $(RELEASE)/src
	cp -r spacegroup $(RELEASE)/src
	cp -r species $(RELEASE)/src
	rm -f $(RELEASE)/src/species/*.in
	cp $(SRC) Makefile protex exciting.bib $(RELEASE)/src
	cd examples; find . -type f -name *.OUT -exec rm -f {} \;
	cd examples; find . -type f -name *~ -exec rm -f {} \;
	cp -r examples $(RELEASE)
	mkdir $(RELEASE)/species
	cp species/*.in $(RELEASE)/species
	cp COPYING $(RELEASE)
	cp README $(RELEASE)
	cp ../Makefile $(RELEASE)
	cp ../setup $(RELEASE)
	mkdir $(RELEASE)/utilities
	cp -r NEdit $(RELEASE)/utilities
	$(MAKE) doc
	cp exciting.pdf docs/exciting
	cd spacegroup;$(MAKE) doc;cp spacegroup.pdf ../docs/spacegroup;$(MAKE) clean
	cp -r docs $(RELEASE)
	tar -C ../release -czf ../release/exciting-$(VERSION).tgz exciting
	cp docs/exciting/exciting.pdf ../release
	cp docs/spacegroup/spacegroup.pdf ../release

lines:
	cat $(SRC) | wc -l

