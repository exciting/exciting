
MAKE = make
AR = ar

include ../make.inc

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC_modules = modmain.f90 modxcifc.f90

SRC_main = main.f90

SRC_routines = \
 autoradmt.f90 findprim.f90 findsym.f90 findsymcrys.f90 genppts.f90 rfinp.f90 \
 rvfcross.f90 seceqn.f90 symrf.f90 symrfir.f90 wavefmt.f90 force.f90 \
 elfplot.f90

SRC_routines_flux = \
 rhonorm.f90 energy.f90 gensdmat.f90 \
 zpotclmt.f90 writegeom.f90 nfftifc.f90 zfftifc.f90 allatoms.f90 gridsize.f90 \
 poteff.f90 genrmesh.f90 readfermi.f90 potcoul.f90 gensfacgp.f90 \
 checkmt.f90 zfinp.f90 match.f90 forcek.f90 writeefg.f90 packeff.f90 \
 findsymlat.f90 genffacg.f90 mixerifc.f90 \
 genlofr.f90 atom.f90 writefermi.f90 writekpts.f90 fsmfield.f90 mossbauer.f90 \
 occupy.f90 writelinen.f90 writeinfo.f90 readinput.f90 readspecies.f90 \
 charge.f90 moment.f90 writesym.f90 genidxlo.f90 gencore.f90 addrhocr.f90 \
 gengvec.f90 genshtmat.f90 plot1d.f90 plot2d.f90 plot3d.f90 updatpos.f90 \
 writeiad.f90 symvect.f90 vecplot.f90 genylmg.f90 linengy.f90 init0.f90 \
 init1.f90 gengpvec.f90 genpmat.f90 ggamt.f90 ggair.f90 genveffig.f90 \
 gencfun.f90 genapwfr.f90 seceqnfv.f90 seceqnsv.f90 seceqnss.f90 getngkmax.f90 \
 dos.f90 rhoinit.f90 potplot.f90 writestate.f90 potxc.f90 zpotcoul.f90 \
 gndstate.f90 rhovalk.f90 readstate.f90 bandstr.f90 writeeval.f90 rhoplot.f90 \
 gwf2val.f90 gwf2cr.f90 rfarray.f90 fermisurf.f90 symrvfir.f90 symrvf.f90 \
 symrfmt.f90 wfplot.f90 hmlrad.f90 olprad.f90 olpistl.f90 olpaa.f90 olpalo.f90 \
 olplolo.f90 hmlistl.f90 hmlaa.f90 hmlalo.f90 hmllolo.f90 init2.f90 \
 writelsj.f90 wavefcr.f90 gradrf.f90 dbxcplot.f90 effmass.f90 projsbf.f90 \
 seceqnhf.f90 writelat.f90 energynn.f90 findigp0.f90 writeengy.f90 \
 writechg.f90 writeforce.f90 genjlgpr.f90 findkpt.f90 putevalfv.f90 \
 getevalfv.f90 putevalsv.f90 getevalsv.f90 putoccsv.f90 getoccsv.f90 \
 putevecfv.f90 getevecfv.f90 putevecsv.f90 getevecsv.f90 delevec.f90 \
 genvmatk.f90 genvmat.f90 genwfsv.f90 gendmat.f90 findsymsite.f90 geomplot.f90 \
 hartfock.f90 genkinmat.f90 energykncr.f90 zfint.f90 rfpack.f90 rfint.f90 \
 writefsm.f90 genbmatk.f90 seceqnit.f90 gdft.f90 genexpiqr.f90 writeexpiqr.f90 \
 elnes.f90 symdmat.f90 genlmirep.f90

SRC_phonon = \
 phonon.f90 dyntask.f90 phcell.f90 phdisp.f90 readdyn.f90 dynqtor.f90 \
 dynrtoq.f90 dynsymapp.f90 dynsym.f90 dyndiag.f90 phdos.f90 sumrule.f90 \
 writephn.f90 phveff.f90 phdveff.f90 writedveff.f90 readdveff.f90 \
 phfext.f90 epcouple.f90 genepmat.f90 writegamma.f90 readgamma.f90 \
 phlwidth.f90 writelambda.f90 alpha2f.f90 phdelete.f90

SRC_optics = \
 writepmat.f90 moke.f90 dielectric.f90

SRC_advanced = \
 vnlrho.f90 vnlrhomt.f90 genwiq2.f90 genjlgq0r.f90 \
 exxengy.f90 exxengyk.f90 oepmain.f90 oepresk.f90 oepvnl.f90 \
 oepvnlk.f90 oepmagmt.f90 oepmag.f90 zrhogp.f90

SRC_LDAU = \
 gendmatlu.f90 genveelu.f90 genvmatlu.f90 writeldapu.f90

SRC_rdmft = \
 rdmft.f90 rdmdkdc.f90 rdmminc.f90 rdmvaryc.f90 rdmdedc.f90 rdmdexcdc.f90 \
 rdmenergy.f90 rdmengyxc.f90 rdmwritededn.f90 rdmwriteengy.f90 rdmminn.f90 \
 rdmvaryn.f90 rdmdedn.f90 rdmdexcdn.f90 rdmvnln.f90 rdmvnlc.f90 rdmdsdn.f90 \
 rdmentropy.f90

SRC_lib = \
 euler.f90 wigner3j.f90 gaunt.f90 gauntyry.f90 r3mm.f90 r3mtm.f90 r3mmt.f90 \
 r3mv.f90 r3mtv.f90 r3cross.f90 r3dist.f90 r3taxi.f90 r3dot.f90 r3minv.f90 \
 r3mdet.f90 r3frac.f90 i3mdet.f90 i3mtv.f90 factnm.f90 factr.f90 hermite.f90 \
 brzint.f90 sphcrd.f90 sphcover.f90 erf.f90 clebgor.f90 \
 sbessel.f90 sbesseldm.f90 genylm.f90 genrlm.f90 zmatinp.f90 \
 lopzflm.f90 sortidx.f90 gcd.f90 zfmtinp.f90 rfmtinp.f90 findband.f90 \
 gradzfmt.f90 gradrfmt.f90 ztorflm.f90 rtozflm.f90 zflmconj.f90 rotzflm.f90 \
 polynom.f90 sdelta.f90 stheta.f90 sdelta_mp.f90 stheta_mp.f90 sdelta_fd.f90 \
 stheta_fd.f90 sdelta_sq.f90 stheta_sq.f90 rdiracint.f90 rdiracdme.f90 \
 rdirac.f90 rschrodint.f90 rschroddme.f90 rschrodapp.f90 reciplat.f90 \
 connect.f90 flushifc.f90 spline.f90 writewiq2.f90 rfinterp.f90 rfmtctof.f90 \
 fderiv.f90 fsmooth.f90 rotaxang.f90 i3minv.f90 axangsu2.f90 z2mm.f90 \
 z2mctm.f90 z2mmct.f90 potnucl.f90 timesec.f90 rfirvec.f90 vecfbz.f90 \
 mixadapt.f90 mixpulay.f90


SRC_lib_flux =

SRC_xc = \
  xc_pzca.f90 xc_pwca.f90 xc_pbe.f90 xc_am05.f90 xc_xalpha.f90 xc_wc06.f90 \
  x_wc06.f90 x_pbe.f90 c_pbe.f90 c_pbe_gcor.f90 xc_vbh.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines) $(SRC_routines_flux) \
 $(SRC_lib) $(SRC_lib_flux) $(SRC_advanced) $(SRC_phonon) $(SRC_optics) \
 $(SRC_xc) $(SRC_LDAU) $(SRC_rdmft)

OBJ = $(SRC:.f90=.o)
EXE = exciting

exciting:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_LPK) $(LIB_FFT) $(LIB_SYS)

blas:
	cd BLAS; $(MAKE); cp blas.a ..

lapack:
	cd LAPACK; $(MAKE); cp lapack.a ..

fft:
	cd fftlib; $(MAKE); cp fftlib.a ..

xc:
	cd xclib; $(MAKE); cp xclib.a ..

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
         exciting.bib depend ../release

VERSION = $$(awk -F"/" '/data version/ {print $$2}' modmain.f90 | sed 's/ //g;s/,/./g')$
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
	cd examples; find . -type f -name *~ -exec rm -f {} \;
	cp -r examples $(RELEASE)
	mkdir $(RELEASE)/species
	cp species/*.in $(RELEASE)/species
	cp COPYING $(RELEASE)
	cp README $(RELEASE)
	cp ../Makefile $(RELEASE)
	cp ../setup $(RELEASE)
	mkdir $(RELEASE)/build/utilities
	cp -r NEdit $(RELEASE)/build/utilities
	$(MAKE) doc
	cp exciting.pdf docs/exciting
	cd spacegroup;$(MAKE) doc;cp spacegroup.pdf ../docs/spacegroup;$(MAKE) clean
	cp -r docs $(RELEASE)
	tar -C ../release -czf ../release/exciting-$(VERSION).tgz exciting
	cp docs/exciting/exciting.pdf ../release
	cp docs/spacegroup/spacegroup.pdf ../release

lines:
	cat $(SRC) | wc -l

