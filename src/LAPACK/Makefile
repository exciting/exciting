include ../../make.inc

#-------------------------------------------------------------------------------
# Suffix Rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<

#-------------------------------------------------------------------------------
# File dependencies
#-------------------------------------------------------------------------------
SRC = \
 dgesv.f   dgetf2.f  dgetrf.f  dgetri.f  dgetrs.f  dladiv.f  dlae2.f  \
 dlaebz.f  dlaev2.f  dlagtf.f  dlagts.f  dlamch.f  dlanst.f  dlansy.f \
 dlapy2.f  dlapy3.f  dlarfb.f  dlarf.f   dlarfg.f  dlarft.f  dlarnv.f \
 dlartg.f  dlaruv.f  dlascl.f  dlaset.f  dlasr.f   dlasrt.f  dlassq.f \
 dlaswp.f  dlatrd.f  dorg2l.f  dorg2r.f  dorgql.f  dorgqr.f  dorgtr.f \
 dstebz.f  dsteqr.f  dsterf.f  dsyev.f   dsytd2.f  dsytrd.f  dtrti2.f \
 dtrtri.f  ieeeck.f  ilaenv.f  iparmq.f  xerbla.f  zgesv.f   zgetf2.f \
 zgetrf.f  zgetri.f  zgetrs.f  zheev.f   zhetd2.f  zhetrd.f  zhpevx.f \
 zhpgst.f  zhpgvx.f  zhptrd.f  zlacgv.f  zladiv.f  zlaev2.f  zlanhe.f \
 zlanhp.f  zlarfb.f  zlarf.f   zlarfg.f  zlarft.f  zlascl.f  zlaset.f \
 zlasr.f   zlassq.f  zlaswp.f  zlatrd.f  zpptrf.f  zstein.f  zsteqr.f \
 ztrti2.f  ztrtri.f  zung2l.f  zung2r.f  zungql.f  zungqr.f  zungtr.f \
 zupgtr.f  zupmtr.f

OBJ = $(SRC:.f=.o)

lapack:	$(OBJ)
	ar -rc lapack.a $(OBJ)

clean:
	rm -f *.o *.mod *~ *.a ifc* *.gcno gmon.out

