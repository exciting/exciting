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
 dasum.f   daxpy.f  dcabs1.f  dcopy.f  ddot.f   dgemm.f  dgemv.f   dger.f \
 dnrm2.f   dscal.f  dspmv.f   dspr2.f  dspr.f   dswap.f  dsymv.f   dsyr2.f \
 dsyr2k.f  dtpmv.f  dtpsv.f   dtrmm.f  dtrmv.f  dtrsm.f  dznrm2.f  idamax.f \
 izamax.f  lsame.f  xerbla.f  zaxpy.f  zcopy.f  zdotc.f  zdotu.f   zdscal.f \
 zgemm.f   zgemv.f  zgerc.f   zgeru.f  zhemm.f  zhemv.f  zher2.f   zher2k.f \
 zherk.f   zhpmv.f  zhpr2.f   zhpr.f   zscal.f  zswap.f  ztpmv.f   ztpsv.f \
 ztrmm.f   ztrmv.f  ztrsm.f   ztrsv.f

OBJ = $(SRC:.f=.o)

blas:	$(OBJ)
	ar -rc libblas.a $(OBJ)

clean:
	rm -f *.o *.mod *~ *.a ifc* *.gcno gmon.out

