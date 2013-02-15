#libraries that are part of the distribution

INTSIMPLLIB=libmsec.a leblaiklib.a liblbfgs.a
INTLIBPATH=./lib
INTLIBS=-lxc -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys  -lFoX_common
INTERNALLIBS=  $(INTSIMPLLIB) -L$(INTLIBPATH) $(INTLIBS)
