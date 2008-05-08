#
#     Coded by Diederik R. Fokkema
#
#     $Id$
# 
#
include ../../../make.inc
FC = $(F77)
FFLAGS=$(F77_OPTS)
LIB = libjdqz.a

OBJS = $(patsubst %.f,%.o,$(wildcard *.f))

$(LIB): $(LIB)($(OBJS)) 
	ranlib $@

$(LIB)(jdqz.o): jdqz.f
	$(FC) -u -O1 -c jdqz.f
	$(AR) rv $(LIB) jdqz.o
	-rm jdqz.o
