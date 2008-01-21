#
#     Coded by Diederik R. Fokkema
#
#     $Id$
# 
#
include ../../../../../make.inc
export SHELL=/bin/sh
export FFLAGS = -u -O
FC = $(F77)
LIB = libjdqz.a

OBJS = $(patsubst %.f,%.o,$(wildcard *.f))

$(LIB): $(LIB)($(OBJS)) 
	ranlib $@

$(LIB)(jdqz.o): jdqz.f
	$(FC) -u -O1 -c jdqz.f
	$(AR) rv $(LIB) jdqz.o
	-rm jdqz.o
