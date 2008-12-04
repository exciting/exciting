include ../make.inc

FC = $(F77MT)
FFLAGS = $(F77_OPTS) 
LD = $(FC)
LDFLAGS = $(F77_OPTS) $(LIBS)
AR = ar
ARFLAGS = -rc

TMPFILES = *.mod
