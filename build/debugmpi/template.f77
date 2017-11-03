include ../make.inc

FC = $(F77)
FFLAGS = $(F77_DEBUGOPTS)
LD = $(FC)
LDFLAGS = $(F77_OPTS) $(LIBS)
AR = ar
ARFLAGS = -rc

TMPFILES = *.mod
