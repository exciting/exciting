export SUFFIX=mpi
include ../make.inc

all:
ifeq ($(BUILDMPI),true)
	$(MAKE)  -f ../Make.common
endif

libs:
ifeq ($(BUILDMPI),true)
	$(MAKE) libs -f ../Make.common
endif

clean:
	$(MAKE) clean -f ../Make.common

cleanlibs:
	$(MAKE) cleanlibs  -f ../Make.common





	
