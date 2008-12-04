include ../make.inc
export SUFFIX=mpismp
ifeq  ($(BUILDSMP),true)
ifeq ($(BUILDMPI),true)
BUILDMPISMP=true
F90=$(F90MT)
F77=$(F77MT)
endif
endif


all:
ifeq ($(BUILDMPISMP),true)
	$(MAKE)  -f ../Make.common
endif

libs:
ifeq ($(BUILDMPISMP),true)
	$(MAKE) libs -f ../Make.common
endif

clean:
	$(MAKE) clean -f ../Make.common

cleanlibs:
	$(MAKE) cleanlibs  -f ../Make.common





	
