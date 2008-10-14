include ../make.inc

ifeq  ($(BUILDSMP),true)
ifeq ($(BUILDMPI),true)
BUILDMPISMP=true
endif
endif


all:
ifeq ($(BUILDMPISMP),true)
	make  -f ../Make.common
endif

libs:
ifeq ($(BUILDMPISMP),true)
	make libs -f ../Make.common
endif

clean:
	make clean -f ../Make.common

cleanlibs:
	make cleanlibs  -f ../Make.common





	
