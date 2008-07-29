include ../make.inc

all:
ifeq ($(BUILDMPI),true)
	make  -f ../Make.common
endif

libs:
ifeq ($(BUILDMPI),true)
	make libs -f ../Make.common

endif

clean:
	make clean -f ../Make.common

cleanlibs:
	make cleanlibs  -f ../Make.common





	
