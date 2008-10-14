include ../make.inc

all:
ifeq ($(BUILDSMP),true)
	make  -f ../Make.common
endif
	

libs:
ifeq ($(BUILDSMP),true)
	make libs -f ../Make.common
endif

clean:
	make clean -f ../Make.common

cleanlibs:
	make cleanlibs  -f ../Make.common




	
