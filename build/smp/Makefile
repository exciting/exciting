include ../make.inc
export SUFFIX=smp
ifeq ($(BUILDSMP),true)
MTCOMPILE=true
endif

all:
ifeq ($(BUILDSMP),true)
	$(MAKE)  -f ../Make.common

endif
	

libs:
ifeq ($(BUILDSMP),true)
	$(MAKE) libs -f ../Make.common
endif

clean:
	$(MAKE) clean -f ../Make.common

cleanlibs:
	$(MAKE) cleanlibs  -f ../Make.common




	
