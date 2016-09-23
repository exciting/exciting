.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90 .m4 .exe .xml

BUILD_TARGETS=wxml_lib wcml_lib sax_lib dom_lib

VPATH=
compile_prefix=/home/Sebastian/exciting/exciting/src/FoX/objs
install_prefix=/usr/local
LIB_DIR=$(compile_prefix)/lib
MOD_DIR=$(compile_prefix)/finclude

FPP=
FC=ifort
RANLIB=ranlib

FFLAGS=-g 
FPPFLAGS= -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DFC_ABORT_ARG -DFC_EOR_LF -DRESTRICTED_ASSOCIATED_BUG
LDFLAGS=

FCFLAGS_free_f90=
FPPFLAGS_free_F90=

INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=mod
MKDIR_P=/bin/mkdir -p
INSTALL=/usr/bin/install -c
OBJEXT=o
EXEEXT=
LIBEXT=a
LINK_O_FLAG=-o 

#INCFLAGS must be set by the user makefile

#Dependency rules are created by autoconf according to whether
#discrete preprocessing is necessary or not.
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

