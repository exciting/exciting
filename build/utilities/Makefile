# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 


.DEFAULT:
	-touch $@
all: a.out
SRC =
OBJ =
clean: neat
	-rm -f .cppdefs $(OBJ) a.out
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
a.out: $(OBJ) 
	$(LD) $(OBJ) -o a.out  $(LDFLAGS)
