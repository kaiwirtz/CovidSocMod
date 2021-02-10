#############################################################################
##  Description Generic Makefile
##    targets:
##      (default)     : compilation and start
##      all           : clean, compilation, and start
##      backup        : backup of all source files (as .tgz)
##      clean         : delete all volatile files (such as *.o)
##
##    Author       Kai Wirtz (kai.wirtz@hzg.de)
#############################################################################
SHELL       = /bin/bash
PRG	        = corona
VERSION	    =
SETUP	      =
APPNAME     = $(PRG)$(VERSION)
PRGNAME     = $(PRG)
DATE    =   `date +%d%m%y`

# Current directory (under which the sissi distribution lies):SiSi2.0
SISSIPATH   = ./SiSi2.0
#SISSIPATH   = ~/tools/sisi
SRCDIR      =
CC_OBJECTS  = callSiSi.cc $(PRG).cc

C_OBJECTS   = model.c simulation.c struct_model.c sim_shell.c data.c tools.c

SOURCEFILES =    $(CC_OBJECTS) $(C_OBJECTS)

# Processors
CC            = g++
CC--          = gcc

DEFINES     = -D__C_PRE
CCINCLUDE   = -I$(SISSIPATH)/include
CCOPTS      = -O3 #  -O3 -Wall
C2OPTS      = -O3  # -Wall -g -O3
LNKOPTS     = -L$(SISSIPATH)/lib
LNKLIB      = -lSiSi2.0 -lm # --warn-common

#SOURCEFILES_1 =	$(addsuffix .cc, $(basename $(CC_OBJECTS)))
OBJECTS 	= $(addsuffix .o, $(basename $(SOURCEFILES)))

.SUFFIXES:         # Delete the default suffixes.
.SUFFIXES: .cc .c .o  # Define our suffix list.

.cc.o:
	$(CC) -c  $(C2OPTS) $(CCINCLUDE)  -o $*.o $*.cc

.c.o:
	$(CC--) -c $(CCOPTS) $(CCINCLUDE) $(DEFINES) -o $*.o $*.c

new:	$(APPNAME) run

$(APPNAME): $(OBJECTS)
	$(CC) $(LNKOPTS) -o ./$(PRGNAME).prg $(OBJECTS) $(SISSIPATH)/libsisi.a

all:  clean  $(APPNAME) run

clean:
	@echo -n "Making all clean ... "; \
	if [ -f *.tgz ]; then \
   mv *.tgz archiv; fi; \
	rm -f *.o *~ *% *.prg *.err pdist.mat *.log "#"*"#"; \
        rm -i *.tgz; \
	rm -f multi/*~ mesocosm/*~
	echo Done.

backup: clean
	@echo -n "Creating Archiv of $(APPNAME) .. "; \
	  tar zcvf $(PRGNAME)$(DATE).tgz *.*  Makefile README *.mat *.bin ; \
	  echo "Done.";
run:
	@echo  "Executing $(PRGNAME).prg $(APPNAME)$(SETUP).sim"
	@./$(PRGNAME).prg $(APPNAME)$(SETUP).sim
