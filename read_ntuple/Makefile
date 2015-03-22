#
#   Makefile for GNUmake (utiliser gmake)
#

#FC         = f77
FC         = g77

#FFLAGS	   = -g
#FFLAGS	   = -O4
FFLAGS	   = -O
#LDFLAGS    =

LINKER	   = $(FC)

CERNLIB       = `/cern/pro/bin/cernlib packlib mathlib kernlib` 
#


%.o : %.f
	$(FC) -c -o $@ $(FFLAGS) $<

MC = lecturn.o  init_paw.o end_paw.o further_param.o remplissage.o


lecture: $(MC) 
	$(LINKER) $(LDFLAGS) -o $@.exe $^ \
	$(CERNLIB)

clean:	
	rm $(MC)

