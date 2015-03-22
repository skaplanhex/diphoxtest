#   
#  Makefile for hadro-production of two large-pT photons (+X)
#                               or one large-pT hadron, one large-pT photon (+X)
#                               or two large-pT hadrons (+X)
# if hadron fragmentation function are chosen in parameter.indat, 
# only the fragmentation parts ("ephad") are compiled
# the proper compilation option is called automatically by "start.pl" 
# 
#  GNUmake (use gmake)
#
#choice of the makefile to create the bases library
# MAKEBASES= Makefile.dec
# MAKEBASES= Makefile.linux
#
# choixpdf determines which set of structure functions for the proton is loaded
# CHOIXPDF = cteq5 : cteq5
# CHOIXPDF = cteq6 : cteq6
# CHOIXPDF = mrs99 : mrs99
# CHOIXPDF = mrst01 : mrst01
# CHOIXPDF = mrst02 : mrst02
# CHOIXPDF = mrst03 : mrst03
# CHOIXPDF = mrst04 : mrst04
# CHOIXPDF = pdflib : pdflib
#
CHOIXPDF	= cteq6
# CHOIXHISTO = ntuple : output as ntuple 
# CHOIXHISTO = histo :  output filled directly in histograms
#
CHOIXHISTO	= histo

# FC	= f77
# FC         = pgf77
FC         = gfortran

FFLAGS	   = -O -m32
#FFLAGS	   = -g

ifeq ($(FC),pgf77)
MAKEBASES= Makefile.pgf_linux
LDFLAGS    = -g77libs
endif
ifeq ($(FC),gfortran)
MAKEBASES= Makefile.linux
LDFLAGS    = 
endif
ifeq ($(FC),f77)
MAKEBASES= Makefile.dec
LDFLAGS    = 
endif
LINKER	   = $(FC)

GRACELDIR     = ../basesv5.1
BASESLIB      = bases

PATHCERNLIB	= /cern/pro/bin

ifeq ($(CHOIXPDF),cteq5)
# avec cteq5
  CERNLIB       = `$(PATHCERNLIB)/cernlib packlib mathlib kernlib` 
endif
ifeq ($(CHOIXPDF),cteq6)
# avec cteq6
  CERNLIB       = `$(PATHCERNLIB)/cernlib packlib mathlib kernlib` 
endif
ifeq ($(CHOIXPDF),mrs99)
# avec mrs99
  CERNLIB       = `$(PATHCERNLIB)/cernlib packlib mathlib kernlib` 
endif
ifeq ($(CHOIXPDF),mrst01)
# avec mrst01
  CERNLIB       = `$(PATHCERNLIB)/cernlib packlib mathlib kernlib` 
endif
ifeq ($(CHOIXPDF),mrst02)
# avec mrst02
  CERNLIB       = `$(PATHCERNLIB)/cernlib packlib mathlib kernlib` 
endif
ifeq ($(CHOIXPDF),mrst03)
# avec mrst03
  CERNLIB       = `$(PATHCERNLIB)/cernlib packlib mathlib kernlib` 
endif
ifeq ($(CHOIXPDF),mrst04)
# avec mrst04
  CERNLIB       = `$(PATHCERNLIB)/cernlib packlib mathlib kernlib` 
endif
ifeq ($(CHOIXPDF),pdflib)
# avec pdflib
  CERNLIB       = `$(PATHCERNLIB)/cernlib pdflib804 packlib mathlib kernlib`
endif
#

ifeq ($(CHOIXPDF),cteq5)
# avec cteq5
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs:../pdfs/cteq5
endif
ifeq ($(CHOIXPDF),cteq6)
# avec cteq6
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs:../pdfs/cteq6
endif
ifeq ($(CHOIXPDF),mrs99)
# avec mrs99
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs:../pdfs/mrs99
endif
ifeq ($(CHOIXPDF),mrst01)
# with mrst01
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs:../pdfs/mrst01
endif
ifeq ($(CHOIXPDF),mrst02)
# with mrst02
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs:../pdfs/mrst02
endif
ifeq ($(CHOIXPDF),mrst03)
# with mrst01
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs:../pdfs/mrst03
endif
ifeq ($(CHOIXPDF),mrst04)
# with mrst01
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs:../pdfs/mrst04
endif
ifeq ($(CHOIXPDF),pdflib)
# avec pdflib
  VPATH = ../src/dir/hcoeff:../src/dir/src: \
  	../src/onef/hcoeff:../src/onef/src: \
	../src/twof/hcoeff:../src/twof/src: \
	../src/histo:../src/histo/perlmod:../src/ntuple: \
	../src/main:../src/miscel:../src/vector: \
	../frag/photon:../frag/hadron: \
	../pdfs
endif

%.o : %.f
	$(FC) -m32 $(FFLAGS) -ffixed-line-length-132 $< -c -o $@

%.o : %.c
	gcc -c -o $@ -m32 $<

MDIR = htermdir.o crossdir.o twophotdn.o dir.o \
	h12dir.o h13dir.o h23dir.o h14dir.o \
	h24dir.o h34dir.o consdir.o

MONE = htermone.o crossone.o twophoton.o onefrag.o \
	h12one.o h13one.o h23one.o h14one.o \
	h24one.o h34one.o consone.o

MTWO = htermtwo.o crosstwo.o twophottn.o twofrag.o \
	h12two.o h13two.o h23two.o h14two.o \
	h24two.o h34two.o constwo.o

ifeq ($(CHOIXPDF),cteq5)
# avec cteq5
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o pftopdg_cteq5.o Ctq5Pdf.o \
	modified_trig.o
endif
ifeq ($(CHOIXPDF),cteq6)
# avec cteq6
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o pftopdg_cteq6.o Ctq6Pdf.o \
	modified_trig.o
endif
ifeq ($(CHOIXPDF),mrs99)
# avec mrs99
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o pftopdg_mrs99.o mrs99.o \
	modified_trig.o
endif
ifeq ($(CHOIXPDF),mrst01)
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o pftopdg_mrst01.o mrst2001.o \
	modified_trig.o
endif
ifeq ($(CHOIXPDF),mrst02)
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o pftopdg_mrst02.o mrst2002.o \
	modified_trig.o
endif
ifeq ($(CHOIXPDF),mrst03)
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o pftopdg_mrst03.o mrst2003c.o \
	modified_trig.o
endif
ifeq ($(CHOIXPDF),mrst04)
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o pftopdg_mrst04.o mrst2004.o \
	modified_trig.o
endif
ifeq ($(CHOIXPDF),pdflib)
# avec pdflib
  MC = alfab.o distfunb.o param.o secteff.o alfa_qed.o \
	choiscale.o charnum.o modified_trig.o
endif
ifeq ($(CHOIXHISTO),ntuple)
# creation d'un ntuple
  MAIN = twophot.o ntuple.o isolement.o doubletosingled.o \
  	doubletosingleo.o doubletosinglet.o 
endif
ifeq ($(CHOIXHISTO),histo)
# resultats directement dans les histo
  MAIN = twophot_histo.o init_paw.o end_paw.o histo_selection.o \
          further_param.o remplissage.o isolement.o \
	  doubletosingled_histo.o \
  	  doubletosingleo_histo.o doubletosinglet_histo.o
endif

FRAG = fonfra_inter2.o distribPert.o distribNonPert_setII.o \
	distribNonPert_setI.o locate.o polin2.o	polint.o fragfun_all.o \
	bkk.o kkp.o pkhff.o grille_had_charged.o

MCV = vecsum.o vecaddmulcons.o vecaddvec.o vecinit.o \
	vecmuladd.o vecmulcons.o vecmulsub.o \
	vecmulvec.o

twophot_blic: $(MAIN) $(MC) $(MDIR) $(MONE) $(MTWO) $(FRAG) $(MCV)
	$(LINKER) $(LDFLAGS) -o $@.exe -m32 $^ \
	-L$(GRACELDIR) -l$(BASESLIB) \
	$(CERNLIB)

phothad_blic: $(MAIN) $(MC) twophotdn_dum.o $(MONE) $(MTWO) $(FRAG) $(MCV)
	$(LINKER) $(LDFLAGS) -o $@.exe -m32 $^ \
	-L$(GRACELDIR) -l$(BASESLIB) \
	$(CERNLIB)

twohad_blic: $(MAIN) $(MC) twophotdn_dum.o twophoton_dum.o $(MTWO) \
	$(FRAG) $(MCV)
	$(LINKER) $(LDFLAGS) -o $@.exe -m32 $^ \
	-L$(GRACELDIR) -l$(BASESLIB) \
	$(CERNLIB)

bases:
	cd ../basesv5.1; gmake -f $(MAKEBASES);

clean:	
	rm $(MAIN) $(MC) $(MDIR) $(MONE) $(MTWO) $(FRAG) $(MCV) \
	twophotdn_dum.o twophoton_dum.o

