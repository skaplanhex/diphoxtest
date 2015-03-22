#
#   Makefile for GNUmake (utiliser gmake sur HP)
#

FC         = f77
FFLAGS	   = -O4
AR         = /bin/ar
LIBBASES   = libvec.a

%.o : %.f	
	$(FC) $(FFLAGS) $< -o $@

MC = vecsum.o vecaddmulcons.o vecaddvec.o vecinit.o \
	vecmuladd.o vecmulcons.o vecmulsub.o \
	vecmulvec.o
 
all: $(MC)
	$(AR) -rs $(LIBBASES) $^;
	mv $(LIBBASES) /users/theorie/guillet/lib/.
