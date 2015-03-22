OBJ = mrs99.o 

#FFLAGS=  -K  +ppu -O -g
FFLAGS=  -K  +ppu +O3   
fort77 = f77

makelib  : $(OBJ)
	ar ruv  libmrslib.a *.o

clean     :  $(OBJ)
	rm *.o libmrslib.a

