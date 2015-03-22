C============================================================================= 
      PROGRAM DEMO
C============================================================================= 
C      Trivial test program for CTEQ PDFs: Fred Olness  6/20/96
C
C      Revised 12 March 1999 by Randall J. Scalise (scalise@phys.psu.edu)
C      for CTEQ5
C============================================================================= 
      Implicit Double Precision (A-H, O-Z)
      Dimension PDF5 (-5:5)
      DATA ISET, X, Q /1,0.1,10./


1     WRITE(6,*) ' ENTER: ISET, X, Q '
      READ (5,*)          ISET, X, Q 

      Call SetCtq5(Iset)

      DO 11 Iparton=-5,5
      PDF5(Iparton) =  Ctq5Pdf (Iparton, X, Q)
11    CONTINUE

      DO 21 Iparton=3,5
      If(PDF5(Iparton) .ne. PDF5(-Iparton))
     >     Write(6,*) ' Error: Sea not symmetric. iparton= ',iparton
21    CONTINUE
      

      WRITE(6,111) ISET, X, Q
111   FORMAT(' ISET=',I5,' X=',F12.8,' Q='F12.2)

      WRITE(6,112) (PDF5(I),I=-2,5)
112   FORMAT(' DBAR ',1PE12.5,/,
     >       ' UBAR ',1PE12.5,/,
     >       ' G    ',1PE12.5,/,
     >       ' U    ',1PE12.5,/,
     >       ' D    ',1PE12.5,/,
     >       ' S    ',1PE12.5,/,
     >       ' C    ',1PE12.5,/,
     >       ' B    ',1PE12.5,/,
     >       ' T Quark not allowed in this set! ')

      GOTO 1
      END

C============================================================================= 
C============================================================================= 
