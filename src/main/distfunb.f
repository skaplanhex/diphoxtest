	SUBROUTINE FSTRU(X,M2,IH,F)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION F(-6:6),DXPDF(-6:6)
        CHARACTER*20 PARM(20)
        COMMON/PDF/VALUE(20),VALUEP(20)
        COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX
C 1=D 2=U 3=S 4=C 5=B 6=T 0=G
C -1=DB -2=UB -3=SB -4=CB -5=BB -6=TB
	  if (x.ge.1.d0.or.M2.gt.q2max) then
	    f(0)=0.d0
	    do i=1,6
	      f(i)=0.d0
	      f(-i)=0.d0
	    enddo
	    return
	  endif
	IF (IH.EQ.0) THEN
c ih = 0 fonctions de distribution dans le proton
	  PARM(1) = 'NPTYPE'
	  PARM(2) = 'NGROUP'
	  PARM(3) = 'NSET'
	  CALL PDFSET(PARM,VALUE)
	  IF (M2.LT.Q2MIN) THEN
	    M2 = Q2MIN
	  ENDIF
	  M = DSQRT(M2)
	  CALL PFTOPDG(X,M,DXPDF)
	  F(0) = DXPDF(0)/X
	  DO I=1,6
	    F(I) = DXPDF(I)/X
	    F(-I) = DXPDF(-I)/X
	  ENDDO
	ELSE IF (IH.EQ.1) THEN
c ih = 1 fonctions de distribution dans l antiproton
	  PARM(1) = 'NPTYPE'
	  PARM(2) = 'NGROUP'
	  PARM(3) = 'NSET'
	  CALL PDFSET(PARM,VALUE)
	  IF (M2.LT.Q2MIN) THEN
	    M2 = Q2MIN
	  ENDIF
	  M = DSQRT(M2)
	  CALL PFTOPDG(X,M,DXPDF)
	  F(0) = DXPDF(0)/X
	  DO I=1,6
	    F(I) = DXPDF(-I)/X
	    F(-I) = DXPDF(I)/X
	  ENDDO
	ELSE IF (IH.EQ.3) THEN
c ih = 3 fonctions de distribution dans le pion
	  PARM(1) = 'NPTYPE'
	  PARM(2) = 'NGROUP'
	  PARM(3) = 'NSET'
	  CALL PDFSET(PARM,VALUEP)
	  IF (M2.LT.Q2MIN) THEN
	    M2 = Q2MIN
	  ENDIF
	  M = DSQRT(M2)
	  CALL PFTOPDG(X,M,DXPDF)
	  F(0) = DXPDF(0)/X
	  DO I=1,6
	    F(I) = DXPDF(-I)/X
	    F(-I) = DXPDF(I)/X
	  ENDDO
	ELSE IF (IH.EQ.4) THEN
c ih = 4 fonctions de distribution dans le Be
	  PARM(1) = 'NPTYPE'
	  PARM(2) = 'NGROUP'
	  PARM(3) = 'NSET'
	  CALL PDFSET(PARM,VALUEP)
	  IF (M2.LT.Q2MIN) THEN
	    M2 = Q2MIN
	  ENDIF
	  M = DSQRT(M2)
	  CALL PFTOPDG(X,M,DXPDF)
	  F(0) = DXPDF(0)/X
	  F(1) = (4.D0*DXPDF(1)+5.D0*DXPDF(2))/(9.D0*X)
	  F(2) = (4.D0*DXPDF(2)+5.D0*DXPDF(1))/(9.D0*X)
	  F(-1) = (4.D0*DXPDF(-1)+5.D0*DXPDF(-2))/(9.D0*X)
	  F(-2) = (4.D0*DXPDF(-2)+5.D0*DXPDF(-1))/(9.D0*X)
	  DO I=3,6
	    F(I) = DXPDF(-I)/X
	    F(-I) = DXPDF(I)/X
	  ENDDO
	ENDIF
	RETURN
	END
