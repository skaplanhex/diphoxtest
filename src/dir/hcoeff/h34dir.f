C coefficients devant p_3.p_4/(p_3.p_5 p_4.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H34D(S12,S13,S14,S15,S23,S24,S25,S34,
     #	S35,S45,VH34)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VH34(J0MAX)
	DO I = 1,J0MAX
	  VH34(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      hh34 = 0
	VH34(1) = HH34 
C	j0 = 2 : qi + g --> ph + ph
      t1 = eqi**2
      t2 = t1**2
      t3 = N**2
      t7 = s12**2
      t8 = s25**2
      hh34 = 4*t2*(t3-1)*s15*(t7+t8)/s14/s13/s34
	VH34(2) = HH34 
	RETURN
	END
