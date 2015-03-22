C coefficients devant p_1.p_4/(p_1.p_5 p_4.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H14D(S12,S13,S14,S15,S23,S24,S25,S34,
     #	S35,S45,VH14)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VH14(J0MAX)
	DO I = 1,J0MAX
	  VH14(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      hh14 = 0
	VH14(1) = HH14 
C	j0 = 2 : qi + g --> ph + ph
      hh14 = 0
	VH14(2) = HH14 
	RETURN
	END

