C coefficients devant p_1.p_3/(p_1.p_5 p_3.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H13D(S12,S13,S14,S15,S23,S24,S25,S34,
     #	S35,S45,VH13)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VH13(J0MAX)
	DO I = 1,J0MAX
	  VH13(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      hh13 = 0
	VH13(1) = HH13 
C	j0 = 2 : qi + g --> ph + ph
      hh13 = 0
	VH13(2) = HH13 
	RETURN
	END

