C coefficients devant p_2.p_4/(p_2.p_5 p_4.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H24D(S12,S13,S14,S15,S23,S24,S25,S34,
     #	S35,S45,VH24)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VH24(J0MAX)
	DO I = 1,J0MAX
	  VH24(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      hh24 = 0
	VH24(1) = HH24 
C	j0 = 2 : qi + g --> ph + ph
      t1 = eqi**2
      t2 = t1**2
      t3 = N**2
      t7 = s13**2
      t8 = s35**2
      hh24 = 4*t2*(t3-1)*s15*(t7+t8)/s14/s12/s24
	VH24(2) = HH24 
	RETURN
	END
