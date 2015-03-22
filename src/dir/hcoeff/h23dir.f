C coefficients devant p_2.p_3/(p_2.p_5 p_3.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H23D(S12,S13,S14,S15,S23,S24,S25,S34,
     #	S35,S45,VH23)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VH23(J0MAX)
	DO I = 1,J0MAX
	  VH23(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      hh23 = 0
	VH23(1) = HH23 
C	j0 = 2 : qi + g --> ph + ph
      t1 = eqi**2
      t2 = t1**2
      t3 = N**2
      t7 = s14**2
      t8 = s45**2
      hh23 = 4*t2*(t3-1)*s15*(t7+t8)/s13/s12/s23
	VH23(2) = HH23 
	RETURN
	END

