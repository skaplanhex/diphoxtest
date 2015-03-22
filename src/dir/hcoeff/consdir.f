C coefficients constants
C ----------------------------------------------------
	SUBROUTINE VEC_CONSD(S12,S13,S14,S15,S23,S24,S25,S34,
     #	S35,S45,VCONS)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VCONS(J0MAX)
	DO I = 1,J0MAX
	  VCONS(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t1 = eqi**2
      t2 = t1**2
      t3 = N**2
      t6 = s15**2
      t7 = s25**2
      hcons = 4*t2*(t3-1)*s12*(t6+t7)/s13/s14/s23/s24
	VCONS(1) = HCONS
C	j0 = 2 : qi + g --> ph + ph
      hcons = 0
	VCONS(2) = HCONS
	RETURN
	END
