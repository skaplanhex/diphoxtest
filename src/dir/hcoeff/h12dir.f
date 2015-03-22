C coefficients devant p_1.p_2/(p_1.p_5 p_2.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H12D(S12,S13,S14,S15,S23,S24,S25,S34,
     #	S35,S45,VH12)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VH12(J0MAX)
	DO I = 1,J0MAX
	  VH12(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t1 = eqi**2
      t2 = t1**2
      t3 = N**2
      t6 = s13**2
      t9 = s23**2
      t12 = s14**2
      t15 = s24**2
      hh12 = 4*t2*(t3-1)*(t6*s13*s23+t9*s23*s13+t12*s14*s24+t15*s24*s14)
     #/s13/s14/s23/s24
	VH12(1) = HH12 
C	j0 = 2 : qi + g --> ph + ph
      hh12 = 0
	VH12(2) = HH12 
	RETURN
	END
