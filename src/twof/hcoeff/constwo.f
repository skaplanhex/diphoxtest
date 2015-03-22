C coefficients constants
C ----------------------------------------------------
	SUBROUTINE VEC_CONST(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	VCONS)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VCONS(J0MAX)
	DO I = 1,J0MAX
	  VCONS(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
	HCONS=0.D0
	VCONS(1) = HCONS
C	j0 = 2 : qi + qk --> qi + g
      hcons = 0
	VCONS(2) = HCONS
C	j0 = 3 : qi + qbk --> qi + qbk
      hcons = 0
	VCONS(3) = HCONS
C	j0 = 4 : qi + qbk --> qi + g
      hcons = 0
	VCONS(4) = HCONS
C	j0 = 5 : qi + qbi --> qk + qbk
      hcons = 0
	VCONS(5) = HCONS
C	j0 = 6 : qi + qbi --> qk + g
      hcons = 0
	VCONS(6) = HCONS
C	j0 = 7 : qi + g --> qi + qk
      hcons = 0
	VCONS(7) = HCONS
C	j0 = 8 : qi + g --> qi + qbk
      hcons = 0
	VCONS(8) = HCONS
C	j0 = 9 : qi + g --> qk + qbk
      hcons = 0
	VCONS(9) = HCONS
C	j0 = 10 : qi + qi --> qi + qi
      hcons = 0
	VCONS(10) = HCONS
C	j0 = 11 : qi + qi --> qi + g
      hcons = 0
	VCONS(11) = HCONS
C	j0 = 12 : qi + qbi --> qi + qbi
      hcons = 0
	VCONS(12) = HCONS
C	j0 = 13 : qi + qbi --> qi + g
      hcons = 0
	VCONS(13) = HCONS
C	j0 = 14 : qi + g --> qi + qi
      hcons = 0
	VCONS(14) = HCONS
C	j0 = 15 : qi + g --> qi + qbi
      hcons = 0
	VCONS(15) = HCONS
C	j0 = 16 : qi + qbi --> g + g
      t1 = N**2
      t3 = s15**2
      t5 = s24*t1
      t7 = s25**2
      t11 = s14*t1
      t15 = s12*t3
      t17 = s12*t7
      t19 = s34*t1
      hcons = (t1-1)*(-t3*s13*t5-t7*s13*t5-t3*s23*t11-t7*s23*t11+t15*s34
     #+t17*s34+t15*t19+t17*t19)/s34/s23/s14/s13/s24/t1
	VCONS(16) = HCONS
C	j0 = 17 : qi + g --> qi + g
      t1 = N**2
      t3 = s15**2
      t5 = s34*t1
      t7 = s35**2
      t11 = s14*t1
      t15 = s13*t3
      t17 = s13*t7
      t19 = s24*t1
      hcons = -(t1-1)*(t3*s12*t5+t7*s12*t5+t3*s23*t11+t7*s23*t11-t15*s24
     #-t17*s24-t15*t19-t17*t19)/s24/s23/s14/s12/s34/t1
	VCONS(17) = HCONS
C	j0 = 18 : qi + g --> g + g
      hcons = 0
	VCONS(18) = HCONS
C	j0 = 19 : g + g --> qi +qbi
      t1 = N**2
      t3 = s35**2
      t5 = s24*t1
      t7 = s45**2
      t11 = s23*t1
      t15 = s34*t3
      t17 = s34*t7
      t19 = s12*t1
      hcons = (t1-1)*(-t3*s13*t5-t7*s13*t5-t3*s14*t11-t7*s14*t11+t15*s12
     #+t17*s12+t15*t19+t17*t19)/s12/s14/s23/s13/s24/t1
	VCONS(19) = HCONS
C	j0 = 20 : g + g --> qi + g
      hcons = 0
	VCONS(20) = HCONS
C	j0 = 21 : g + g --> g + g
      hcons = 0
	VCONS(21) = HCONS
	RETURN
	END
