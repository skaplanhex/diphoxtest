C coefficients devant p_1.p_4/(p_1.p_5 p_4.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H14T(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45
     #	,VH14)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VH14(J0MAX)
	DO I = 1,J0MAX
	  VH14(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
	HH14=4.D0*CF*(2.D0*N*CF-1.D0)*(S12**2+S34**2+S14**2+S23**2)/
     #	(2.D0*S13*S24)
	VH14(1) = HH14 
C	j0 = 2 : qi + qk --> qi + g
      hh14 = -2.D0/3.D0*CF/s13/s14*s15/s24*(s12**2+s35**2+s15**2+s23**2)
	VH14(2) = HH14 
C	j0 = 3 : qi + qbk --> qi + qbk
	HH14=8.D0*CF*(S12**2+S34**2+S14**2+S23**2)/(2.D0*S13*S24)
	VH14(3) = HH14 
C	j0 = 4 : qi + qbk --> qi + g
      hh14 = -2.D0/3.D0*CF/s13/s14*s15/s24*(s15**2+s23**2+s12**2+s35**2)
	VH14(4) = HH14 
C	j0 = 5 : qi + qbi --> qk + qbk
	HH14=8.D0*CF*(S13**2+S24**2+S14**2+S23**2)/(2.D0*S12*S34)
	VH14(5) = HH14 
C	j0 = 6 : qi + qbi --> qk + g
      hh14 = -2.D0/3.D0*CF/s12/s14*s15/s34*(s15**2+s23**2+s13**2+s25**2)
	VH14(6) = HH14 
C	j0 = 7 : qi + g --> qi + qk
      hh14 = 2.D0/3.D0*s15*CF*(2*CF*N*s14*s23+2*s34*s12-s13*s24-s14*s23)
     #*(s15**2+s34**2+s14**2+s35**2)/s14/s13/s12/s24/s23
	VH14(7) = HH14 
C	j0 = 8 : qi + g --> qi + qbk
      hh14 = 2.D0/3.D0*s15*CF*(2*s14*s23-s13*s24+2*CF*N*s34*s12-s34*s12)
     #*(s14**2+s35**2+s15**2+s34**2)/s14/s13/s12/s24/s23
	VH14(8) = HH14 
C	j0 = 9 : qi + g --> qk + qbk
      hh14 = 2.D0/3.D0*s45*CF*(2*CF*N*s13*s24+2*s14*s23-s13*s24-s34*s12)
     #*(s14**2+s35**2+s13**2+s45**2)/s14/s34/s12/s23/s24
	VH14(9) = HH14 
C	j0 = 10 : qi + qi --> qi + qi
	AH14=4.D0*CF*(2.D0*N*CF-1.D0)*(S12**2+S34**2+S14**2+S23**2)/
     #	     (2.D0*S13*S24)
	BH14=-4.D0*CF*(S12**2+S34**2+S13**2+S24**2)/(2.D0*S14*S23)
	CH14=8.D0*CF/N*((S12**2+S34**2)*
     #       (S12*S34-S13*S24-S14*S23)/(4.D0*S13*S24*S14*S23))	
	VH14(10) = AH14+BH14+CH14  
C	j0 = 11 : qi + qi --> qi + g
      ah14 = -2.D0/3.D0*CF/s13/s14*s15/s24*(s12**2+s35**2+s15**2+s23**2)
      bh14 = 2.D0/3.D0*CF*(s12**2+s35**2+s13**2+s25**2)*(2*CF*N*s13*s45*
     #s24-s13*s45*s24-s23*s14*s45+2*s12*s45*s34-s15*s24*s34+6*CF*N*s25*s
     #14*s34-3*s25*s14*s34+6*s35*s14*s24)/s14**2/s23/s34/s24
      ch14 = -2*CF*(s12*s35-s13*s25-s15*s23)*(s12**2+s35**2)*(-s34+2*s24
     #*CF*N+2*s24)/s13/s23/s14/s24/N/s34
	VH14(11) = AH14+BH14+CH14  
C	j0 = 12 : qi + qbi --> qi + qbi
	AH14=8.D0*CF*(S12**2+S34**2+S14**2+S23**2)/(2.D0*S13*S24)
	BH14=8.D0*CF*(S14**2+S23**2+S13**2+S24**2)/(2.D0*S12*S34)
	CH14=-16.D0*CF*(CF+1.D0/N)*((S14**2+S23**2)*
     #       (S14*S23-S13*S24-S12*S34)/(4.D0*S13*S24*S12*S34))	
	VH14(12) = AH14+BH14+CH14
C	j0 = 13 : qi + qbi --> qi + g
      ah14 = -2.D0/3.D0*CF/s13/s14*s15/s24*(s15**2+s23**2+s12**2+s35**2)
      bh14 = -2.D0/3.D0*CF/s12/s14*s15/s34*(s15**2+s23**2+s13**2+s25**2)
      ch14 = 0
	VH14(13) = AH14+BH14+CH14
C	j0 = 14 : qi + g --> qi + qi
      ah14 = 2.D0/3.D0*s15*CF*(2*CF*N*s14*s23+2*s34*s12-s13*s24-s14*s23)
     #*(s15**2+s34**2+s14**2+s35**2)/s14/s13/s12/s24/s23
      bh14 = 0
      ch14 = 0
	VH14(14) = AH14+BH14+CH14
C	j0 = 15 : qi + g --> qi + qbi
      ah14 = 2.D0/3.D0*s15*CF*(2*s14*s23-s13*s24+2*CF*N*s34*s12-s34*s12)
     #*(s14**2+s35**2+s15**2+s34**2)/s14/s13/s12/s24/s23
      bh14 = 2.D0/3.D0*s45*CF*(2*CF*N*s13*s24-s13*s24-s34*s12+2*s14*s23)
     #*(s14**2+s35**2+s13**2+s45**2)/s14/s34/s12/s23/s24
      ch14 = 2*CF*(-s14*s35+s13*s45+s15*s34)*(s14**2+s35**2)*(2*CF*N*s14
     #*s23+2*s14*s23-s13*s24-s34*s12+2*s12*s24*CF*N+2*s12*s24)/s14/s13/s
     #34/s12/s24/s23/N
	VH14(15) = AH14+BH14+CH14
C	j0 = 16 : qi + qbi --> g + g
      t1 = N**2
      t3 = s13**2
      t5 = s23**2
      t6 = s23*t5
      t10 = s14**2
      t11 = t10**2
      t13 = s34*s24
      t15 = s24**2
      t18 = s12*s34
      t20 = s15**2
      t21 = t20**2
      t24 = s14*s23*s24
      t26 = s15*t20
      t28 = t18*s14
      t30 = s25**2
      t42 = s14*s13*s24
      t49 = s34*s14
      t53 = t15**2
      t59 = t3**2
      t63 = s25*t30
      t68 = 3*t1*t3*t6*s14*s24-3*t11*s12*t13-3*t10*s24*t15*t18+t21*t1*t2
     #4-t26*s25*t28-t20*t30*t18*s24-t21*s12*t13+t20*t1*t30*t24+t26*t1*s2
     #5*t42+3*t1*t11*t15*s13-3*s13*t3*s12*t49*s23+3*t1*t10*t53*s13-3*t6*
     #s12*t49*s13+3*t1*t59*t24+s15*t1*t63*t42-s15*t63*t28
      hh14 = (t1-1)*t68/t10/s12/s34/s13/s23/s24/3
	VH14(16) = HH14 
C	j0 = 17 : qi + g --> qi + g
      t1 = N**2
      t3 = s12**2
      t9 = s14**2
      t10 = t9**2
      t12 = s34**2
      t16 = t12**2
      t19 = s23**2
      t20 = s23*t19
      t22 = s13*s24
      t23 = t22*s14
      t26 = s35**2
      t27 = s35*t26
      t30 = s14*s12*s34
      t32 = s15**2
      t33 = t32**2
      t35 = s24*s34
      t39 = s14*s34*s23
      t44 = s15*t32
      t64 = t3**2
      t67 = -3*s12*t3*s13*s24*s14*s23+3*t1*t10*t12*s12+3*t1*t9*t16*s12-3
     #*t20*s12*t23+s15*t1*t27*t30-t33*s13*t35+t33*t1*t39-t32*t26*t22*s34
     #+t44*t1*s35*t30-s15*t27*t23+t32*t1*t26*t39-3*t10*s13*t35-t44*s35*t
     #23-3*t9*s34*t12*t22+3*t1*t3*t20*s14*s34+3*t1*t64*t39
      hh14 = (t1-1)*t67/t9/s12/s13/s24/s23/s34/3
	VH14(17) = HH14 
C	j0 = 18 : qi + g --> g + g
      t1 = N**2
      t3 = t1*s45
      t4 = s12**2
      t5 = t4**2
      t16 = s35**2
      t18 = t1*s35*t16
      t19 = s13**2
      t23 = t3*t4
      t24 = s25**2
      t25 = t24*s24
      t28 = t1*s35
      t29 = t19**2
      t33 = s45**2
      t35 = t1*s45*t33
      t36 = s14**2
      t40 = t36**2
      t55 = s12*t4
      t66 = s25*t24
      t70 = t3*t5*s23*s14+3*t1*t5*s25*s34*s14+t3*t5*s24*s13+3*t18*t19*s2
     #4*s14+t23*t25*s13+3*t28*t29*s24*s14+t35*t36*s24*s13+t3*t40*s24*s13
     #+t3*t29*s23*s14+t3*t40*s34*s12-s15*t24*s24*s23*s12*s14-s15*t55*s24
     #*s23*s14-s15*t16*s34*s23*s14*s13+3*t1*t4*t66*s34*s14
      t75 = s24*s34
      t81 = s34*s12*s14
      t90 = s13*t19
      t122 = t35*t36*s34*s12-s15*t5*t75-s15*t29*t75+3*t18*s13*t81+t3*t29
     #*s34*s12-s15*t19*t16*s24*s34+3*t28*t90*t81-s15*t4*t25*s34+3*t1*s25
     #*t55*s24*s14*s13+t23*t24*s23*s14+3*t1*t66*s24*s12*s14*s13+t3*t19*t
     #16*s34*s12-s15*t90*s34*s23*s14+t1*t16*t19*s45*s23*s14
      hh14 = (t1-1)*(t70+t122)/s24/s34/s23/s12/t36/s13/3
	VH14(18) = HH14 
C	j0 = 19 : g + g --> qi +qbi
      t1 = N**2
      t3 = s13**2
      t5 = s14**2
      t7 = s34*s12
      t9 = s45**2
      t10 = s35**2
      t14 = s35*t10
      t16 = t7*s14
      t18 = t9**2
      t20 = s12*s13
      t23 = t5**2
      t26 = t3**2
      t30 = s45*t9
      t35 = s24**2
      t36 = t35**2
      t39 = s23*s14*s13
      t44 = s23**2
      t45 = s23*t44
      t64 = s24*s13*s14
      t69 = -3*s13*t3*t5*t7-t9*t10*t7*s13-s45*t14*t16-t18*s34*t20+3*t1*t
     #3*t23*s24+3*t1*t26*t5*s24-t30*s35*t16-3*t23*s34*t20+3*t1*t36*t39+t
     #10*t1*t9*t39+3*t1*t45*t35*s14*s13+t1*t18*t39-3*t45*s34*s12*s24*s14
     #-3*s24*t35*s34*s12*s23*s14+t14*t1*s45*t64+s35*t1*t30*t64
      hh14 = (t1-1)*t69/t5/s34/s12/s24/s23/s13/3
	VH14(19) = HH14 
C	j0 = 20 : g + g --> qi + g
      t1 = N**2
      t3 = s25**2
      t4 = s35*t3
      t6 = s24*s12
      t7 = t6*s14
      t10 = s25*t3*t1
      t12 = s23*s24
      t16 = s34*s12
      t17 = t16*s14
      t19 = s13**2
      t20 = s15**2
      t21 = s15*t20
      t24 = t1*s24*s14
      t26 = s23**2
      t27 = s23*t26
      t29 = s25*t27*t1
      t31 = s13*s24
      t32 = t31*s14
      t34 = s35*t27
      t36 = t19**2
      t39 = t26**2
      t49 = t1*s12*s14
      t56 = -3*t4*s23*t7+3*t10*s13*t12*s14+3*t10*s23*t17+t19*t21*t24+3*t
     #29*t17+3*t29*t32-3*t34*t7+s15*t36*t24+s15*t39*t24+s45*t27*t1*t31*s
     #12+s45*t26*t3*t49+s45*t3*t1*s13*s23*t6
      t57 = s34**2
      t58 = t57**2
      t63 = s45**2
      t64 = s45*t63
      t69 = s24*s34
      t100 = s45*t58*t49+s45*t39*t49+t57*t64*t49+t64*t1*s13*t69*s14+s15*
     #t27*t1*t69*s12+s15*t26*t3*t24+s15*s13*t19*t1*t17+s15*t3*t1*t12*t16
     #+t21*t1*s13*t17-3*t4*t1*t12*s12*s14-3*t34*t1*t7+s45*s34*t57*t1*t32
      t113 = s14**2
      hh14 = -(t1-1)*(t56+t100)/t1/s13/s23/s24/s34/s12/t113/3
	VH14(20) = HH14 
C	j0 = 21 : g + g --> g + g
      t1 = N**2
      t5 = s14**2
      t6 = t5**2
      t7 = s14*t6
      t10 = s12**2
      t11 = t10**2
      t15 = s34**2
      t16 = t15**2
      t22 = s13**2
      t23 = t22**2
      t27 = s25**2
      t28 = s25*t27
      t29 = s15*t28
      t30 = s34*s12
      t33 = s35**2
      t34 = s35*t33
      t35 = s15*t34
      t39 = s24**2
      t40 = t39**2
      t44 = s45*t34
      t48 = s45*t28
      t54 = s14*s12
      t57 = s45**2
      t58 = t57**2
      t60 = s13*s24
      t64 = 3*t7*s34*s12+3*s12*t11*s34*s14+3*s34*t16*s14*s12+3*t7*s13*s2
     #4+3*s13*t23*s14*s24+t29*t30*s24+t35*s34*s13*s24+3*s24*t40*s14*s13+
     #t44*s14*s13*s23+t48*s12*s13*s24+t44*t30*s13+t48*t54*s23+t58*s14*t6
     #0+t58*s34*t54
      t65 = s15**2
      t66 = t65**2
      t71 = t27**2
      t78 = s23**2
      t79 = t78**2
      t90 = t33**2
      t101 = t66*s14*t60+t66*s34*t54+3*t71*s14*t60+3*t71*s34*t54+3*t23*s
     #34*t54+3*t79*s34*t54+3*t40*s34*t54+3*t11*s14*t60+3*t79*s14*t60+3*t
     #16*s14*t60+3*t90*s34*t54+t35*s34*s14*s23+t29*s14*s24*s23+3*t90*s14
     #*t60
      hh14 = 2.D0/3.D0*(t1-1)*N*t1*(t64+t101)/t5/s34/s12/s13/s24/s23
	VH14(21) = HH14 
	RETURN
	END

