C coefficients devant p_3.p_4/(p_3.p_5 p_4.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H34O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,VH34)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VH34(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VH34(I) = 0.D0
	ENDDO
	VC = (N**2-1.D0)
C	qi + qk --> qi + ph
      hh34 = 0
C	j0 = 1 : D + U --> D + ph
	VH34(1) = HH34 
C	j0 = 2 : D + Dp --> D + ph
	VH34(2) = HH34 
C	j0 = 3 : U + D --> U + ph
	VH34(3) = HH34 
C	j0 = 4 : U + Up --> U + ph
	VH34(4) = HH34 
C	qi + qbk --> qi + ph
      hh34 = 0
C	j0 = 5 : D + Ub --> D + ph
	VH34(5) = HH34 
C	j0 = 6 : D + Dpb --> D + ph
	VH34(6) = HH34 
C	j0 = 7 : U + Db --> U + ph
	VH34(7) = HH34 
C	j0 = 8 : U + Upb --> U + ph
	VH34(8) = HH34 
C	qi + qbi --> qk + ph
C	j0 = 9 : D + Db --> U + ph
	EQI = CQI(9)
	EQK = CQK(9)
      hh34 = 2*Vc*(s15**2+s25**2+s13**2+s23**2)*(-eqi*eqk*s15*s24*s34+eq
     #i*eqk*s25*s14*s34+eqi*eqk*s13*s45*s24-eqi*eqk*s23*s45*s14+eqk**2*s
     #35*s14*s24+eqi**2*s12*s45*s34)/s34**2/s12/s14/s24
	VH34(9) = HH34 
C	j0 = 10 : D + Db --> Dp + ph
	EQI = CQI(10)
	EQK = CQK(10)
      hh34 = 2*Vc*(s15**2+s25**2+s13**2+s23**2)*(-eqi*eqk*s15*s24*s34+eq
     #i*eqk*s25*s14*s34+eqi*eqk*s13*s45*s24-eqi*eqk*s23*s45*s14+eqk**2*s
     #35*s14*s24+eqi**2*s12*s45*s34)/s34**2/s12/s14/s24
	VH34(10) = HH34 
C	j0 = 11 : U + Ub --> D + ph
	EQI = CQI(11)
	EQK = CQK(11)
      hh34 = 2*Vc*(s15**2+s25**2+s13**2+s23**2)*(-eqi*eqk*s15*s24*s34+eq
     #i*eqk*s25*s14*s34+eqi*eqk*s13*s45*s24-eqi*eqk*s23*s45*s14+eqk**2*s
     #35*s14*s24+eqi**2*s12*s45*s34)/s34**2/s12/s14/s24
	VH34(11) = HH34 
C	j0 = 12 : U + Ub --> Up + ph
	EQI = CQI(12)
	EQK = CQK(12)
      hh34 = 2*Vc*(s15**2+s25**2+s13**2+s23**2)*(-eqi*eqk*s15*s24*s34+eq
     #i*eqk*s25*s14*s34+eqi*eqk*s13*s45*s24-eqi*eqk*s23*s45*s14+eqk**2*s
     #35*s14*s24+eqi**2*s12*s45*s34)/s34**2/s12/s14/s24
	VH34(12) = HH34 
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      hh34 = 0
	VH34(13) = HH34 
C	j0 = 14 : qi + qbi --> qi + ph
      hh34 = 2*Vc*eqi**2*(N*s35*s12*s15**2+N*s35**3*s12+N*s35*s12**3+N*s
     #35*s12*s23**2+N*s25*s13*s15**2+N*s25**3*s13+N*s25*s13*s23**2+N*s25
     #*s13**3-s15**3*s23+s15**2*s35*s12+s15**2*s25*s13-s23**3*s15+s23**2
     #*s35*s12+s23**2*s25*s13)*(s14+s24)/N/s12/s13/s14/s34/s24
	VH34(14) = HH34 
C	j0 = 15 : qi + g --> qi + ph
      hh34 = 0
	VH34(15) = HH34 
C	j0 = 16 : qi + g --> g + ph
      hh34 = -2*Vc*eqi**2*(s14**3*s45+s14*s45**3+s12**3*s25+s12*s25**3+s
     #13**3*s35+s13*s35**3)*(-N*s13+2*s23*CF-s23*N)/s34/s23/s12/s14/s13
	VH34(16) = HH34 
C	j0 = 17 : qi + qbi --> g + ph
      hh34 = 0
	VH34(17) = HH34 
C	j0 = 18 : g + g --> qi + ph
      hh34 = 0
	VH34(18) = HH34 
	RETURN
	END
