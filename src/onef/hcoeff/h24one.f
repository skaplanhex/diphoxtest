C coefficients devant p_2.p_4/(p_2.p_5 p_4.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H24O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,VH24)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VH24(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VH24(I) = 0.D0
	ENDDO
	VC = (N**2-1.D0)
C	qi + qk --> qi + ph
C	j0 = 1 : D + U --> D + ph
	EQI = CQI(1)
	EQK = CQK(1)
      hh24 = 2*Vc*(s12**2+s23**2+s15**2+s35**2)*(-eqi*eqk*s12*s34*s45+eq
     #i*eqk*s23*s14*s45+eqi*eqk*s15*s24*s34-eqi*eqk*s35*s24*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s24*s45)/s24**2/s13/s14/s34
	VH24(1) = HH24 
C	j0 = 2 : D + Dp --> D + ph
	EQI = CQI(2)
	EQK = CQK(2)
      hh24 = 2*Vc*(s12**2+s23**2+s15**2+s35**2)*(-eqi*eqk*s12*s34*s45+eq
     #i*eqk*s23*s14*s45+eqi*eqk*s15*s24*s34-eqi*eqk*s35*s24*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s24*s45)/s24**2/s13/s14/s34
	VH24(2) = HH24 
C	j0 = 3 : U + D --> U + ph
	EQI = CQI(3)
	EQK = CQK(3)
      hh24 = 2*Vc*(s12**2+s23**2+s15**2+s35**2)*(-eqi*eqk*s12*s34*s45+eq
     #i*eqk*s23*s14*s45+eqi*eqk*s15*s24*s34-eqi*eqk*s35*s24*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s24*s45)/s24**2/s13/s14/s34
	VH24(3) = HH24 
C	j0 = 4 : U + Up --> U + ph
	EQI = CQI(4)
	EQK = CQK(4)
      hh24 = 2*Vc*(s12**2+s23**2+s15**2+s35**2)*(-eqi*eqk*s12*s34*s45+eq
     #i*eqk*s23*s14*s45+eqi*eqk*s15*s24*s34-eqi*eqk*s35*s24*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s24*s45)/s24**2/s13/s14/s34
	VH24(4) = HH24 
C	qi + qbk --> qi + ph
C	j0 = 5 : D + Ub --> D + ph
	EQI = CQI(5)
	EQK = CQK(5)
      hh24 = 2*Vc*(s15**2+s35**2+s12**2+s23**2)*(-eqi*eqk*s15*s34*s24+eq
     #i*eqk*s35*s14*s24+eqi*eqk*s12*s45*s34-eqi*eqk*s23*s45*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s45*s24)/s24**2/s13/s14/s34
	VH24(5) = HH24 
C	j0 = 6 : D + Dpb --> D + ph
	EQI = CQI(6)
	EQK = CQK(6)
      hh24 = 2*Vc*(s15**2+s35**2+s12**2+s23**2)*(-eqi*eqk*s15*s34*s24+eq
     #i*eqk*s35*s14*s24+eqi*eqk*s12*s45*s34-eqi*eqk*s23*s45*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s45*s24)/s24**2/s13/s14/s34
	VH24(6) = HH24 
C	j0 = 7 : U + Db --> U + ph
	EQI = CQI(7)
	EQK = CQK(7)
      hh24 = 2*Vc*(s15**2+s35**2+s12**2+s23**2)*(-eqi*eqk*s15*s34*s24+eq
     #i*eqk*s35*s14*s24+eqi*eqk*s12*s45*s34-eqi*eqk*s23*s45*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s45*s24)/s24**2/s13/s14/s34
	VH24(7) = HH24 
C	j0 = 8 : U + Upb --> U + ph
	EQI = CQI(8)
	EQK = CQK(8)
      hh24 = 2*Vc*(s15**2+s35**2+s12**2+s23**2)*(-eqi*eqk*s15*s34*s24+eq
     #i*eqk*s35*s14*s24+eqi*eqk*s12*s45*s34-eqi*eqk*s23*s45*s14+eqk**2*s
     #25*s14*s34+eqi**2*s13*s45*s24)/s24**2/s13/s14/s34
	VH24(8) = HH24 
C	qi + qbi --> qk + ph
      hh24 = 0
C	j0 = 9 : D + Db --> U + ph
	VH24(9) = HH24 
C	j0 = 10 : D + Db --> Dp + ph
	VH24(10) = HH24 
C	j0 = 11 : U + Ub --> D + ph
	VH24(11) = HH24 
C	j0 = 12 : U + Ub --> Up + ph
	VH24(12) = HH24 
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      hh24 = -2*Vc*eqi**2*(N*s23*s15*s12**2+N*s23**3*s15+N*s23*s15**3+N*
     #s23*s15*s35**2+N*s25*s13*s12**2+N*s25**3*s13+N*s25*s13*s35**2+N*s2
     #5*s13**3-s12**3*s35+s12**2*s23*s15+s12**2*s25*s13-s35**3*s12+s35**
     #2*s23*s15+s35**2*s25*s13)*(-s34+s14)/N/s23/s13/s24/s14/s34
	VH24(13) = HH24 
C	j0 = 14 : qi + qbi --> qi + ph
      hh24 = 2*Vc*eqi**2*(N*s35*s12*s15**2+N*s35**3*s12+N*s35*s12**3+N*s
     #35*s12*s23**2+N*s25*s13*s15**2+N*s25**3*s13+N*s25*s13*s23**2+N*s25
     #*s13**3-s15**3*s23+s15**2*s35*s12+s15**2*s25*s13-s23**3*s15+s23**2
     #*s35*s12+s23**2*s25*s13)*(s14-s34)/N/s12/s13/s14/s34/s24
	VH24(14) = HH24 
C	j0 = 15 : qi + g --> qi + ph
      hh24 = 0
	VH24(15) = HH24 
C	j0 = 16 : qi + g --> g + ph
      hh24 = 2*Vc*eqi**2*(s14**3*s45+s14*s45**3+s12**3*s25+s12*s25**3+s1
     #3**3*s35+s13*s35**3)*(N*s12+2*s23*CF-s23*N)/s24/s23/s13/s14/s12
	VH24(16) = HH24 
C	j0 = 17 : qi + qbi --> g + ph
      hh24 = 0
	VH24(17) = HH24 
C	j0 = 18 : g + g --> qi + ph
      hh24 = 2*Vc*eqi**2*(s34**3*s45+s34*s45**3+s23**3*s25+s23*s25**3+s1
     #3**3*s15+s13*s15**3)*(N*s23+2*s12*CF-s12*N)/s24/s12/s13/s34/s23
	VH24(18) = HH24 
	RETURN
	END
