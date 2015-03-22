C coefficients devant p_1.p_4/(p_1.p_5 p_4.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H14O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,VH14)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VH14(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VH14(I) = 0.D0
	ENDDO
	VC = (N**2-1.D0)
C	qi + qk --> qi + ph
      hh14 = 0
C	j0 = 1 : D + U --> D + ph
	VH14(1) = HH14 
C	j0 = 2 : D + Dp --> D + ph
	VH14(2) = HH14 
C	j0 = 3 : U + D --> U + ph
	VH14(3) = HH14 
C	j0 = 4 : U + Up --> U + ph
	VH14(4) = HH14 
C	qi + qbk --> qi + ph
      hh14 = 0
C	j0 = 5 : D + Ub --> D + ph
	VH14(5) = HH14 
C	j0 = 6 : D + Dpb --> D + ph
	VH14(6) = HH14 
C	j0 = 7 : U + Db --> U + ph
	VH14(7) = HH14 
C	j0 = 8 : U + Upb --> U + ph
	VH14(8) = HH14 
C	qi + qbi --> qk + ph
      hh14 = 0
C	j0 = 9 : D + Db --> U + ph
	VH14(9) = HH14 
C	j0 = 10 : D + Db --> Dp + ph
	VH14(10) = HH14 
C	j0 = 11 : U + Ub --> D + ph
	VH14(11) = HH14 
C	j0 = 12 : U + Ub --> Up + ph
	VH14(12) = HH14 
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      hh14 = -2*Vc*eqi**2*(N*s23*s15*s12**2+N*s23**3*s15+N*s23*s15**3+N*
     #s23*s15*s35**2+N*s25*s13*s12**2+N*s25**3*s13+N*s25*s13*s35**2+N*s2
     #5*s13**3-s12**3*s35+s12**2*s23*s15+s12**2*s25*s13-s35**3*s12+s35**
     #2*s23*s15+s35**2*s25*s13)*(-s34+s24)/N/s23/s13/s24/s14/s34
	VH14(13) = HH14 
C	j0 = 14 : qi + qbi --> qi + ph
      hh14 = 0
	VH14(14) = HH14 
C	j0 = 15 : qi + g --> qi + ph
      hh14 = 0
	VH14(15) = HH14 
C	j0 = 16 : qi + g --> g + ph
      hh14 = 0
	VH14(16) = HH14 
C	j0 = 17 : qi + qbi --> g + ph
      hh14 = 0
	VH14(17) = HH14 
C	j0 = 18 : g + g --> qi + ph
      hh14 = 2*Vc*eqi**2*(s34**3*s45+s34*s45**3+s23**3*s25+s23*s25**3+s1
     #3**3*s15+s13*s15**3)*(N*s13+2*s12*CF-s12*N)/s14/s12/s23/s34/s13
	VH14(18) = HH14 
	RETURN
	END

