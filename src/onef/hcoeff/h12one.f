C coefficients devant p_1.p_2/(p_1.p_5 p_2.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H12O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,VH12)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VH12(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VH12(I) = 0.D0
	ENDDO
	VC = (N**2-1.D0)
C	qi + qk --> qi + ph
      hh12 = 0
C	j0 = 1 : D + U --> D + ph
	VH12(1) = HH12 
C	j0 = 2 : D + Dp --> D + ph
	VH12(2) = HH12 
C	j0 = 3 : U + D --> U + ph
	VH12(3) = HH12 
C	j0 = 4 : U + Up --> U + ph
	VH12(4) = HH12 
C	qi + qbk --> qi + ph
      hh12 = 0
C	j0 = 5 : D + Ub --> D + ph
	VH12(5) = HH12 
C	j0 = 6 : D + Dpb --> D + ph
	VH12(6) = HH12 
C	j0 = 7 : U + Db --> U + ph
	VH12(7) = HH12 
C	j0 = 8 : U + Upb --> U + ph
	VH12(8) = HH12 
C	qi + qbi --> qk + ph
      hh12 = 0
C	j0 = 9 : D + Db --> U + ph
	VH12(9) = HH12 
C	j0 = 10 : D + Db --> Dp + ph
	VH12(10) = HH12 
C	j0 = 11 : U + Ub --> D + ph
	VH12(11) = HH12 
C	j0 = 12 : U + Ub --> Up + ph
	VH12(12) = HH12 
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      hh12 = -2*Vc*eqi**2*(N*s23*s15*s12**2+N*s23**3*s15+N*s23*s15**3+N*
     #s23*s15*s35**2+N*s25*s13*s12**2+N*s25**3*s13+N*s25*s13*s35**2+N*s2
     #5*s13**3-s12**3*s35+s12**2*s23*s15+s12**2*s25*s13-s35**3*s12+s35**
     #2*s23*s15+s35**2*s25*s13)*(s12*s34-s23*s14-s24*s13-s24*s14)/s12/N/
     #s23/s13/s24/s14/s34
	VH12(13) = HH12 
C	j0 = 14 : qi + qbi --> qi + ph
      hh12 = 0
	VH12(14) = HH12 
C	j0 = 15 : qi + g --> qi + ph
      hh12 = 2/s12*Vc*eqi**2*N/s14/s34/s23*(s14**3*s34+s14*s34**3+s12**3
     #*s23+s12*s23**3+s15**3*s35+s15*s35**3)
	VH12(15) = HH12 
C	j0 = 16 : qi + g --> g + ph
      hh12 = 0
	VH12(16) = HH12 
C	j0 = 17 : qi + qbi --> g + ph
      hh12 = 2*Vc*eqi**2*(2*CF-N)*(s14**3*s24+s14*s24**3+s13**3*s23+s13*
     #s23**3+s15**3*s25+s15*s25**3)/s14/s13/s24/s23
	VH12(17) = HH12 
C	j0 = 18 : g + g --> qi + ph
      hh12 = -2*Vc*eqi**2*(2*CF-N)*(s34**3*s45+s34*s45**3+s23**3*s25+s23
     #*s25**3+s13**3*s15+s13*s15**3)/s12/s23/s13/s34
	VH12(18) = HH12 
	RETURN
	END
