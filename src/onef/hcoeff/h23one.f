C coefficients devant p_2.p_3/(p_2.p_5 p_3.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H23O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,VH23)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VH23(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VH23(I) = 0.D0
	ENDDO
	VC = (N**2-1.D0)
C	qi + qk --> qi + ph
      hh23 = 0
C	j0 = 1 : D + U --> D + ph
	VH23(1) = HH23 
C	j0 = 2 : D + Dp --> D + ph
	VH23(2) = HH23 
C	j0 = 3 : U + D --> U + ph
	VH23(3) = HH23 
C	j0 = 4 : U + Up --> U + ph
	VH23(4) = HH23 
C	qi + qbk --> qi + ph
      hh23 = 0
C	j0 = 5 : D + Ub --> D + ph
	VH23(5) = HH23 
C	j0 = 6 : D + Dpb --> D + ph
	VH23(6) = HH23 
C	j0 = 7 : U + Db --> U + ph
	VH23(7) = HH23 
C	j0 = 8 : U + Upb --> U + ph
	VH23(8) = HH23 
C	qi + qbi --> qk + ph
      hh23 = 0
C	j0 = 9 : D + Db --> U + ph
	VH23(9) = HH23 
C	j0 = 10 : D + Db --> Dp + ph
	VH23(10) = HH23 
C	j0 = 11 : U + Ub --> D + ph
	VH23(11) = HH23 
C	j0 = 12 : U + Ub --> Up + ph
	VH23(12) = HH23 
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      hh23 = 0
	VH23(13) = HH23 
C	j0 = 14 : qi + qbi --> qi + ph
      hh23 = 2*Vc*eqi**2*(N*s35*s12*s15**2+N*s35**3*s12+N*s35*s12**3+N*s
     #35*s12*s23**2+N*s25*s13*s15**2+N*s25**3*s13+N*s25*s13*s23**2+N*s25
     #*s13**3-s15**3*s23+s15**2*s35*s12+s15**2*s25*s13-s23**3*s15+s23**2
     #*s35*s12+s23**2*s25*s13)*(-s23*s14+s12*s34+s24*s13-s24*s34)/s23/N/
     #s12/s13/s14/s34/s24
	VH23(14) = HH23 
C	j0 = 15 : qi + g --> qi + ph
      hh23 = 2/s12*Vc*eqi**2*N/s14/s34/s23*(s14**3*s34+s14*s34**3+s12**3
     #*s23+s12*s23**3+s15**3*s35+s15*s35**3)
	VH23(15) = HH23 
C	j0 = 16 : qi + g --> g + ph
      hh23 = 2*Vc*eqi**2*(2*CF-N)*(s14**3*s45+s14*s45**3+s12**3*s25+s12*
     #s25**3+s13**3*s35+s13*s35**3)/s23/s12/s13/s14
	VH23(16) = HH23 
C	j0 = 17 : qi + qbi --> g + ph
      hh23 = 2/s13*Vc*eqi**2*N/s14/s24/s23*(s14**3*s24+s14*s24**3+s13**3
     #*s23+s13*s23**3+s15**3*s25+s15*s25**3)
	VH23(17) = HH23 
C	j0 = 18 : g + g --> qi + ph
      hh23 = 0
	VH23(18) = HH23 
	RETURN
	END

