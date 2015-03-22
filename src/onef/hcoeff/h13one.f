C coefficients devant p_1.p_3/(p_1.p_5 p_3.p_5)
C ----------------------------------------------------
	SUBROUTINE VEC_H13O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,VH13)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VH13(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VH13(I) = 0.D0
	ENDDO
	VC = (N**2-1.D0)
C	qi + qk --> qi + ph
      hh13 = 0
C	j0 = 1 : D + U --> D + ph
	VH13(1) = HH13 
C	j0 = 2 : D + Dp --> D + ph
	VH13(2) = HH13 
C	j0 = 3 : U + D --> U + ph
	VH13(3) = HH13 
C	j0 = 4 : U + Up --> U + ph
	VH13(4) = HH13 
C	qi + qbk --> qi + ph
      hh13 = 0
C	j0 = 5 : D + Ub --> D + ph
	VH13(5) = HH13 
C	j0 = 6 : D + Dpb --> D + ph
	VH13(6) = HH13 
C	j0 = 7 : U + Db --> U + ph
	VH13(7) = HH13 
C	j0 = 8 : U + Upb --> U + ph
	VH13(8) = HH13 
C	qi + qbi --> qk + ph
      hh13 = 0
C	j0 = 9 : D + Db --> U + ph
	VH13(9) = HH13 
C	j0 = 10 : D + Db --> Dp + ph
	VH13(10) = HH13 
C	j0 = 11 : U + Ub --> D + ph
	VH13(11) = HH13 
C	j0 = 12 : U + Ub --> Up + ph
	VH13(12) = HH13 
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      hh13 = 0
	VH13(13) = HH13 
C	j0 = 14 : qi + qbi --> qi + ph
      hh13 = 0
	VH13(14) = HH13 
C	j0 = 15 : qi + g --> qi + ph
      hh13 = 2*Vc*eqi**2*(2*CF-N)*(s14**3*s34+s14*s34**3+s12**3*s23+s12*
     #s23**3+s15**3*s35+s15*s35**3)/s14/s12/s34/s23
	VH13(15) = HH13 
C	j0 = 16 : qi + g --> g + ph
      hh13 = 0
	VH13(16) = HH13 
C	j0 = 17 : qi + qbi --> g + ph
      hh13 = 2/s13*Vc*eqi**2*N/s14/s24/s23*(s14**3*s24+s14*s24**3+s13**3
     #*s23+s13*s23**3+s15**3*s25+s15*s25**3)
	VH13(17) = HH13 
C	j0 = 18 : g + g --> qi + ph
      hh13 = 0
	VH13(18) = HH13 
	RETURN
	END

