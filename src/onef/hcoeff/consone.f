C coefficients constants
C ----------------------------------------------------
	SUBROUTINE VEC_CONSO(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,VCONS)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VCONS(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VCONS(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      hcons = 0
C	j0 = 1 : D + U --> D + ph
	VCONS(1) = HCONS
C	j0 = 2 : D + Dp --> D + ph
	VCONS(2) = HCONS
C	j0 = 3 : U + D --> U + ph
	VCONS(3) = HCONS
C	j0 = 4 : U + Up --> U + ph
	VCONS(4) = HCONS
C	qi + qbk --> qi + ph
      hcons = 0
C	j0 = 5 : D + Ub --> D + ph
	VCONS(5) = HCONS
C	j0 = 6 : D + Dpb --> D + ph
	VCONS(6) = HCONS
C	j0 = 7 : U + Db --> U + ph
	VCONS(7) = HCONS
C	j0 = 8 : U + Upb --> U + ph
	VCONS(8) = HCONS
C	qi + qbi --> qk + ph
      hcons = 0
C	j0 = 9 : D + Db --> U + ph
	VCONS(9) = HCONS
C	j0 = 10 : D + Db --> Dp + ph
	VCONS(10) = HCONS
C	j0 = 11 : U + Ub --> D + ph
	VCONS(11) = HCONS
C	j0 = 12 : U + Ub --> Up + ph
	VCONS(12) = HCONS
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      hcons = 0
	VCONS(13) = HCONS
C	j0 = 14 : qi + qbi --> qi + ph
      hcons = 0
	VCONS(14) = HCONS
C	j0 = 15 : qi + g --> qi + ph
      hcons = 0
	VCONS(15) = HCONS
C	j0 = 16 : qi + g --> g + ph
      hcons = 0
	VCONS(16) = HCONS
C	j0 = 17 : qi + qbi --> g + ph
      hcons = 0
	VCONS(17) = HCONS
C	j0 = 18 : g + g --> qi + ph
      hcons = 0
	VCONS(18) = HCONS
	RETURN
	END
