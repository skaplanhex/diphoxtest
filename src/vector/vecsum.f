C	
	DOUBLE PRECISION FUNCTION VEC_DSUM(U,IDIM)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM)
	TEMP = 0.D0
	DO I=1,IDIM
	  TEMP = U(I) + TEMP
	ENDDO
	VEC_DSUM = TEMP
	RETURN
	END