C	
	SUBROUTINE VEC_DADD_VECTOR(U,V,IDIM,R)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM),V(IDIM),R(IDIM)
	DO I=1,IDIM
	  R(I) = U(I) + V(I)
	ENDDO
	RETURN
	END
C	
	SUBROUTINE VEC_DMULT_VECTOR(U,V,IDIM,R)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM),V(IDIM),R(IDIM)
	DO I=1,IDIM
	  R(I) = U(I)*V(I)
	ENDDO
	RETURN
	END
C	
	SUBROUTINE VEC_DMULT_CONSTANT(U,IDIM,A,R)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM),R(IDIM)
	DO I=1,IDIM
	  R(I) = U(I)*A
	ENDDO
	RETURN
	END
C	
	SUBROUTINE VEC_DMULT_ADD(U,V,IDIM,A,R)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM),V(IDIM),R(IDIM)
	DO I=1,IDIM
	  R(I) = A*V(I) + U(I)
	ENDDO
	RETURN
	END
C	
	SUBROUTINE VEC_DMULT_SUB(U,V,IDIM,A,R)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM),V(IDIM),R(IDIM)
	DO I=1,IDIM
	  R(I) = A*V(I) - U(I)
	ENDDO
	RETURN
	END
C	
	SUBROUTINE VEC_DADD_MULT_CONSTANT(U,V,IDIM,A,R)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM),V(IDIM),R(IDIM)
	DO I=1,IDIM
	  R(I) = ( U(I) + V(I) )*A
	ENDDO
	RETURN
	END
C	
	SUBROUTINE VEC_DINIT(U,IDIM,A)
	IMPLICIT DOUBLE PRECISION (A-H,L-Z)
	DIMENSION U(IDIM)
	DO I=1,IDIM
	  U(I) = A
	ENDDO
	RETURN
	END
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
