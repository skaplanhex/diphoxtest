C*****************************************************
	SUBROUTINE SPART15ZO(SC,TC,UC,Z1,M,PTM,PART15Z)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART15Z(K0MAX),V115Z(J0MAX),V215Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_15ZO(SC,TC,UC,Z1,M,PTM,CQI,CQK,V115Z)
	CALL VEC_25ZO(SC,UC,TC,Z1,M,PTM,CQI,CQK,V215Z)
	DO I=1,J0MAX
	  PART15Z(I) = V115Z(I)
	  PART15Z(I+18) = V215Z(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART15LO(SC,TC,UC,Z1,PART15L)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART15L(K0MAX),V115L(J0MAX),V215L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_15LO(SC,TC,UC,Z1,CQI,CQK,V115L)
	CALL VEC_25LO(SC,UC,TC,Z1,CQI,CQK,V215L)
	DO I=1,J0MAX
	  PART15L(I) = V115L(I)
	  PART15L(I+18) = V215L(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART15DO(SC,TC,UC,M,PTM,PART15D)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART15D(K0MAX),V115D(J0MAX),V215D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_15DO(SC,TC,UC,M,PTM,CQI,CQK,V115D)
	CALL VEC_25DO(SC,UC,TC,M,PTM,CQI,CQK,V215D)
	DO I=1,J0MAX
	  PART15D(I) = V115D(I)
	  PART15D(I+18) = V215D(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART25ZO(SC,TC,UC,Z2,M,PTM,PART25Z)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART25Z(K0MAX),V125Z(J0MAX),V225Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_25ZO(SC,TC,UC,Z2,M,PTM,CQI,CQK,V125Z)
	CALL VEC_15ZO(SC,UC,TC,Z2,M,PTM,CQI,CQK,V225Z)
	DO I=1,J0MAX
	  PART25Z(I) = V125Z(I)
	  PART25Z(I+18) = V225Z(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART25LO(SC,TC,UC,Z2,PART25L)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART25L(K0MAX),V125L(J0MAX),V225L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_25LO(SC,TC,UC,Z2,CQI,CQK,V125L)
	CALL VEC_15LO(SC,UC,TC,Z2,CQI,CQK,V225L)
	DO I=1,J0MAX
	  PART25L(I) = V125L(I)
	  PART25L(I+18) = V225L(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART25DO(SC,TC,UC,M,PTM,PART25D)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART25D(K0MAX),V125D(J0MAX),V225D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_25DO(SC,TC,UC,M,PTM,CQI,CQK,V125D)
	CALL VEC_15DO(SC,UC,TC,M,PTM,CQI,CQK,V225D)
	DO I=1,J0MAX
	  PART25D(I) = V125D(I)
	  PART25D(I+18) = V225D(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART35RO(SC,TC,UC,Z3,R,PART35R)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART35R(K0MAX),V135R(J0MAX),V235R(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_35RO(SC,TC,UC,Z3,R,CQI,CQK,V135R)
	CALL VEC_35RO(SC,UC,TC,Z3,R,CQI,CQK,V235R)
	DO I=1,J0MAX
	  PART35R(I) = V135R(I)
	  PART35R(I+18) = V235R(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART35ZO(SC,TC,UC,Z3,MF,PT3,PART35Z)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART35Z(K0MAX),V135Z(J0MAX),V235Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_35ZO(SC,TC,UC,Z3,MF,PT3,CQI,CQK,V135Z)
	CALL VEC_35ZO(SC,UC,TC,Z3,MF,PT3,CQI,CQK,V235Z)
	DO I=1,J0MAX
	  PART35Z(I) = V135Z(I)
	  PART35Z(I+18) = V235Z(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART35LO(SC,TC,UC,Z3,PART35L)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART35L(K0MAX),V135L(J0MAX),V235L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_35LO(SC,TC,UC,Z3,CQI,CQK,V135L)
	CALL VEC_35LO(SC,UC,TC,Z3,CQI,CQK,V235L)
	DO I=1,J0MAX
	  PART35L(I) = V135L(I)
	  PART35L(I+18) = V235L(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART35DO(SC,TC,UC,MF,PT3,PART35D)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART35D(K0MAX),V135D(J0MAX),V235D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_35DO(SC,TC,UC,MF,PT3,CQI,CQK,V135D)
	CALL VEC_35DO(SC,UC,TC,MF,PT3,CQI,CQK,V235D)
	DO I=1,J0MAX
	  PART35D(I) = V135D(I)
	  PART35D(I+18) = V235D(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART45RO(SC,TC,UC,Z4,R,PART45R)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART45R(K0MAX),V145R(J0MAX),V245R(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_45RO(SC,TC,UC,Z4,R,CQI,CQK,V145R)
	CALL VEC_45RO(SC,UC,TC,Z4,R,CQI,CQK,V245R)
	DO I=1,J0MAX
	  PART45R(I) = V145R(I)
	  PART45R(I+18) = V245R(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART45ZO(SC,TC,UC,Z4,MF,PT4,PART45Z)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART45Z(K0MAX),V145Z(J0MAX),V245Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_45ZO(SC,TC,UC,Z4,MF,PT4,CQI,CQK,V145Z)
	CALL VEC_45ZO(SC,UC,TC,Z4,MF,PT4,CQI,CQK,V245Z)
	DO I=1,J0MAX
	  PART45Z(I) = V145Z(I)
	  PART45Z(I+18) = V245Z(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART45LO(SC,TC,UC,Z4,PART45L)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART45L(K0MAX),V145L(J0MAX),V245L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_45LO(SC,TC,UC,Z4,CQI,CQK,V145L)
	CALL VEC_45LO(SC,UC,TC,Z4,CQI,CQK,V245L)
	DO I=1,J0MAX
	  PART45L(I) = V145L(I)
	  PART45L(I+18) = V245L(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPART45DO(SC,TC,UC,MF,PT4,PART45D)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PART45D(K0MAX),V145D(J0MAX),V245D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_45DO(SC,TC,UC,MF,PT4,CQI,CQK,V145D)
	CALL VEC_45DO(SC,UC,TC,MF,PT4,CQI,CQK,V245D)
	DO I=1,J0MAX
	  PART45D(I) = V145D(I)
	  PART45D(I+18) = V245D(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPARTVIO(S,SC,TC,UC,MU,PT,PTM,PARTVI)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PARTVI(K0MAX),V1VI(J0MAX),V2VI(J0MAX)
	DIMENSION CQI(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	CALL VEC_VIO(S,SC,TC,UC,MU,PT,PTM,CQI,V1VI)
	CALL VEC_VIO(S,SC,UC,TC,MU,PT,PTM,CQI,V2VI)
	DO I=1,J0MAX
	  PARTVI(I) = V1VI(I)
	  PARTVI(I+18) = V2VI(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SPARTVI2O(S,SC,TC,UC,FI,PARTVI2)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION PARTVI2(K0MAX),V1VI(J0MAX),V2VI(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	CALL VEC_VI2O(S,SC,TC,UC,FI,CQI,CQK,V1VI)
	CALL VEC_VI2O(S,SC,UC,TC,FI,CQI,CQK,V2VI)
	DO I=1,J0MAX
	  PARTVI2(I) = V1VI(I)
	  PARTVI2(I+18) = V2VI(I) 
	ENDDO 
	RETURN
	END	
C*****************************************************
	SUBROUTINE SH12O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,H12)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION H12(K0MAX)
	DIMENSION V1H(J0MAX),V2H(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	ZERO = 0.D0
	CALL VEC_DINIT(H12,K0MAX,ZERO)	
	CALL VEC_H12O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,V1H)
	CALL VEC_H12O(S12,S23,S24,S25,S13,S14,S15,S34,S35,S45,
     #	CQI,CQK,V2H)
	DO I=1,J0MAX
	  H12(I) = V1H(I)
	  H12(I+18) = V2H(I) 
	ENDDO 
	RETURN
	END
C
	SUBROUTINE SH13O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,H13)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION H13(K0MAX)
	DIMENSION V1H(J0MAX),V2H(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	ZERO = 0.D0
	CALL VEC_DINIT(H13,K0MAX,ZERO)	
	CALL VEC_H13O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,V1H)
	CALL VEC_H23O(S12,S23,S24,S25,S13,S14,S15,S34,S35,S45,
     #	CQI,CQK,V2H)
	DO I=1,J0MAX
	  H13(I) = V1H(I)
	  H13(I+18) = V2H(I) 
	ENDDO 
	RETURN
	END
C
	SUBROUTINE SH23O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,H23)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION H23(K0MAX)
	DIMENSION V1H(J0MAX),V2H(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	ZERO = 0.D0
	CALL VEC_DINIT(H23,K0MAX,ZERO)	
	CALL VEC_H23O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,V1H)
	CALL VEC_H13O(S12,S23,S24,S25,S13,S14,S15,S34,S35,S45,
     #	CQI,CQK,V2H)
	DO I=1,J0MAX
	  H23(I) = V1H(I)
	  H23(I+18) = V2H(I) 
	ENDDO 
	RETURN
	END
C
	SUBROUTINE SH14O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,H14)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION H14(K0MAX)
	DIMENSION V1H(J0MAX),V2H(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	ZERO = 0.D0
	CALL VEC_DINIT(H14,K0MAX,ZERO)	
	CALL VEC_H14O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,V1H)
	CALL VEC_H24O(S12,S23,S24,S25,S13,S14,S15,S34,S35,S45,
     #	CQI,CQK,V2H)
	DO I=1,J0MAX
	  H14(I) = V1H(I)
	  H14(I+18) = V2H(I) 
	ENDDO 
	RETURN
	END
C
	SUBROUTINE SH24O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,H24)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION H24(K0MAX)
	DIMENSION V1H(J0MAX),V2H(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	ZERO = 0.D0
	CALL VEC_DINIT(H24,K0MAX,ZERO)	
	CALL VEC_H24O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,V1H)
	CALL VEC_H14O(S12,S23,S24,S25,S13,S14,S15,S34,S35,S45,
     #	CQI,CQK,V2H)
	DO I=1,J0MAX
	  H24(I) = V1H(I)
	  H24(I+18) = V2H(I) 
	ENDDO 
	RETURN
	END
C
	SUBROUTINE SH34O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,H34)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION H34(K0MAX)
	DIMENSION V1H(J0MAX),V2H(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	ZERO = 0.D0
	CALL VEC_DINIT(H34,K0MAX,ZERO)	
	CALL VEC_H34O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,V1H)
	CALL VEC_H34O(S12,S23,S24,S25,S13,S14,S15,S34,S35,S45,
     #	CQI,CQK,V2H)
	DO I=1,J0MAX
	  H34(I) = V1H(I)
	  H34(I+18) = V2H(I) 
	ENDDO 
	RETURN
	END
C
	SUBROUTINE SCONSO(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,CONS)
	IMPLICIT REAL*8 (A-H,L-Z)
	DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION CONS(K0MAX)
	DIMENSION V1H(J0MAX),V2H(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DATA (CQI(I),I=1,J0MAX)/2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,2*-0.3333333333333333D0,
     #	2*0.6666666666666667D0,6*0.3333333333333333D0/
	DATA (CQK(I),I=1,J0MAX)/0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,
     #  0.6666666666666667D0,2*-0.3333333333333333D0,
     #	0.6666666666666667D0,0.6666666666666667D0,
     #	2*-0.3333333333333333D0,0.6666666666666667D0,	
     #  6*0.3333333333333333D0/
	ZERO = 0.D0
	CALL VEC_DINIT(CONS,K0MAX,ZERO)	
	CALL VEC_CONSO(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,V1H)
	CALL VEC_CONSO(S12,S23,S24,S25,S13,S14,S15,S34,S35,S45,
     #	CQI,CQK,V2H)
	DO I=1,J0MAX
	  CONS(I) = V1H(I)
	  CONS(I+18) = V2H(I) 
	ENDDO 
	RETURN
	END
C*****************************************************
	SUBROUTINE BORNO(SC,TC,UC,N,RBORN)
	IMPLICIT REAL*8 (A-H,L-Z)
	PARAMETER (K0MAX=36,J0MAX=18)
	DIMENSION RBORN(K0MAX),V1B(J0MAX),V2B(J0MAX)
	CALL VEC_BO(SC,TC,UC,N,V1B)
	CALL VEC_BO(SC,UC,TC,N,V2B)
	DO I=1,J0MAX
	  RBORN(I) = V1B(I)
	  RBORN(I+18) = V2B(I) 
	ENDDO 
	RETURN
	END	
