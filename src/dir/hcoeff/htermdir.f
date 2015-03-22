C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 1
C*****************************************************
	SUBROUTINE VEC_15ZD(SC,TC,UC,Z1,M,PTM,VPART15Z)
C PARTIE EN 1/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART15Z(J0MAX)
	DO I = 1,J0MAX
	  VPART15Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = eqi**4*F(SC,TC,UC,N)*(ANM4_QQ(z1)+2*A4_QQ(z1)*dlog(ptm)-2*A4_
     #QQ(z1)*dlog(M)-FAQQ(z1))
	VPART15Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART15Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_15LD(SC,TC,UC,Z1,VPART15L)
C PARTIE EN LOG(1-Z1)/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART15L(J0MAX)
	DO I = 1,J0MAX
	  VPART15L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = -eqi**4*F(SC,TC,UC,N)*FBQQ(z1)
	VPART15L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART15L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_15DD(SC,TC,UC,M,PTM,VPART15D)
C PARTIE EN DELTA(1-Z1)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART15D(J0MAX)
	DO I = 1,J0MAX
	  VPART15D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = -eqi**4*F(SC,TC,UC,N)*(-2*B_QQ()*dlog(ptm)+2*B_QQ()*dlog(M)+F
     #CQQ(1.D0))
	VPART15D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART15D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 2
C*****************************************************
	SUBROUTINE VEC_25ZD(SC,TC,UC,Z2,M,PTM,VPART25Z)
C PARTIE EN 1/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	COMMON/AURENCHE/IAUREN
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART25Z(J0MAX)
	DO I = 1,J0MAX
	  VPART25Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = eqi**4*F(SC,TC,UC,N)*(ANM4_QQ(z2)+2*A4_QQ(z2)*dlog(ptm)-2*A4_
     #QQ(z2)*dlog(M)-FAQQ(z2))
	VPART25Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
	IF (IAUREN.EQ.0) THEN
c convention 2*(1-epsilon) degrees de polaristion pour le gluon
      t0 = -eqi**4*F(SC,TC,UC,N)*(N-1)*(N+1)*(-ANM4_QG(z2)-2*A4_QG(z2)*d
     #log(ptm)+2*A4_QG(z2)*dlog(M)+FAQG(z2))/N
	ELSE IF (IAUREN.EQ.1) THEN
c convention 2 degrees de polaristion pour le gluon
      t0 = eqi**4*F(SC,TC,UC,N)*(N-1)*(N+1)*(ANM4_QG(z2)+A4_QG(z2)+2*A4_
     #QG(z2)*dlog(ptm)-2*A4_QG(z2)*dlog(M)-FAQG(z2))/N
	ENDIF
	VPART25Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_25LD(SC,TC,UC,Z2,VPART25L)
C PARTIE EN LOG(1-Z2)/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART25L(J0MAX)
	DO I = 1,J0MAX
	  VPART25L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = -eqi**4*F(SC,TC,UC,N)*FBQQ(z2)
	VPART25L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -eqi**4*F(SC,TC,UC,N)*FBQG(z2)*(N-1)*(N+1)/N
	VPART25L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_25DD(SC,TC,UC,M,PTM,VPART25D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART25D(J0MAX)
	DO I = 1,J0MAX
	  VPART25D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = -eqi**4*F(SC,TC,UC,N)*(-2*B_QQ()*dlog(ptm)+2*B_QQ()*dlog(M)+F
     #CQQ(1.D0))
	VPART25D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -eqi**4*F(SC,TC,UC,N)*FCQG(1.D0)*(N-1)*(N+1)/N
	VPART25D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 3
C*****************************************************
	SUBROUTINE VEC_35RD(SC,TC,UC,Z3,R,VPART35R)
C PARTIE EN LOG(R)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35R(J0MAX)
	DO I = 1,J0MAX
	  VPART35R(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART35R(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -eqi**4*A4_PQ(z3)*dlog(R**2)/(1-z3)*E(TC,SC,UC,N)
	VPART35R(2) = T0
	RETURN
	END
	SUBROUTINE VEC_35ZD(SC,TC,UC,Z3,MF,PT3,VPART35Z)
C PARTIE EN 1/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35Z(J0MAX)
	DO I = 1,J0MAX
	  VPART35Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART35Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -eqi**4*E(TC,SC,UC,N)*(ANM4_PQ(z3)+2*A4_PQ(z3)*dlog(pt3)-2*A4
     #_PQ(z3)*dlog(Mf)-DAPQ(z3))
	VPART35Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_35LD(SC,TC,UC,Z3,VPART35L)
C PARTIE EN LOG(1-Z3)/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35L(J0MAX)
	DO I = 1,J0MAX
	  VPART35L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART35L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -eqi**4*E(TC,SC,UC,N)*(-DBPQ(z3)+2*A4_PQ(z3))
	VPART35L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_35DD(SC,TC,UC,MF,PT3,VPART35D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35D(J0MAX)
	DO I = 1,J0MAX
	  VPART35D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART35D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = eqi**4*E(TC,SC,UC,N)*DCPQ(1.D0)
	VPART35D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 4
C*****************************************************
	SUBROUTINE VEC_45RD(SC,TC,UC,Z4,R,VPART45R)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45R(J0MAX)
	DO I = 1,J0MAX
	  VPART45R(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART45R(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -eqi**4*A4_PQ(z4)*dlog(R**2)/(1-z4)*E(UC,SC,TC,N)
	VPART45R(2) = T0
	RETURN
	END
	SUBROUTINE VEC_45ZD(SC,TC,UC,Z4,MF,PT4,VPART45Z)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45Z(J0MAX)
	DO I = 1,J0MAX
	  VPART45Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART45Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = eqi**4*E(UC,SC,TC,N)*(-ANM4_PQ(z4)-2*A4_PQ(z4)*dlog(pt4)+2*A4
     #_PQ(z4)*dlog(Mf)+DAPQ(z4))
	VPART45Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_45LD(SC,TC,UC,Z4,VPART45L)
C PARTIE EN LOG(1-Z4)/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45L(J0MAX)
	DO I = 1,J0MAX
	  VPART45L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART45L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -eqi**4*E(UC,SC,TC,N)*(-DBPQ(z4)+2*A4_PQ(z4))
	VPART45L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_45DD(SC,TC,UC,MF,PT4,VPART45D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45D(J0MAX)
	DO I = 1,J0MAX
	  VPART45D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qbi --> ph + ph
      t0 = 0
	VPART45D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = eqi**4*E(UC,SC,TC,N)*DCPQ(1.D0)
	VPART45D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE LA PARTIE VIRTUELLE
C*****************************************************
	SUBROUTINE VEC_VID(S,SC,TC,UC,MU,PT,PTM,VPARTVI)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPARTVI(J0MAX)
	PI=4.D0*DATAN(1.D0)
	DO I = 1,J0MAX
	  VPARTVI(I) = 0.D0
	ENDDO
	MU2=MU**2
C	j0 = 1 : qi + qbi --> ph + ph
C termes provenant des parties collineaires
      teraj = -2*B_QQ()*dlog(ptm**2/S)*eqi**4*F(SC,TC,UC,N)+(N**2-1)/N*d
     #log(SC/S)*F(SC,TC,UC,N)*eqi**4*dlog(ptm**2/S)-A4_QQ(1.D0)*eqi**4*F
     #(SC,TC,UC,N)*dlog(ptm**2/S)**2/2
C termes finies d'Ellis et Sexton
      teres = 4*eqi**4*(N**2-1)*(SC**2*(-Pi**2+dlog(SC/S)**2)+(TC**2+SC*
     #*2)*dlog(-TC/S)**2+(-2*TC**2-2*SC**2)*dlog(SC/S)*dlog(-TC/S)-2*TC*
     #UC*dlog(SC/S)+UC*(3*UC+2*TC)*dlog(-TC/S)+Pi**2*SC**2+Pi**2*TC**2/2
     #+Pi**2*UC**2/2-7.D0/2.D0*TC**2-7.D0/2.D0*UC**2)/UC/TC+4*eqi**4*(N*
     #*2-1)*(SC**2*(-Pi**2+dlog(SC/S)**2)+(UC**2+SC**2)*dlog(-UC/S)**2+(
     #-2*UC**2-2*SC**2)*dlog(SC/S)*dlog(-UC/S)-2*TC*UC*dlog(SC/S)+TC*(3*
     #TC+2*UC)*dlog(-UC/S)+Pi**2*SC**2+Pi**2*TC**2/2+Pi**2*UC**2/2-7.D0/
     #2.D0*UC**2-7.D0/2.D0*TC**2)/TC/UC
	VPARTVI(1)=(TERAJ+TERES)
C	j0 = 2 : qi + g --> ph + ph
      teraj = 0
      teres = 0
	VPARTVI(2)=(TERAJ+TERES)
	RETURN
	END
	SUBROUTINE VEC_VI2D(S,SC,TC,UC,YS,FI,VPARTVI2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2)
	DIMENSION VPARTVI2(J0MAX),HH34(J0MAX)
	PI=4.D0*DATAN(1.D0)
	DO I = 1,J0MAX
	  VPARTVI2(I) = 0.D0
	  HH34(I) = 0.D0
	ENDDO
C-------------------------------------------------------------------
	YS = DLOG(TC/UC)/2.D0
	s12 = sc/2.d0
	s13 = -tc/2.d0
	s23 = -uc/2.d0
	s14 = s23
	s24 = s13
	s34 = s12
	s15 = 0.d0
	s25 = 0.d0
	s35 = 0.d0
	s45 = 0.d0
	CALL VEC_H34D(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,HH34)
	DW=DCOSH(2.D0*YS)+DCOS(2.D0*FI)
	XJ1=DLOG(DSIN(FI))
	XJ2=DSIN(2.D0*FI)*DLOG(DSIN(FI))*DATAN(DSIN(FI)/(1.D0-DCOS(FI)))
	A34=4.D0*XJ2+DSINH(2.D0*YS)*2.D0*YS*XJ1
	FFG=A34/DW+DLOG(2.D0)*DLOG(4.D0*DCOSH(YS)**2)
	DO I = 1,J0MAX
	  VPARTVI2(I) = HH34(I)*(FFG+FFG)/PI
	ENDDO
	RETURN
	END
C*****************************************************
C	
	SUBROUTINE STRFRAD(X1,IH1,X2,IH2,SFD)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/GDF/JNF
	COMMON/SCALE/M,MF,MU
	DIMENSION F1(-6:6),F2(-6:6),QCHA(6)
	PARAMETER (K0MAX=4,J0MAX=2)
	DIMENSION SF(K0MAX),SFD(K0MAX)
	DATA QCHA/-1.D0,2.D0,-1.D0,2.D0,-1.D0,2.D0/
	CALL FSTRU(X1,M*M,IH1, F1)
	CALL FSTRU(X2,M*M,IH2, F2)
C initialisation du tableau SF(K0MAX)
	ZERO = 0.D0
	CALL VEC_DINIT(SF,K0MAX,ZERO)
	CALL VEC_DINIT(SFD,K0MAX,ZERO)
C on commence les processus
	DO I=1,JNF
C	j0 = 1 : qi + qbi --> ph + ph
	  SF(1) = QCHA(I)**4*(F1(I)*F2(-I)) + SF(1)
	  SF(1+2) = QCHA(I)**4*(F1(-I)*F2(I)) + SF(1+2)
C	j0 = 2 : qi + g --> ph + ph
	  SF(2) = QCHA(I)**4*(F1(I)*F2(0) + 	    
     #             F1(-I)*F2(0)) + SF(2)
	  SF(2+2) = QCHA(I)**4*(F1(0)*F2(I) + 	    
     #               F1(0)*F2(-I)) + SF(2+2)
	ENDDO 	    
C on divise tout par x1 x2
	XX = 1.D0/(X1*X2)
	CALL VEC_DMULT_CONSTANT(SF,K0MAX,XX,SFD)	    
	RETURN
	END
