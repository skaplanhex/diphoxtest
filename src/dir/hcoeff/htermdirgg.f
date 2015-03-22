c partie pour la boite g g --> photon photon
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 1
C*****************************************************
	SUBROUTINE VEC_15ZGD(SC,TC,UC,Z1,M,PTM,VPART15Z)
C PARTIE EN 1/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART15Z(J0MAX)
	DO I = 1,J0MAX
	  VPART15Z(I) = 0.D0
	ENDDO
c l'equivalent de eqi**4*F(s,t,u) pour gg est la fonction ggpp(s,t,u)
C	j0 = 1 : g + g --> ph + ph
      t0 = GGPP(SC,TC,UC)*(ANM4_GG(z1)+2*A4_GG(z1)*dlog(ptm)-2*A4_
     #GG(z1)*dlog(M)-FAGG(z1))
	VPART15Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = GGPP(SC,TC,UC)*(ANM4_GQ(z1)+2*A4_GQ(z1)*dlog(ptm)-2*A4_
     #GQ(z1)*dlog(M)-FAGQ(z1))*N/(N**2-1.D0)
	VPART15Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_15LGD(SC,TC,UC,Z1,VPART15L)
C PARTIE EN LOG(1-Z1)/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART15L(J0MAX)
	DO I = 1,J0MAX
	  VPART15L(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = -GGPP(SC,TC,UC)*FBGG(z1)
	VPART15L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -GGPP(SC,TC,UC)*FBGQ(z1)*N/(N**2-1.D0)
 	VPART15L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_15DGD(SC,TC,UC,M,PTM,VPART15D)
C PARTIE EN DELTA(1-Z1)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART15D(J0MAX)
	DO I = 1,J0MAX
	  VPART15D(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = -GGPP(SC,TC,UC)*(-2*B_GG()*dlog(ptm)+2*B_GG()*dlog(M)+F
     #CGG(1.D0))
	VPART15D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = -GGPP(SC,TC,UC)*FCGQ(1.D0)*N/(N**2-1.D0)
	VPART15D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 2
C*****************************************************
	SUBROUTINE VEC_25ZGD(SC,TC,UC,Z2,M,PTM,VPART25Z)
C PARTIE EN 1/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	COMMON/AURENCHE/IAUREN
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART25Z(J0MAX)
	DO I = 1,J0MAX
	  VPART25Z(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = GGPP(SC,TC,UC)*(ANM4_GG(z2)+2*A4_GG(z2)*dlog(ptm)-2*A4_
     #GG(z2)*dlog(M)-FAGG(z2))
	VPART25Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART25Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_25LGD(SC,TC,UC,Z2,VPART25L)
C PARTIE EN LOG(1-Z2)/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART25L(J0MAX)
	DO I = 1,J0MAX
	  VPART25L(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = -GGPP(SC,TC,UC)*FBGG(z2)
	VPART25L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART25L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_25DGD(SC,TC,UC,M,PTM,VPART25D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART25D(J0MAX)
	DO I = 1,J0MAX
	  VPART25D(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = -GGPP(SC,TC,UC)*(-2*B_GG()*dlog(ptm)+2*B_GG()*dlog(M)+F
     #CGG(1.D0))
	VPART25D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART25D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 3
C*****************************************************
	SUBROUTINE VEC_35RGD(SC,TC,UC,Z3,R,VPART35R)
C PARTIE EN LOG(R)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35R(J0MAX)
	DO I = 1,J0MAX
	  VPART35R(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART35R(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART35R(2) = T0
	RETURN
	END
	SUBROUTINE VEC_35ZGD(SC,TC,UC,Z3,MF,PT3,VPART35Z)
C PARTIE EN 1/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35Z(J0MAX)
	DO I = 1,J0MAX
	  VPART35Z(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART35Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART35Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_35LGD(SC,TC,UC,Z3,VPART35L)
C PARTIE EN LOG(1-Z3)/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35L(J0MAX)
	DO I = 1,J0MAX
	  VPART35L(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART35L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART35L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_35DGD(SC,TC,UC,MF,PT3,VPART35D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART35D(J0MAX)
	DO I = 1,J0MAX
	  VPART35D(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART35D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART35D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 4
C*****************************************************
	SUBROUTINE VEC_45RGD(SC,TC,UC,Z4,R,VPART45R)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45R(J0MAX)
	DO I = 1,J0MAX
	  VPART45R(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART45R(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART45R(2) = T0
	RETURN
	END
	SUBROUTINE VEC_45ZGD(SC,TC,UC,Z4,MF,PT4,VPART45Z)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45Z(J0MAX)
	DO I = 1,J0MAX
	  VPART45Z(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART45Z(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART45Z(2) = T0
	RETURN
	END
	SUBROUTINE VEC_45LGD(SC,TC,UC,Z4,VPART45L)
C PARTIE EN LOG(1-Z4)/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45L(J0MAX)
	DO I = 1,J0MAX
	  VPART45L(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART45L(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART45L(2) = T0
	RETURN
	END
	SUBROUTINE VEC_45DGD(SC,TC,UC,MF,PT4,VPART45D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPART45D(J0MAX)
	DO I = 1,J0MAX
	  VPART45D(I) = 0.D0
	ENDDO
C	j0 = 1 : g + g --> ph + ph
      t0 = 0
	VPART45D(1) = T0
C	j0 = 2 : qi + g --> ph + ph
      t0 = 0
	VPART45D(2) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE LA PARTIE VIRTUELLE
C*****************************************************
	SUBROUTINE VEC_VIGD(S,SC,TC,UC,MU,PT,PTM,VPARTVI)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2,EQI=0.333333333333333)
	DIMENSION VPARTVI(J0MAX)
	PI=4.D0*DATAN(1.D0)
	DO I = 1,J0MAX
	  VPARTVI(I) = 0.D0
	ENDDO
	MU2=MU**2
C	j0 = 1 : g + g --> ph + ph
C termes provenant des parties collineaires
      teraj = -2*B_GG()*dlog(ptm**2/S)*GGPP(SC,TC,UC)+A4_GG(1.D0)*d
     #log(SC/S)*GGPP(SC,TC,UC)*dlog(ptm**2/S)-A4_GG(1.D0)*GGPP
     #(SC,TC,UC)*dlog(ptm**2/S)**2/2
C termes finies d'Ellis et Sexton (ils ne sont pas connus dnas ce cas)
      teres = 0
	VPARTVI(1)=(TERAJ+TERES)
C	j0 = 2 : qi + g --> ph + ph
      teraj = 0
      teres = 0
	VPARTVI(2)=(TERAJ+TERES)
	RETURN
	END
	SUBROUTINE VEC_VI2GD(S,SC,TC,UC,YS,FI,VPARTVI2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=2)
	DIMENSION VPARTVI2(J0MAX)
	DO I = 1,J0MAX
	  VPARTVI2(I) = 0.D0
	ENDDO
	RETURN
	END
C*****************************************************
C	
	SUBROUTINE STRFRAGD(X1,IH1,X2,IH2,SFD)
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
C	j0 = 1 : g + g --> ph + ph
	SF(1) = (F1(0)*F2(0))/2.D0
	SF(1+2) = SF(1)
C	j0 = 2 : qi + g --> ph + ph
	DO I=1,JNF
	  SF(2) = (F1(I)*F2(0) + 	    
     #             F1(-I)*F2(0)) + SF(2)
	  SF(2+2) = (F1(0)*F2(I) + 	    
     #               F1(0)*F2(-I)) + SF(2+2)
	ENDDO 	    
C on divise tout par x1 x2
	XX = 1.D0/(X1*X2)
	CALL VEC_DMULT_CONSTANT(SF,K0MAX,XX,SFD)	    
	RETURN
	END
