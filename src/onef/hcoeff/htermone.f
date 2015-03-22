c bug dans les parties virtuelles corrige le 26/07/96
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 1
C*****************************************************
	SUBROUTINE VEC_15ZO(SC,TC,UC,Z1,M,PTM,CQI,CQK,VPART15Z)
C PARTIE EN 1/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART15Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART15Z(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART15Z(I) = T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART15Z(I) = T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = 0
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART15Z(I) = T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = eqi**2*E(UC,SC,TC,N)*(-ANM4_GQ(z1)-2*A4_GQ(z1)*dlog(ptm)+2*A4
     #_GQ(z1)*dlog(M)+FAGQ(z1))/CF/2
	VPART15Z(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = 0
	VPART15Z(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = -eqi**2*E(TC,SC,UC,N)*(ANM4_QQ(z1)+2*A4_QQ(z1)*dlog(ptm)-2*A4
     #_QQ(z1)*dlog(M)-FAQQ(z1))
	VPART15Z(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = 0
	VPART15Z(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*(-ANM4_QQ(z1)-2*dlog(ptm)*A4_QQ(z1)+2*d
     #log(M)*A4_QQ(z1)+FAQQ(z1))
	VPART15Z(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = -eqi**2*E(TC,SC,UC,N)*(N-1)*(N+1)*(ANM4_QG(z1)+2*A4_QG(z1)*dl
     #og(ptm)-2*A4_QG(z1)*dlog(M)-FAQG(z1))/N
	VPART15Z(18) = T0
C
C
	RETURN
	END
	SUBROUTINE VEC_15LO(SC,TC,UC,Z1,CQI,CQK,VPART15L)
C PARTIE EN LOG(1-Z1)/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART15L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART15L(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART15L(I) = T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART15L(I) = T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = 0
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART15L(I) = T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = eqi**2*E(UC,SC,TC,N)*FBGQ(z1)/CF/2
	VPART15L(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = 0
	VPART15L(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FBQQ(z1)
	VPART15L(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = 0
	VPART15L(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*FBQQ(z1)
	VPART15L(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FBQG(z1)*(N-1)*(N+1)/N
	VPART15L(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_15DO(SC,TC,UC,M,PTM,CQI,CQK,VPART15D)
C PARTIE EN DELTA(1-Z1)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART15D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART15D(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART15D(I) = T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART15D(I) = T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = 0
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART15D(I) = T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = eqi**2*E(UC,SC,TC,N)*FCGQ(1.D0)/CF/2
	VPART15D(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = 0
	VPART15D(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = -eqi**2*E(TC,SC,UC,N)*(2*B_QQ()*dlog(ptm)-2*B_QQ()*dlog(M)-FC
     #QQ(1.D0))
	VPART15D(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = 0
	VPART15D(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = eqi**2*E(SC,TC,UC,N)*(2*B_QQ()*dlog(ptm)-2*B_QQ()*dlog(M)-FCQ
     #Q(1.D0))
	VPART15D(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FCQG(1.D0)*(N-1)*(N+1)/N
	VPART15D(18) = T0
C
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 2
C*****************************************************
	SUBROUTINE VEC_25ZO(SC,TC,UC,Z2,M,PTM,CQI,CQK,VPART25Z)
C PARTIE EN 1/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART25Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART25Z(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = -E(TC,SC,UC,N)*(ANM4_GQ(z2)+2*A4_GQ(z2)*dlog(ptm)-2*A4
     #_GQ(z2)*dlog(M)-FAGQ(z2))/CF/2
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART25Z(I) = CQI(I)**2*T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = -E(TC,SC,UC,N)*(ANM4_GQ(z2)+2*A4_GQ(z2)*dlog(ptm)-2*A4
     #_GQ(z2)*dlog(M)-FAGQ(z2))/CF/2
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART25Z(I) = CQI(I)**2*T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = 0
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART25Z(I) = T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*(-ANM4_GQ(z2)-2*A4_GQ(z2)*dlog(ptm)+2*A4
     #_GQ(z2)*dlog(M)+FAGQ(z2))/CF/2
	VPART25Z(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*(-ANM4_GQ(z2)-2*A4_GQ(z2)*dlog(ptm)+2*A4
     #_GQ(z2)*dlog(M)+FAGQ(z2))/CF/2
	VPART25Z(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*(-ANM4_GG(z2)-2*dlog(ptm)*A4_GG(z2)+2*dl
     #og(M)*A4_GG(z2)+FAGG(z2))
	VPART25Z(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = eqi**2*E(SC,TC,UC,N)*(N-1)*(N+1)*(ANM4_QG(z2)+2*A4_QG(z2)*dlo
     #g(ptm)-2*A4_QG(z2)*dlog(M)-FAQG(z2))/N
	VPART25Z(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = eqi**2*E(SC,TC,UC,N)*(ANM4_QQ(z2)+2*dlog(ptm)*A4_QQ(z2)-2*dlo
     #g(M)*A4_QQ(z2)-FAQQ(z2))
	VPART25Z(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = -eqi**2*E(UC,SC,TC,N)*(N-1)*(N+1)*(ANM4_QG(z2)+2*A4_QG(z2)*dl
     #og(ptm)-2*A4_QG(z2)*dlog(M)-FAQG(z2))/N
	VPART25Z(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_25LO(SC,TC,UC,Z2,CQI,CQK,VPART25L)
C PARTIE EN LOG(1-Z2)/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART25L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART25L(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = E(TC,SC,UC,N)*FBGQ(z2)/CF/2
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART25L(I) = CQI(I)**2*T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = E(TC,SC,UC,N)*FBGQ(z2)/CF/2
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART25L(I) = CQI(I)**2*T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = 0
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART25L(I) = T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FBGQ(z2)/CF/2
	VPART25L(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FBGQ(z2)/CF/2
	VPART25L(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FBGG(z2)
	VPART25L(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*FBQG(z2)*(N-1)*(N+1)/N
	VPART25L(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*FBQQ(z2)
	VPART25L(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = eqi**2*E(UC,SC,TC,N)*FBQG(z2)*(N-1)*(N+1)/N
	VPART25L(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_25DO(SC,TC,UC,M,PTM,CQI,CQK,VPART25D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART25D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART25D(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = E(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART25D(I) = CQI(I)**2*T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = E(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART25D(I) = CQI(I)**2*T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = 0
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART25D(I) = T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
	VPART25D(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
	VPART25D(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*(-2*dlog(ptm)*B_GG()+2*dlog(M)*B_GG()+FC
     #GG(1.D0))
	VPART25D(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*FCQG(1.D0)*(N-1)*(N+1)/N
	VPART25D(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = eqi**2*E(SC,TC,UC,N)*(2*B_QQ()*dlog(ptm)-2*B_QQ()*dlog(M)-FCQ
     #Q(1.D0))
	VPART25D(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = eqi**2*E(UC,SC,TC,N)*FCQG(1.D0)*(N-1)*(N+1)/N
	VPART25D(18) = T0
C
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 3
C*****************************************************
	SUBROUTINE VEC_35RO(SC,TC,UC,Z3,R,CQI,CQK,VPART35R)
C PARTIE EN LOG(R)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART35R(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART35R(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART35R(I) = T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART35R(I) = T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = A4_QG(z3)*dlog(R**2)/(1-z3)*E(SC,TC,UC,N)
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART35R(I) = CQI(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = 0
	VPART35R(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = A4_QG(z3)*dlog(R**2)/(1-z3)*eqi**2*E(SC,TC,UC,N)
	VPART35R(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = -A4_QQ(z3)*dlog(R**2)/(1-z3)*eqi**2*E(TC,SC,UC,N)
	VPART35R(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = -A4_GQ(z3)*dlog(R**2)/(1-z3)*eqi**2*E(TC,SC,UC,N)
	VPART35R(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = A4_GG(z3)*dlog(R**2)/(1-z3)*eqi**2*E(SC,TC,UC,N)
	VPART35R(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = 0
	VPART35R(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_35ZO(SC,TC,UC,Z3,MF,PT3,CQI,CQK,VPART35Z)
C PARTIE EN 1/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART35Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART35Z(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART35Z(I) = T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART35Z(I) = T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = E(SC,TC,UC,N)*(ANM4_QG(z3)+2*dlog(pt3)*A4_QG(z3)-2*dlo
     #g(Mf)*A4_QG(z3)-DAQG(z3))
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART35Z(I) = CQI(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = 0
	VPART35Z(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*(-ANM4_QG(z3)-2*dlog(pt3)*A4_QG(z3)+2*d
     #log(Mf)*A4_QG(z3)+DAQG(z3))
	VPART35Z(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = -eqi**2*E(TC,SC,UC,N)*(ANM4_QQ(z3)+2*dlog(pt3)*A4_QQ(z3)-2*dl
     #og(Mf)*A4_QQ(z3)-DAQQ(z3))
	VPART35Z(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = -eqi**2*E(TC,SC,UC,N)*(ANM4_GQ(z3)+2*A4_GQ(z3)*dlog(pt3)-2*A4
     #_GQ(z3)*dlog(Mf)-DAGQ(z3))
	VPART35Z(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*(-ANM4_GG(z3)-2*A4_GG(z3)*dlog(pt3)+2*A
     #4_GG(z3)*dlog(Mf)+DAGG(z3))
	VPART35Z(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = 0
	VPART35Z(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_35LO(SC,TC,UC,Z3,CQI,CQK,VPART35L)
C PARTIE EN LOG(1-Z3)/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART35L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART35L(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART35L(I) = T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART35L(I) = T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = E(SC,TC,UC,N)*(-DBQG(z3)+2*A4_QG(z3))
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART35L(I) = CQI(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = 0
	VPART35L(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*(DBQG(z3)-2*A4_QG(z3))
	VPART35L(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*(DBQQ(z3)-2*A4_QQ(z3))
	VPART35L(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = eqi**2*E(TC,SC,UC,N)*(DBGQ(z3)-2*A4_GQ(z3))
	VPART35L(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*(DBGG(z3)-2*A4_GG(z3))
	VPART35L(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = 0
	VPART35L(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_35DO(SC,TC,UC,MF,PT3,CQI,CQK,VPART35D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART35D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART35D(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART35D(I) = T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART35D(I) = T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = -E(SC,TC,UC,N)*DCQG(1.D0)
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART35D(I) = CQI(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = 0
	VPART35D(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = -eqi**2*E(SC,TC,UC,N)*DCQG(1.D0)
	VPART35D(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = eqi**2*E(TC,SC,UC,N)*(-2*dlog(pt3)*B_QQ()+2*dlog(Mf)*B_QQ()+D
     #CQQ(1.D0))
	VPART35D(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = eqi**2*E(TC,SC,UC,N)*DCGQ(1.D0)
	VPART35D(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = eqi**2*E(SC,TC,UC,N)*(2*B_GG()*dlog(pt3)-2*B_GG()*dlog(Mf)-DC
     #GG(1.D0))
	VPART35D(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = 0
	VPART35D(18) = T0
C
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 4
C*****************************************************
	SUBROUTINE VEC_45RO(SC,TC,UC,Z4,R,CQI,CQK,VPART45R)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART45R(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART45R(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = A4_PQ(z4)*dlog(R**2)/(1-z4)*A(SC,TC,UC,N)
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART45R(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = A4_PQ(z4)*dlog(R**2)/(1-z4)*A(UC,TC,SC,N)
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART45R(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = A4_PQ(z4)*dlog(R**2)/(1-z4)*A(UC,SC,TC,N)
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART45R(I) = CQK(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = eqi**2*A4_PQ(z4)*dlog(R**2)/(1-z4)*(A(SC,TC,UC,N)+A(SC,UC,TC,
     #N)+B(SC,TC,UC,N))
	VPART45R(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = eqi**2*A4_PQ(z4)*dlog(R**2)/(1-z4)*(A(UC,TC,SC,N)+A(UC,SC,TC,
     #N)+B(UC,TC,SC,N))
	VPART45R(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = 0
	VPART45R(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = -eqi**2*A4_PQ(z4)*dlog(R**2)/(1-z4)*C(UC,SC,TC,N)
	VPART45R(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = 0
	VPART45R(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = eqi**2*A4_PQ(z4)*dlog(R**2)/(1-z4)*C(SC,TC,UC,N)
	VPART45R(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_45ZO(SC,TC,UC,Z4,MF,PT4,CQI,CQK,VPART45Z)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART45Z(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART45Z(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = A(SC,TC,UC,N)*(ANM4_PQ(z4)+2*dlog(pt4)*A4_PQ(z4)-2*dlo
     #g(Mf)*A4_PQ(z4)-DAPQ(z4))
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART45Z(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = A(UC,TC,SC,N)*(ANM4_PQ(z4)+2*dlog(pt4)*A4_PQ(z4)-2*dlo
     #g(Mf)*A4_PQ(z4)-DAPQ(z4))
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART45Z(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = A(UC,SC,TC,N)*(ANM4_PQ(z4)+2*A4_PQ(z4)*dlog(pt4)-2*A4_
     #PQ(z4)*dlog(Mf)-DAPQ(z4))
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART45Z(I) = CQK(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = -eqi**2*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))*(-2*dlog(
     #pt4)*A4_PQ(z4)-ANM4_PQ(z4)+DAPQ(z4)+2*dlog(Mf)*A4_PQ(z4))
	VPART45Z(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = -eqi**2*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))*(2*A4_PQ(
     #z4)*dlog(Mf)-2*A4_PQ(z4)*dlog(pt4)-ANM4_PQ(z4)+DAPQ(z4))
	VPART45Z(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = 0
	VPART45Z(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = -eqi**2*C(UC,SC,TC,N)*(ANM4_PQ(z4)+2*A4_PQ(z4)*dlog(pt4)-2*A4
     #_PQ(z4)*dlog(Mf)-DAPQ(z4))
	VPART45Z(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = 0
	VPART45Z(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = -eqi**2*C(SC,TC,UC,N)*(-ANM4_PQ(z4)-2*A4_PQ(z4)*dlog(pt4)+2*A
     #4_PQ(z4)*dlog(Mf)+DAPQ(z4))
	VPART45Z(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_45LO(SC,TC,UC,Z4,CQI,CQK,VPART45L)
C PARTIE EN LOG(1-Z4)/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART45L(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART45L(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = A(SC,TC,UC,N)*(-DBPQ(z4)+2*A4_PQ(z4))
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART45L(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = A(UC,TC,SC,N)*(-DBPQ(z4)+2*A4_PQ(z4))
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART45L(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = A(UC,SC,TC,N)*(-DBPQ(z4)+2*A4_PQ(z4))
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART45L(I) = CQK(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = -eqi**2*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))*(-2*A4_PQ
     #(z4)+DBPQ(z4))
	VPART45L(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = eqi**2*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))*(-DBPQ(z4)
     #+2*A4_PQ(z4))
	VPART45L(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = 0
	VPART45L(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = -eqi**2*C(UC,SC,TC,N)*(-DBPQ(z4)+2*A4_PQ(z4))
	VPART45L(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = 0
	VPART45L(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = -eqi**2*C(SC,TC,UC,N)*(DBPQ(z4)-2*A4_PQ(z4))
	VPART45L(18) = T0
C
	RETURN
	END
	SUBROUTINE VEC_45DO(SC,TC,UC,MF,PT4,CQI,CQK,VPART45D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPART45D(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
	DO I = 1,J0MAX
	  VPART45D(I) = 0.D0
	ENDDO
C	qi + qk --> qi + ph
      t0 = -A(SC,TC,UC,N)*DCPQ(1.D0)
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPART45D(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbk --> qi + ph
      t0 = -A(UC,TC,SC,N)*DCPQ(1.D0)
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPART45D(I) = CQK(I)**2*T0
	ENDDO
C	qi + qbi --> qk + ph
      t0 = -A(UC,SC,TC,N)*DCPQ(1.D0)
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPART45D(I) = CQK(I)**2*T0
	ENDDO
C
	EQI = CQI(13)
C	j0 = 13 : qi + qi --> qi + ph
      t0 = -eqi**2*DCPQ(1.D0)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N)
     #)
	VPART45D(13) = T0
C	j0 = 14 : qi + qbi --> qi + ph
      t0 = -eqi**2*DCPQ(1.D0)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N)
     #)
	VPART45D(14) = T0
C	j0 = 15 : qi + g --> qi + ph
      t0 = 0
	VPART45D(15) = T0
C	j0 = 16 : qi + g --> g + ph
      t0 = eqi**2*C(UC,SC,TC,N)*DCPQ(1.D0)
	VPART45D(16) = T0
C	j0 = 17 : qi + qbi --> g + ph
      t0 = 0
	VPART45D(17) = T0
C	j0 = 18 : g + g --> qi + ph
      t0 = -eqi**2*C(SC,TC,UC,N)*DCPQ(1.D0)
	VPART45D(18) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE LA PARTIE VIRTUELLE
C*****************************************************
	SUBROUTINE VEC_VIO(S,SC,TC,UC,MU,PT,PTM,CQI,VPARTVI)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPARTVI(J0MAX)
	DIMENSION CQI(J0MAX)
	PI=4.D0*DATAN(1.D0)
	DO I = 1,J0MAX
	  VPARTVI(I) = 0.D0
	ENDDO
	MU2=MU**2
	EQI = CQI(13)
C	qi + qk --> qi + ph
      teraj = 0
      teres = 0
C	j0 = 1 : D + U --> D + ph
C	j0 = 2 : D + Dp --> D + ph
C	j0 = 3 : U + D --> U + ph
C	j0 = 4 : U + Up --> U + ph
	DO I=1,4
	  VPARTVI(I) = (TERAJ+TERES)
	ENDDO
C	qi + qbk --> qi + ph
      teraj = 0
      teres = 0
C	j0 = 5 : D + Ub --> D + ph
C	j0 = 6 : D + Dpb --> D + ph
C	j0 = 7 : U + Db --> U + ph
C	j0 = 8 : U + Upb --> U + ph
	DO I=5,8
	  VPARTVI(I) = (TERAJ+TERES)
	ENDDO
C	qi + qbi --> qk + ph
      teraj = 0
      teres = 0
C	j0 = 9 : D + Db --> U + ph
C	j0 = 10 : D + Db --> Dp + ph
C	j0 = 11 : U + Ub --> D + ph
C	j0 = 12 : U + Ub --> Up + ph
	DO I=9,12
	  VPARTVI(I) = (TERAJ+TERES)
	ENDDO
C	j0 = 13 : qi + qi --> qi + ph
      teraj = 0
      teres = 0
	VPARTVI(13)=(TERAJ+TERES)
C	j0 = 14 : qi + qbi --> qi + ph
      teraj = 0
      teres = 0
	VPARTVI(14)=(TERAJ+TERES)
C	j0 = 15 : qi + g --> qi + ph
      teraj = (B_QQ()+B_GG())*dlog(ptm**2/S)*eqi**2*E(TC,SC,UC,N)+B_QQ()
     #*dlog(pt**2/S)*eqi**2*E(TC,SC,UC,N)-((N**2-1)/N*dlog(-TC/S)+N*(dlo
     #g(-UC/S)+dlog(SC/S)-dlog(-TC/S)))*E(TC,UC,SC,N)*eqi**2*dlog(ptm**2
     #/S)-A4_QQ(1.D0)*eqi**2*E(TC,SC,UC,N)*(dlog(pt**2/S)**2/4-dlog(pt**
     #2/S)*dlog(ptm**2/S)/2)+(A4_QQ(1.D0)+A4_GG(1.D0))*eqi**2*E(TC,SC,UC
     #,N)*dlog(ptm**2/S)**2/4
      s1 = 4*(-11.D0/6.D0*N*dlog(MU2/S)+2.D0/3.D0*GTR*dlog(MU2/S))*(N**2
     #-1)/SC/UC*(UC**2+SC**2)*eqi**2
      s3 = (-2*N**2+2)*eqi**2*(-1/N*(TC**2*dlog(-TC/S)**2+(UC**2+TC**2)*
     #dlog(-UC/S)**2+(-2*UC**2-2*TC**2)*dlog(-TC/S)*dlog(-UC/S)-2*UC*SC*
     #dlog(-TC/S)+SC*(3*SC+2*UC)*dlog(-UC/S)+Pi**2*TC**2+Pi**2*UC**2/2+P
     #i**2*SC**2/2-7.D0/2.D0*UC**2-7.D0/2.D0*SC**2)+N*(3*SC**2*dlog(-UC/
     #S)+(Pi**2/2-7.D0/2.D0-dlog(-UC/S)*dlog(SC/S))*(UC**2+SC**2)))/SC/U
     #C
      s4 = (-2*N**2+2)*eqi**2*(-1/N*(TC**2*dlog(-TC/S)**2+(SC**2+TC**2)*
     #(-Pi**2+dlog(SC/S)**2)+(-2*SC**2-2*TC**2)*dlog(-TC/S)*dlog(SC/S)-2
     #*UC*SC*dlog(-TC/S)+UC*(3*UC+2*SC)*dlog(SC/S)+Pi**2*TC**2+Pi**2*UC*
     #*2/2+Pi**2*SC**2/2-7.D0/2.D0*SC**2-7.D0/2.D0*UC**2)+N*(3*UC**2*dlo
     #g(SC/S)+(Pi**2/2-7.D0/2.D0-dlog(-UC/S)*dlog(SC/S))*(UC**2+SC**2)))
     #/UC/SC
      s2 = s3+s4
      teres = s1+s2
	VPARTVI(15)=(TERAJ+TERES)
C	j0 = 16 : qi + g --> g + ph
      teraj = 0
      teres = 0
	VPARTVI(16)=(TERAJ+TERES)
C	j0 = 17 : qi + qbi --> g + ph
      teraj = -2*B_QQ()*dlog(ptm**2/S)*eqi**2*E(SC,TC,UC,N)-B_GG()*dlog(
     #pt**2/S)*eqi**2*E(SC,TC,UC,N)+((N**2-1)/N*dlog(SC/S)+N*(dlog(-UC/S
     #)+dlog(-TC/S)-dlog(SC/S)))*E(SC,UC,TC,N)*eqi**2*dlog(ptm**2/S)+A4_
     #GG(1.D0)*eqi**2*E(SC,TC,UC,N)*(dlog(pt**2/S)**2/4-dlog(pt**2/S)*dl
     #og(ptm**2/S)/2)-A4_QQ(1.D0)*eqi**2*E(SC,TC,UC,N)*dlog(ptm**2/S)**2
     #/2
      s1 = 4*(11.D0/6.D0*N*dlog(MU2/S)-2.D0/3.D0*GTR*dlog(MU2/S))*(N**2-
     #1)/TC/UC*(UC**2+TC**2)*eqi**2
      s3 = (2*N**2-2)*eqi**2*(-1/N*(SC**2*(-Pi**2+dlog(SC/S)**2)+(UC**2+
     #SC**2)*dlog(-UC/S)**2+(-2*UC**2-2*SC**2)*dlog(SC/S)*dlog(-UC/S)-2*
     #UC*TC*dlog(SC/S)+TC*(3*TC+2*UC)*dlog(-UC/S)+Pi**2*SC**2+Pi**2*UC**
     #2/2+Pi**2*TC**2/2-7.D0/2.D0*UC**2-7.D0/2.D0*TC**2)+N*(3*TC**2*dlog
     #(-UC/S)+(Pi**2/2-7.D0/2.D0-dlog(-UC/S)*dlog(-TC/S))*(UC**2+TC**2))
     #)/TC/UC
      s4 = (2*N**2-2)*eqi**2*(-1/N*(SC**2*(-Pi**2+dlog(SC/S)**2)+(TC**2+
     #SC**2)*dlog(-TC/S)**2+(-2*TC**2-2*SC**2)*dlog(SC/S)*dlog(-TC/S)-2*
     #UC*TC*dlog(SC/S)+UC*(3*UC+2*TC)*dlog(-TC/S)+Pi**2*SC**2+Pi**2*UC**
     #2/2+Pi**2*TC**2/2-7.D0/2.D0*TC**2-7.D0/2.D0*UC**2)+N*(3*UC**2*dlog
     #(-TC/S)+(Pi**2/2-7.D0/2.D0-dlog(-UC/S)*dlog(-TC/S))*(UC**2+TC**2))
     #)/UC/TC
      s2 = s3+s4
      teres = s1+s2
	VPARTVI(17)=(TERAJ+TERES)
C	j0 = 18 : g + g --> qi + ph
      teraj = 0
      teres = 0
	VPARTVI(18)=(TERAJ+TERES)
C
	RETURN
	END
	SUBROUTINE VEC_VI2O(S,SC,TC,UC,FI,CQI,CQK,VPARTVI2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=18)
	DIMENSION VPARTVI2(J0MAX),HH34(J0MAX)
	DIMENSION CQI(J0MAX),CQK(J0MAX)
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
	CALL VEC_H34O(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45,
     #	CQI,CQK,HH34)
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
C pour production d'un photon	
	SUBROUTINE STRFRAO(X1,IH1,X2,IH2,X3,IH3,SFD)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/GDF/JNF
	COMMON/SCALE/M,MF,MU
	DIMENSION F1(-6:6),F2(-6:6),D3(-6:6)
	PARAMETER (K0MAX=36,J0MAX=18)
c qcharge est la valeur absolue de la charge des quarks fois 3 au carre
	DIMENSION SF(K0MAX),SFD(K0MAX),QCHARGE(6)
	DATA QCHARGE/1.D0,4.D0,1.D0,4.D0,1.D0,4.D0/
	CALL FSTRU(X1,M*M,IH1, F1)
	CALL FSTRU(X2,M*M,IH2, F2)
	CALL DFRAG(X3,MF*MF,IH3, D3)
C initialisation du tableau SF(K0MAX)
	ZERO = 0.D0
	CALL VEC_DINIT(SF,K0MAX,ZERO)
	CALL VEC_DINIT(SFD,K0MAX,ZERO)
C on commence les processus
	DO I=1,JNF,2
	  DO K=2,JNF,2
C	qi + qk --> qi + ph
C	j0 = 1 : D + U --> D + ph 
	    SF(1) = F1(I)*F2(K)*D3(I) + 	    
     #              F1(-I)*F2(-K)*D3(-I) + SF(1)
	    SF(1+18) = F1(K)*F2(I)*D3(I) + 	    
     #                 F1(-K)*F2(-I)*D3(-I) + SF(1+18)
C	j0 = 3 : U + D --> U + ph
	    SF(3) = F1(K)*F2(I)*D3(K) + 	    
     #              F1(-K)*F2(-I)*D3(-K) + SF(3)
	    SF(3+18) = F1(I)*F2(K)*D3(K) + 	    
     #                 F1(-I)*F2(-K)*D3(-K) + SF(3+18)
C	qi + qbk --> qi + ph
C	j0 = 5 : D + Ub --> D + ph
	    SF(5) = F1(I)*F2(-K)*D3(I) + 	    
     #              F1(-I)*F2(K)*D3(-I) + SF(5)
	    SF(5+18) = F1(-K)*F2(I)*D3(I) + 	    
     #                 F1(K)*F2(-I)*D3(-I) + SF(5+18)
C	j0 = 7 : U + Db --> U + ph
	    SF(7) = F1(K)*F2(-I)*D3(K) + 	    
     #              F1(-K)*F2(I)*D3(-K) + SF(7)
	    SF(7+18) = F1(-I)*F2(K)*D3(K) + 	    
     #                 F1(I)*F2(-K)*D3(-K) + SF(7+18)
C	qi + qbi --> qk + ph
C	j0 = 9 : D + Db --> U + ph
	      SF(9) = F1(I)*F2(-I)*D3(K) + 
     #                F1(-I)*F2(I)*D3(-K) + SF(9)
	      SF(9+18) = F1(-I)*F2(I)*D3(K) + 
     #                   F1(I)*F2(-I)*D3(-K) + SF(9+18)
C	j0 = 11 : U + Ub --> D + ph
	      SF(11) = F1(K)*F2(-K)*D3(I) + 	    
     #                F1(-K)*F2(K)*D3(-I) + SF(11)
	      SF(11+18) = F1(-K)*F2(K)*D3(I) + 	    
     #                   F1(K)*F2(-K)*D3(-I) + SF(11+18)
	  ENDDO
	ENDDO
	DO I=1,JNF,2
	  DO K=I+2,JNF,2
C	j0 = 2 : D + Dp --> D + ph
	    SF(2) = F1(I)*F2(K)*D3(I) + F1(K)*F2(I)*D3(K) +	    
     #              F1(-I)*F2(-K)*D3(-I) + F1(-K)*F2(-I)*D3(-K) + SF(2)
	    SF(2+18) = F1(K)*F2(I)*D3(I) + F1(I)*F2(K)*D3(K) + 	    
     #                 F1(-K)*F2(-I)*D3(-I) + F1(-I)*F2(-K)*D3(-K) + 
     #	               SF(2+18)
C	j0 = 4 : U + Up --> U + ph
	    SF(4) = F1(I+1)*F2(K+1)*D3(I+1) + F1(K+1)*F2(I+1)*D3(K+1) +	    
     #              F1(-I-1)*F2(-K-1)*D3(-I-1) + 
     #	            F1(-K-1)*F2(-I-1)*D3(-K-1) + SF(4)
	    SF(4+18) = F1(K+1)*F2(I+1)*D3(I+1) + 
     #	               F1(I+1)*F2(K+1)*D3(K+1) +	    
     #                 F1(-K-1)*F2(-I-1)*D3(-I-1) +
     #	               F1(-I-1)*F2(-K-1)*D3(-K-1) + SF(4+18)
C	j0 = 6 : D + Dpb --> D + ph
	    SF(6) = F1(I)*F2(-K)*D3(I) + F1(-K)*F2(I)*D3(-K) +	    
     #              F1(-I)*F2(K)*D3(-I) + F1(K)*F2(-I)*D3(K) + SF(6)
	    SF(6+18) = F1(-K)*F2(I)*D3(I) + F1(I)*F2(-K)*D3(-K) +	    
     #                 F1(K)*F2(-I)*D3(-I) + F1(-I)*F2(K)*D3(K) + 
     #	               SF(6+18)
C	j0 = 8 : U + Upb --> U + ph
	    SF(8) = F1(I+1)*F2(-K-1)*D3(I+1) + 
     #	            F1(-K-1)*F2(I+1)*D3(-K-1) + 	    
     #              F1(-I-1)*F2(K+1)*D3(-I-1) + 
     #	            F1(K+1)*F2(-I-1)*D3(K+1) + SF(8)
	    SF(8+18) = F1(-K-1)*F2(I+1)*D3(I+1) + 
     #	               F1(I+1)*F2(-K-1)*D3(-K-1) +	    
     #                 F1(K+1)*F2(-I-1)*D3(-I-1) +
     #                 F1(-I-1)*F2(K+1)*D3(K+1) + SF(8+18)
C	j0 = 10 : D + Db --> Dp + ph
	      SF(10) = F1(I)*F2(-I)*D3(K) + 
     #                F1(-I)*F2(I)*D3(-K) + SF(10)
	      SF(10+18) = F1(-I)*F2(I)*D3(K) + 
     #                   F1(I)*F2(-I)*D3(-K) + SF(10+18)
C	j0 = 12 : U + Ub --> Up + ph
	      SF(12) = F1(K+1)*F2(-K-1)*D3(I+1) + 	    
     #                F1(-K-1)*F2(K+1)*D3(-I-1) + SF(12)
	      SF(12+18) = F1(-K-1)*F2(K+1)*D3(I+1) + 	    
     #                   F1(K+1)*F2(-K-1)*D3(-I-1) + SF(12+18)
	  ENDDO
	ENDDO
	DO I=1,JNF
C	j0 = 13 : qi + qi --> qi + ph
	  SF(13) = QCHARGE(I)*(F1(I)*F2(I)*D3(I) + 	    
     #             F1(-I)*F2(-I)*D3(-I))/2.D0 + SF(13)
C	j0 = 14 : qi + qbi --> qi + ph
	  SF(14) = QCHARGE(I)*(F1(I)*F2(-I)*D3(I) + 	    
     #            F1(-I)*F2(I)*D3(-I)) + SF(14)
	  SF(14+18) = QCHARGE(I)*(F1(-I)*F2(I)*D3(I) + 	    
     #              F1(I)*F2(-I)*D3(-I)) + SF(14+18)
C	j0 = 15 : qi + g --> qi + ph
	  SF(15) = QCHARGE(I)*(F1(I)*F2(0)*D3(I) + 	    
     #            F1(-I)*F2(0)*D3(-I)) + SF(15)
	  SF(15+18) = QCHARGE(I)*(F1(0)*F2(I)*D3(I) + 	    
     #               F1(0)*F2(-I)*D3(-I)) + SF(15+18)
C	j0 = 16 : qi + g --> g + ph
	  SF(16) = QCHARGE(I)*(F1(I)*F2(0)*D3(0) + 	    
     #              F1(-I)*F2(0)*D3(0)) + SF(16)
	  SF(16+18) = QCHARGE(I)*(F1(0)*F2(I)*D3(0) + 	    
     #              F1(0)*F2(-I)*D3(0)) + SF(16+18)
C	j0 = 17 : qi + qbi --> g + ph
	  SF(17) = QCHARGE(I)*F1(I)*F2(-I)*D3(0) + SF(17)
	  SF(17+18) = QCHARGE(I)*F1(-I)*F2(I)*D3(0) + SF(17+18)
C	j0 = 18 : g + g --> qi + ph
	  SF(18) = QCHARGE(I)*(F1(0)*F2(0)*D3(I) + 	    
     #             F1(0)*F2(0)*D3(-I))/2.d0 + SF(18)
	ENDDO 	    
C
     	SF(13+18) = SF(13)
C
     	SF(18+18) = SF(18)
C
C on divise tout par x1 x2
	XX = 1.D0/(X1*X2)
	CALL VEC_DMULT_CONSTANT(SF,K0MAX,XX,SFD)	    
	RETURN
	END
C
