C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 1
C*****************************************************
	SUBROUTINE VEC_15ZT(SC,TC,UC,Z1,M,PTM,VPART15Z)
C PARTIE EN 1/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART15Z(J0MAX)
	DO I = 1,J0MAX
	  VPART15Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A(SC,TC,UC,N)*(ANM4_QQ(z1)+2*A4_QQ(z1)*dlog(ptm)-2*A4_QQ(z1)*
     #dlog(M)-FAQQ(z1))
	VPART15Z(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = 0
	VPART15Z(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = A(UC,TC,SC,N)*(ANM4_QQ(z1)+2*dlog(ptm)*A4_QQ(z1)-2*dlog(M)*A4
     #_QQ(z1)-FAQQ(z1))
	VPART15Z(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = 0
	VPART15Z(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A(UC,SC,TC,N)*(ANM4_QQ(z1)+2*A4_QQ(z1)*dlog(ptm)-2*A4_QQ(z1)*
     #dlog(M)-FAQQ(z1))
	VPART15Z(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = 0
	VPART15Z(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 0
	VPART15Z(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 0
	VPART15Z(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = -C(SC,TC,UC,N)*(-ANM4_GQ(z1)-2*A4_GQ(z1)*dlog(ptm)+2*A4_GQ(z1
     #)*dlog(M)+FAGQ(z1))/CF/2
	VPART15Z(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = (ANM4_QQ(z1)+2*A4_QQ(z1)*dlog(ptm)-2*A4_QQ(z1)*dlog(M)-FAQQ(z
     #1))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART15Z(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = C(UC,SC,TC,N)*(-ANM4_GQ(z1)-2*A4_GQ(z1)*dlog(ptm)+2*A4_GQ(z1)
     #*dlog(M)+FAGQ(z1))/CF/2
	VPART15Z(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(-ANM4_QQ(z1)-2*A4_QQ(z1)*dlog(ptm)+2*A4_QQ(z1)*dlog(M)+FAQQ
     #(z1))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART15Z(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = 0
	VPART15Z(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = 0
	VPART15Z(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -C(SC,TC,UC,N)*(-ANM4_GQ(z1)-2*A4_GQ(z1)*dlog(ptm)+2*A4_GQ(z1
     #)*dlog(M)+FAGQ(z1))/CF/2
	VPART15Z(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-ANM4_QQ(z1)-2*A4_QQ(z1)*dlog(ptm)+2*A4_QQ(z1
     #)*dlog(M)+FAQQ(z1))
	VPART15Z(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = -C(TC,SC,UC,N)*(ANM4_QQ(z1)+2*A4_QQ(z1)*dlog(ptm)-2*A4_QQ(z1)
     #*dlog(M)-FAQQ(z1))
	VPART15Z(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = -N*Dg(SC,TC,UC,N)*(-ANM4_GQ(z1)-2*A4_GQ(z1)*dlog(ptm)+2*A4_GQ
     #(z1)*dlog(M)+FAGQ(z1))/(N-1)/(N+1)
	VPART15Z(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = C(SC,TC,UC,N)*(ANM4_GG(z1)+2*dlog(ptm)*A4_GG(z1)-2*dlog(M)*A4
     #_GG(z1)-FAGG(z1))
	VPART15Z(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = C(TC,SC,UC,N)*(N-1)*(N+1)*(-ANM4_QG(z1)-2*A4_QG(z1)*dlog(ptm)
     #+2*A4_QG(z1)*dlog(M)+FAQG(z1))/N
	VPART15Z(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*(-ANM4_GG(z1)-2*dlog(ptm)*A4_GG(z1)+2*dlog(M)
     #*A4_GG(z1)+FAGG(z1))
	VPART15Z(21) = T0
	RETURN
	END
	SUBROUTINE VEC_15LT(SC,TC,UC,Z1,VPART15L)
C PARTIE EN LOG(1-Z1)/(1-Z1)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART15L(J0MAX)
	DO I = 1,J0MAX
	  VPART15L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = -A(SC,TC,UC,N)*FBQQ(z1)
	VPART15L(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = 0
	VPART15L(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = -A(UC,TC,SC,N)*FBQQ(z1)
	VPART15L(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = 0
	VPART15L(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = -A(UC,SC,TC,N)*FBQQ(z1)
	VPART15L(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = 0
	VPART15L(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 0
	VPART15L(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 0
	VPART15L(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = -C(SC,TC,UC,N)*FBGQ(z1)/CF/2
	VPART15L(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = -FBQQ(z1)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART15L(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = C(UC,SC,TC,N)*FBGQ(z1)/CF/2
	VPART15L(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -FBQQ(z1)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART15L(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = 0
	VPART15L(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = 0
	VPART15L(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -C(SC,TC,UC,N)*FBGQ(z1)/CF/2
	VPART15L(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*FBQQ(z1)
	VPART15L(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = C(TC,SC,UC,N)*FBQQ(z1)
	VPART15L(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = -N*Dg(SC,TC,UC,N)*FBGQ(z1)/(N-1)/(N+1)
	VPART15L(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = -C(SC,TC,UC,N)*FBGG(z1)
	VPART15L(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = C(TC,SC,UC,N)*FBQG(z1)*(N-1)*(N+1)/N
	VPART15L(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*FBGG(z1)
	VPART15L(21) = T0
	RETURN
	END
	SUBROUTINE VEC_15DT(SC,TC,UC,M,PTM,VPART15D)
C PARTIE EN DELTA(1-Z1)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART15D(J0MAX)
	DO I = 1,J0MAX
	  VPART15D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = -A(SC,TC,UC,N)*(-2*dlog(ptm)*B_QQ()+2*dlog(M)*B_QQ()+FCQQ(1.D
     #0))
	VPART15D(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = 0
	VPART15D(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = -A(UC,TC,SC,N)*(-2*dlog(ptm)*B_QQ()+2*dlog(M)*B_QQ()+FCQQ(1.D
     #0))
	VPART15D(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = 0
	VPART15D(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A(UC,SC,TC,N)*(2*B_QQ()*dlog(ptm)-2*B_QQ()*dlog(M)-FCQQ(1.D0)
     #)
	VPART15D(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = 0
	VPART15D(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 0
	VPART15D(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 0
	VPART15D(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = -C(SC,TC,UC,N)*FCGQ(1.D0)/CF/2
	VPART15D(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = (A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))*(-FCQQ(1.D0)-2*B_
     #QQ()*dlog(M)+2*B_QQ()*dlog(ptm))
	VPART15D(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = C(UC,SC,TC,N)*FCGQ(1.D0)/CF/2
	VPART15D(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(-2*B_QQ()*dlog(ptm)+2*B_QQ()*dlog(M)+FCQQ(1.D0))*(A(UC,TC,S
     #C,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART15D(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = 0
	VPART15D(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = 0
	VPART15D(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -C(SC,TC,UC,N)*FCGQ(1.D0)/CF/2
	VPART15D(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-2*B_QQ()*dlog(ptm)+2*B_QQ()*dlog(M)+FCQQ(1.D
     #0))
	VPART15D(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = C(TC,SC,UC,N)*(-2*dlog(ptm)*B_QQ()+2*dlog(M)*B_QQ()+FCQQ(1.D0
     #))
	VPART15D(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = -N*Dg(SC,TC,UC,N)*FCGQ(1.D0)/(N-1)/(N+1)
	VPART15D(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = C(SC,TC,UC,N)*(2*B_GG()*dlog(ptm)-2*B_GG()*dlog(M)-FCGG(1.D0)
     #)
	VPART15D(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = C(TC,SC,UC,N)*FCQG(1.D0)*(N-1)*(N+1)/N
	VPART15D(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*(-2*B_GG()*dlog(ptm)+2*B_GG()*dlog(M)+FCGG(1.
     #D0))
	VPART15D(21) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 2
C*****************************************************
	SUBROUTINE VEC_25ZT(SC,TC,UC,Z2,M,PTM,VPART25Z)
C PARTIE EN 1/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART25Z(J0MAX)
	DO I = 1,J0MAX
	  VPART25Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = -A(SC,TC,UC,N)*(-ANM4_QQ(z2)-2*dlog(ptm)*A4_QQ(z2)+2*dlog(M)*
     #A4_QQ(z2)+FAQQ(z2))
	VPART25Z(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = -C(TC,SC,UC,N)*(ANM4_GQ(z2)+2*A4_GQ(z2)*dlog(ptm)-2*A4_GQ(z2)
     #*dlog(M)-FAGQ(z2))/CF/2
	VPART25Z(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = A(UC,TC,SC,N)*(ANM4_QQ(z2)+2*A4_QQ(z2)*dlog(ptm)-2*A4_QQ(z2)*
     #dlog(M)-FAQQ(z2))
	VPART25Z(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = -C(TC,SC,UC,N)*(ANM4_GQ(z2)+2*A4_GQ(z2)*dlog(ptm)-2*A4_GQ(z2)
     #*dlog(M)-FAGQ(z2))/CF/2
	VPART25Z(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A(UC,SC,TC,N)*(ANM4_QQ(z2)+2*A4_QQ(z2)*dlog(ptm)-2*A4_QQ(z2)*
     #dlog(M)-FAQQ(z2))
	VPART25Z(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = 0
	VPART25Z(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 2*CF*A(SC,TC,UC,N)*(ANM4_QG(z2)+2*A4_QG(z2)*dlog(ptm)-2*A4_QG
     #(z2)*dlog(M)-FAQG(z2))
	VPART25Z(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 2*CF*A(UC,TC,SC,N)*(ANM4_QG(z2)+2*A4_QG(z2)*dlog(ptm)-2*A4_QG
     #(z2)*dlog(M)-FAQG(z2))
	VPART25Z(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 2*CF*A(UC,SC,TC,N)*(ANM4_QG(z2)+2*A4_QG(z2)*dlog(ptm)-2*A4_QG
     #(z2)*dlog(M)-FAQG(z2))
	VPART25Z(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = (ANM4_QQ(z2)+2*A4_QQ(z2)*dlog(ptm)-2*A4_QQ(z2)*dlog(M)-FAQQ(z
     #2))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART25Z(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = -C(TC,SC,UC,N)*(ANM4_GQ(z2)+2*A4_GQ(z2)*dlog(ptm)-2*A4_GQ(z2)
     #*dlog(M)-FAGQ(z2))/CF/2
	VPART25Z(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(-ANM4_QQ(z2)-2*A4_QQ(z2)*dlog(ptm)+2*A4_QQ(z2)*dlog(M)+FAQQ
     #(z2))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART25Z(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = -C(TC,SC,UC,N)*(ANM4_GQ(z2)+2*A4_GQ(z2)*dlog(ptm)-2*A4_GQ(z2)
     #*dlog(M)-FAGQ(z2))/CF/2
	VPART25Z(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = -2*CF*(-ANM4_QG(z2)-2*A4_QG(z2)*dlog(ptm)+2*A4_QG(z2)*dlog(M)
     #+FAQG(z2))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART25Z(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -2*CF*(-ANM4_QG(z2)-2*A4_QG(z2)*dlog(ptm)+2*A4_QG(z2)*dlog(M)
     #+FAQG(z2))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART25Z(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-ANM4_QQ(z2)-2*A4_QQ(z2)*dlog(ptm)+2*A4_QQ(z2
     #)*dlog(M)+FAQQ(z2))
	VPART25Z(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = -C(TC,SC,UC,N)*(ANM4_GG(z2)+2*A4_GG(z2)*dlog(ptm)-2*A4_GG(z2)
     #*dlog(M)-FAGG(z2))
	VPART25Z(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = C(SC,TC,UC,N)*(N-1)*(N+1)*(ANM4_QG(z2)+2*A4_QG(z2)*dlog(ptm)-
     #2*A4_QG(z2)*dlog(M)-FAQG(z2))/N
	VPART25Z(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = C(SC,TC,UC,N)*(ANM4_GG(z2)+2*dlog(ptm)*A4_GG(z2)-2*dlog(M)*A4
     #_GG(z2)-FAGG(z2))
	VPART25Z(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = -C(UC,SC,TC,N)*(N-1)*(N+1)*(ANM4_QG(z2)+2*A4_QG(z2)*dlog(ptm)
     #-2*A4_QG(z2)*dlog(M)-FAQG(z2))/N
	VPART25Z(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = Dg(SC,TC,UC,N)*(ANM4_GG(z2)+2*dlog(ptm)*A4_GG(z2)-2*dlog(M)*A
     #4_GG(z2)-FAGG(z2))
	VPART25Z(21) = T0
	RETURN
	END
	SUBROUTINE VEC_25LT(SC,TC,UC,Z2,VPART25L)
C PARTIE EN LOG(1-Z2)/(1-Z2)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART25L(J0MAX)
	DO I = 1,J0MAX
	  VPART25L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = -A(SC,TC,UC,N)*FBQQ(z2)
	VPART25L(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = C(TC,SC,UC,N)*FBGQ(z2)/CF/2
	VPART25L(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = -A(UC,TC,SC,N)*FBQQ(z2)
	VPART25L(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = C(TC,SC,UC,N)*FBGQ(z2)/CF/2
	VPART25L(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = -A(UC,SC,TC,N)*FBQQ(z2)
	VPART25L(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = 0
	VPART25L(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = -2*CF*A(SC,TC,UC,N)*FBQG(z2)
	VPART25L(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = -2*CF*A(UC,TC,SC,N)*FBQG(z2)
	VPART25L(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = -2*CF*A(UC,SC,TC,N)*FBQG(z2)
	VPART25L(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = -FBQQ(z2)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART25L(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = C(TC,SC,UC,N)*FBGQ(z2)/CF/2
	VPART25L(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -FBQQ(z2)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART25L(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = C(TC,SC,UC,N)*FBGQ(z2)/CF/2
	VPART25L(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = -2*CF*FBQG(z2)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART25L(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -2*CF*FBQG(z2)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART25L(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*FBQQ(z2)
	VPART25L(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = C(TC,SC,UC,N)*FBGG(z2)
	VPART25L(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = -C(SC,TC,UC,N)*FBQG(z2)*(N-1)*(N+1)/N
	VPART25L(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = -C(SC,TC,UC,N)*FBGG(z2)
	VPART25L(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = C(UC,SC,TC,N)*FBQG(z2)*(N-1)*(N+1)/N
	VPART25L(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*FBGG(z2)
	VPART25L(21) = T0
	RETURN
	END
	SUBROUTINE VEC_25DT(SC,TC,UC,M,PTM,VPART25D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART25D(J0MAX)
	DO I = 1,J0MAX
	  VPART25D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = -A(SC,TC,UC,N)*(-2*dlog(ptm)*B_QQ()+2*dlog(M)*B_QQ()+FCQQ(1.D
     #0))
	VPART25D(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = C(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
	VPART25D(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = -A(UC,TC,SC,N)*(-2*dlog(ptm)*B_QQ()+2*dlog(M)*B_QQ()+FCQQ(1.D
     #0))
	VPART25D(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = C(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
	VPART25D(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A(UC,SC,TC,N)*(2*B_QQ()*dlog(ptm)-2*B_QQ()*dlog(M)-FCQQ(1.D0)
     #)
	VPART25D(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = 0
	VPART25D(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = -2*CF*A(SC,TC,UC,N)*FCQG(1.D0)
	VPART25D(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = -2*CF*A(UC,TC,SC,N)*FCQG(1.D0)
	VPART25D(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = -2*CF*A(UC,SC,TC,N)*FCQG(1.D0)
	VPART25D(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = (A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))*(-FCQQ(1.D0)-2*B_
     #QQ()*dlog(M)+2*B_QQ()*dlog(ptm))
	VPART25D(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = C(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
	VPART25D(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(-2*B_QQ()*dlog(ptm)+2*B_QQ()*dlog(M)+FCQQ(1.D0))*(A(UC,TC,S
     #C,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART25D(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = C(TC,SC,UC,N)*FCGQ(1.D0)/CF/2
	VPART25D(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = -2*CF*FCQG(1.D0)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART25D(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -2*CF*FCQG(1.D0)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART25D(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-2*B_QQ()*dlog(ptm)+2*B_QQ()*dlog(M)+FCQQ(1.D
     #0))
	VPART25D(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = -C(TC,SC,UC,N)*(2*B_GG()*dlog(ptm)-2*B_GG()*dlog(M)-FCGG(1.D0
     #))
	VPART25D(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = -C(SC,TC,UC,N)*FCQG(1.D0)*(N-1)*(N+1)/N
	VPART25D(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = C(SC,TC,UC,N)*(2*B_GG()*dlog(ptm)-2*B_GG()*dlog(M)-FCGG(1.D0)
     #)
	VPART25D(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = C(UC,SC,TC,N)*FCQG(1.D0)*(N-1)*(N+1)/N
	VPART25D(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*(-2*B_GG()*dlog(ptm)+2*B_GG()*dlog(M)+FCGG(1.
     #D0))
	VPART25D(21) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 3
C*****************************************************
	SUBROUTINE VEC_35RT(SC,TC,UC,Z3,R,VPART35R)
C PARTIE EN LOG(R)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART35R(J0MAX)
	DO I = 1,J0MAX
	  VPART35R(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A4_QQ(z3)*dlog(R**2)/(1-z3)*A(SC,TC,UC,N)
	VPART35R(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = 0
	VPART35R(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = A4_QQ(z3)*dlog(R**2)/(1-z3)*A(UC,TC,SC,N)
	VPART35R(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = 0
	VPART35R(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A4_QQ(z3)*dlog(R**2)/(1-z3)*A(UC,SC,TC,N)
	VPART35R(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = A4_QG(z3)*dlog(R**2)/(1-z3)*C(SC,TC,UC,N)
	VPART35R(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 0
	VPART35R(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 0
	VPART35R(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART35R(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = A4_QQ(z3)*dlog(R**2)/(1-z3)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC
     #,TC,UC,N))
	VPART35R(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = 0
	VPART35R(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = A4_QQ(z3)*dlog(R**2)/(1-z3)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC
     #,TC,SC,N))
	VPART35R(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = A4_QG(z3)*dlog(R**2)/(1-z3)*C(SC,TC,UC,N)
	VPART35R(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = -A4_QG(z3)*dlog(R**2)/(1-z3)*C(UC,SC,TC,N)
	VPART35R(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = 0
	VPART35R(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = A4_GG(z3)*dlog(R**2)/(1-z3)*C(SC,TC,UC,N)
	VPART35R(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = -A4_QQ(z3)*dlog(R**2)/(1-z3)*C(TC,SC,UC,N)
	VPART35R(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = -A4_GQ(z3)*dlog(R**2)/(1-z3)*C(TC,SC,UC,N)
	VPART35R(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = A4_QQ(z3)*dlog(R**2)/(1-z3)*C(SC,TC,UC,N)
	VPART35R(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = A4_QG(z3)*dlog(R**2)/(1-z3)*Dg(SC,TC,UC,N)
	VPART35R(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = A4_GG(z3)*dlog(R**2)/(1-z3)*Dg(SC,TC,UC,N)
	VPART35R(21) = T0
	RETURN
	END
	SUBROUTINE VEC_35ZT(SC,TC,UC,Z3,MF,PT3,VPART35Z)
C PARTIE EN 1/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART35Z(J0MAX)
	DO I = 1,J0MAX
	  VPART35Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A(SC,TC,UC,N)*(ANM4_QQ(z3)+2*A4_QQ(z3)*dlog(pt3)-2*A4_QQ(z3)*
     #dlog(Mf)-DAQQ(z3))
	VPART35Z(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = 0
	VPART35Z(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = A(UC,TC,SC,N)*(ANM4_QQ(z3)+2*A4_QQ(z3)*dlog(pt3)-2*A4_QQ(z3)*
     #dlog(Mf)-DAQQ(z3))
	VPART35Z(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = 0
	VPART35Z(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A(UC,SC,TC,N)*(ANM4_QQ(z3)+2*A4_QQ(z3)*dlog(pt3)-2*A4_QQ(z3)*
     #dlog(Mf)-DAQQ(z3))
	VPART35Z(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = -C(SC,TC,UC,N)*(-ANM4_QG(z3)-2*dlog(pt3)*A4_QG(z3)+2*dlog(Mf)
     #*A4_QG(z3)+DAQG(z3))
	VPART35Z(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 0
	VPART35Z(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 0
	VPART35Z(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART35Z(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = (-DAQQ(z3)+ANM4_QQ(z3)-2*A4_QQ(z3)*dlog(Mf)+2*A4_QQ(z3)*dlog(
     #pt3))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART35Z(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = 0
	VPART35Z(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(-2*A4_QQ(z3)*dlog(pt3)+2*A4_QQ(z3)*dlog(Mf)-ANM4_QQ(z3)+DAQ
     #Q(z3))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART35Z(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = C(SC,TC,UC,N)*(ANM4_QG(z3)+2*dlog(pt3)*A4_QG(z3)-2*dlog(Mf)*A
     #4_QG(z3)-DAQG(z3))
	VPART35Z(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = C(UC,SC,TC,N)*(-ANM4_QG(z3)-2*dlog(pt3)*A4_QG(z3)+2*dlog(Mf)*
     #A4_QG(z3)+DAQG(z3))
	VPART35Z(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = 0
	VPART35Z(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-ANM4_GG(z3)-2*dlog(pt3)*A4_GG(z3)+2*dlog(Mf)
     #*A4_GG(z3)+DAGG(z3))
	VPART35Z(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = -C(TC,SC,UC,N)*(ANM4_QQ(z3)+2*A4_QQ(z3)*dlog(pt3)-2*A4_QQ(z3)
     #*dlog(Mf)-DAQQ(z3))
	VPART35Z(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = C(TC,SC,UC,N)*(-ANM4_GQ(z3)-2*A4_GQ(z3)*dlog(pt3)+2*A4_GQ(z3)
     #*dlog(Mf)+DAGQ(z3))
	VPART35Z(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = -C(SC,TC,UC,N)*(-ANM4_QQ(z3)-2*dlog(pt3)*A4_QQ(z3)+2*dlog(Mf)
     #*A4_QQ(z3)+DAQQ(z3))
	VPART35Z(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = Dg(SC,TC,UC,N)*(ANM4_QG(z3)+2*A4_QG(z3)*dlog(pt3)-2*A4_QG(z3)
     #*dlog(Mf)-DAQG(z3))
	VPART35Z(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*(-ANM4_GG(z3)-2*A4_GG(z3)*dlog(pt3)+2*A4_GG(z
     #3)*dlog(Mf)+DAGG(z3))
	VPART35Z(21) = T0
	RETURN
	END
	SUBROUTINE VEC_35LT(SC,TC,UC,Z3,VPART35L)
C PARTIE EN LOG(1-Z3)/(1-Z3)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART35L(J0MAX)
	DO I = 1,J0MAX
	  VPART35L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A(SC,TC,UC,N)*(-DBQQ(z3)+2*A4_QQ(z3))
	VPART35L(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = 0
	VPART35L(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = A(UC,TC,SC,N)*(-DBQQ(z3)+2*A4_QQ(z3))
	VPART35L(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = 0
	VPART35L(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A(UC,SC,TC,N)*(-DBQQ(z3)+2*A4_QQ(z3))
	VPART35L(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = -C(SC,TC,UC,N)*(DBQG(z3)-2*A4_QG(z3))
	VPART35L(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 0
	VPART35L(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 0
	VPART35L(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART35L(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = -(-2*A4_QQ(z3)+DBQQ(z3))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC
     #,UC,N))
	VPART35L(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = 0
	VPART35L(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(DBQQ(z3)-2*A4_QQ(z3))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,
     #SC,N))
	VPART35L(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = -C(SC,TC,UC,N)*(DBQG(z3)-2*A4_QG(z3))
	VPART35L(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = -C(UC,SC,TC,N)*(-DBQG(z3)+2*A4_QG(z3))
	VPART35L(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = 0
	VPART35L(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(DBGG(z3)-2*A4_GG(z3))
	VPART35L(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = -C(TC,SC,UC,N)*(-DBQQ(z3)+2*A4_QQ(z3))
	VPART35L(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = C(TC,SC,UC,N)*(DBGQ(z3)-2*A4_GQ(z3))
	VPART35L(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = -C(SC,TC,UC,N)*(DBQQ(z3)-2*A4_QQ(z3))
	VPART35L(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = -Dg(SC,TC,UC,N)*(DBQG(z3)-2*A4_QG(z3))
	VPART35L(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*(DBGG(z3)-2*A4_GG(z3))
	VPART35L(21) = T0
	RETURN
	END
	SUBROUTINE VEC_35DT(SC,TC,UC,MF,PT3,VPART35D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART35D(J0MAX)
	DO I = 1,J0MAX
	  VPART35D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = -A(SC,TC,UC,N)*(-2*dlog(pt3)*B_QQ()+2*dlog(Mf)*B_QQ()+DCQQ(1.
     #D0))
	VPART35D(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = 0
	VPART35D(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = -A(UC,TC,SC,N)*(-2*dlog(pt3)*B_QQ()+2*dlog(Mf)*B_QQ()+DCQQ(1.
     #D0))
	VPART35D(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = 0
	VPART35D(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A(UC,SC,TC,N)*(2*B_QQ()*dlog(pt3)-2*B_QQ()*dlog(Mf)-DCQQ(1.D0
     #))
	VPART35D(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = -C(SC,TC,UC,N)*DCQG(1.D0)
	VPART35D(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = 0
	VPART35D(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = 0
	VPART35D(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART35D(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = -(DCQQ(1.D0)-2*B_QQ()*dlog(pt3)+2*B_QQ()*dlog(Mf))*(A(SC,TC,U
     #C,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART35D(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = 0
	VPART35D(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(2*B_QQ()*dlog(Mf)-2*B_QQ()*dlog(pt3)+DCQQ(1.D0))*(A(UC,TC,S
     #C,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART35D(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = -C(SC,TC,UC,N)*DCQG(1.D0)
	VPART35D(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = C(UC,SC,TC,N)*DCQG(1.D0)
	VPART35D(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = 0
	VPART35D(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-2*B_GG()*dlog(pt3)+2*B_GG()*dlog(Mf)+DCGG(1.
     #D0))
	VPART35D(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = C(TC,SC,UC,N)*(-2*dlog(pt3)*B_QQ()+2*dlog(Mf)*B_QQ()+DCQQ(1.D
     #0))
	VPART35D(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = C(TC,SC,UC,N)*DCGQ(1.D0)
	VPART35D(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = -C(SC,TC,UC,N)*(-2*B_QQ()*dlog(pt3)+2*B_QQ()*dlog(Mf)+DCQQ(1.
     #D0))
	VPART35D(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = -Dg(SC,TC,UC,N)*DCQG(1.D0)
	VPART35D(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = Dg(SC,TC,UC,N)*(2*B_GG()*dlog(pt3)-2*B_GG()*dlog(Mf)-DCGG(1.D
     #0))
	VPART35D(21) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE 5 COLLINEAIRE A 4
C*****************************************************
	SUBROUTINE VEC_45RT(SC,TC,UC,Z4,R,VPART45R)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART45R(J0MAX)
	DO I = 1,J0MAX
	  VPART45R(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A4_QQ(z4)*dlog(R**2)/(1-z4)*A(SC,TC,UC,N)
	VPART45R(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = A4_GQ(z4)*dlog(R**2)/(1-z4)*A(SC,TC,UC,N)
	VPART45R(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = A4_QQ(z4)*dlog(R**2)/(1-z4)*A(UC,TC,SC,N)
	VPART45R(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = A4_GQ(z4)*dlog(R**2)/(1-z4)*A(UC,TC,SC,N)
	VPART45R(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = A4_QQ(z4)*dlog(R**2)/(1-z4)*A(UC,SC,TC,N)
	VPART45R(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = A4_GQ(z4)*dlog(R**2)/(1-z4)*A(UC,SC,TC,N)
	VPART45R(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = -A4_QG(z4)*dlog(R**2)/(1-z4)*C(TC,SC,UC,N)
	VPART45R(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = -A4_QG(z4)*dlog(R**2)/(1-z4)*C(TC,SC,UC,N)
	VPART45R(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART45R(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = A4_QQ(z4)*dlog(R**2)/(1-z4)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC
     #,TC,UC,N))
	VPART45R(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = A4_GQ(z4)*dlog(R**2)/(1-z4)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC
     #,TC,UC,N))
	VPART45R(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = A4_QQ(z4)*dlog(R**2)/(1-z4)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC
     #,TC,SC,N))
	VPART45R(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = A4_GQ(z4)*dlog(R**2)/(1-z4)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC
     #,TC,SC,N))
	VPART45R(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = -A4_QG(z4)*dlog(R**2)/(1-z4)*C(TC,SC,UC,N)
	VPART45R(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -A4_QG(z4)*dlog(R**2)/(1-z4)*C(TC,SC,UC,N)
	VPART45R(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = A4_GG(z4)*dlog(R**2)/(1-z4)*C(SC,TC,UC,N)
	VPART45R(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = -A4_GG(z4)*dlog(R**2)/(1-z4)*C(TC,SC,UC,N)
	VPART45R(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = -A4_GQ(z4)*dlog(R**2)/(1-z4)*C(UC,SC,TC,N)
	VPART45R(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = A4_QQ(z4)*dlog(R**2)/(1-z4)*C(SC,TC,UC,N)
	VPART45R(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = A4_GQ(z4)*dlog(R**2)/(1-z4)*C(SC,TC,UC,N)
	VPART45R(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = A4_GG(z4)*dlog(R**2)/(1-z4)*Dg(SC,TC,UC,N)
	VPART45R(21) = T0
	RETURN
	END
	SUBROUTINE VEC_45ZT(SC,TC,UC,Z4,MF,PT4,VPART45Z)
C PARTIE EN 1/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART45Z(J0MAX)
	DO I = 1,J0MAX
	  VPART45Z(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A(SC,TC,UC,N)*(ANM4_QQ(z4)+2*A4_QQ(z4)*dlog(pt4)-2*A4_QQ(z4)*
     #dlog(Mf)-DAQQ(z4))
	VPART45Z(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = A(SC,TC,UC,N)*(ANM4_GQ(z4)+2*A4_GQ(z4)*dlog(pt4)-2*A4_GQ(z4)*
     #dlog(Mf)-DAGQ(z4))
	VPART45Z(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = -A(UC,TC,SC,N)*(-ANM4_QQ(z4)-2*A4_QQ(z4)*dlog(pt4)+2*A4_QQ(z4
     #)*dlog(Mf)+DAQQ(z4))
	VPART45Z(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = A(UC,TC,SC,N)*(ANM4_GQ(z4)+2*A4_GQ(z4)*dlog(pt4)-2*A4_GQ(z4)*
     #dlog(Mf)-DAGQ(z4))
	VPART45Z(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = -A(UC,SC,TC,N)*(-ANM4_QQ(z4)-2*A4_QQ(z4)*dlog(pt4)+2*A4_QQ(z4
     #)*dlog(Mf)+DAQQ(z4))
	VPART45Z(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = -A(UC,SC,TC,N)*(-ANM4_GQ(z4)-2*dlog(pt4)*A4_GQ(z4)+2*dlog(Mf)
     #*A4_GQ(z4)+DAGQ(z4))
	VPART45Z(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = C(TC,SC,UC,N)*(-ANM4_QG(z4)-2*A4_QG(z4)*dlog(pt4)+2*A4_QG(z4)
     #*dlog(Mf)+DAQG(z4))
	VPART45Z(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = C(TC,SC,UC,N)*(-ANM4_QG(z4)-2*A4_QG(z4)*dlog(pt4)+2*A4_QG(z4)
     #*dlog(Mf)+DAQG(z4))
	VPART45Z(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART45Z(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = -(DAQQ(z4)+2*A4_QQ(z4)*dlog(Mf)-2*A4_QQ(z4)*dlog(pt4)-ANM4_QQ
     #(z4))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART45Z(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = -(-ANM4_GQ(z4)-2*A4_GQ(z4)*dlog(pt4)+2*A4_GQ(z4)*dlog(Mf)+DAG
     #Q(z4))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART45Z(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(DAQQ(z4)+2*A4_QQ(z4)*dlog(Mf)-2*A4_QQ(z4)*dlog(pt4)-ANM4_QQ
     #(z4))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART45Z(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = -(-ANM4_GQ(z4)-2*A4_GQ(z4)*dlog(pt4)+2*A4_GQ(z4)*dlog(Mf)+DAG
     #Q(z4))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART45Z(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = C(TC,SC,UC,N)*(-ANM4_QG(z4)-2*dlog(pt4)*A4_QG(z4)+2*dlog(Mf)*
     #A4_QG(z4)+DAQG(z4))
	VPART45Z(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -C(TC,SC,UC,N)*(ANM4_QG(z4)+2*dlog(pt4)*A4_QG(z4)-2*dlog(Mf)*
     #A4_QG(z4)-DAQG(z4))
	VPART45Z(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-ANM4_GG(z4)-2*A4_GG(z4)*dlog(pt4)+2*A4_GG(z4
     #)*dlog(Mf)+DAGG(z4))
	VPART45Z(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = C(TC,SC,UC,N)*(-ANM4_GG(z4)-2*A4_GG(z4)*dlog(pt4)+2*A4_GG(z4)
     #*dlog(Mf)+DAGG(z4))
	VPART45Z(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = C(UC,SC,TC,N)*(-ANM4_GQ(z4)-2*A4_GQ(z4)*dlog(pt4)+2*A4_GQ(z4)
     #*dlog(Mf)+DAGQ(z4))
	VPART45Z(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = -C(SC,TC,UC,N)*(-ANM4_QQ(z4)-2*A4_QQ(z4)*dlog(pt4)+2*A4_QQ(z4
     #)*dlog(Mf)+DAQQ(z4))
	VPART45Z(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = -C(SC,TC,UC,N)*(-ANM4_GQ(z4)-2*A4_GQ(z4)*dlog(pt4)+2*A4_GQ(z4
     #)*dlog(Mf)+DAGQ(z4))
	VPART45Z(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = Dg(SC,TC,UC,N)*(ANM4_GG(z4)+2*dlog(pt4)*A4_GG(z4)-2*dlog(Mf)*
     #A4_GG(z4)-DAGG(z4))
	VPART45Z(21) = T0
	RETURN
	END
	SUBROUTINE VEC_45LT(SC,TC,UC,Z4,VPART45L)
C PARTIE EN LOG(1-Z4)/(1-Z4)+
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART45L(J0MAX)
	DO I = 1,J0MAX
	  VPART45L(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A(SC,TC,UC,N)*(-DBQQ(z4)+2*A4_QQ(z4))
	VPART45L(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = A(SC,TC,UC,N)*(-DBGQ(z4)+2*A4_GQ(z4))
	VPART45L(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = A(UC,TC,SC,N)*(-DBQQ(z4)+2*A4_QQ(z4))
	VPART45L(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = A(UC,TC,SC,N)*(-DBGQ(z4)+2*A4_GQ(z4))
	VPART45L(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = -A(UC,SC,TC,N)*(DBQQ(z4)-2*A4_QQ(z4))
	VPART45L(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = -A(UC,SC,TC,N)*(DBGQ(z4)-2*A4_GQ(z4))
	VPART45L(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = C(TC,SC,UC,N)*(DBQG(z4)-2*A4_QG(z4))
	VPART45L(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = C(TC,SC,UC,N)*(DBQG(z4)-2*A4_QG(z4))
	VPART45L(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART45L(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = -(-2*A4_QQ(z4)+DBQQ(z4))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC
     #,UC,N))
	VPART45L(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = -(DBGQ(z4)-2*A4_GQ(z4))*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,
     #UC,N))
	VPART45L(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(DBQQ(z4)-2*A4_QQ(z4))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,
     #SC,N))
	VPART45L(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = -(DBGQ(z4)-2*A4_GQ(z4))*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,
     #SC,N))
	VPART45L(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = C(TC,SC,UC,N)*(DBQG(z4)-2*A4_QG(z4))
	VPART45L(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = -C(TC,SC,UC,N)*(-DBQG(z4)+2*A4_QG(z4))
	VPART45L(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(DBGG(z4)-2*A4_GG(z4))
	VPART45L(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = C(TC,SC,UC,N)*(DBGG(z4)-2*A4_GG(z4))
	VPART45L(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = C(UC,SC,TC,N)*(DBGQ(z4)-2*A4_GQ(z4))
	VPART45L(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = C(SC,TC,UC,N)*(-DBQQ(z4)+2*A4_QQ(z4))
	VPART45L(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = -C(SC,TC,UC,N)*(DBGQ(z4)-2*A4_GQ(z4))
	VPART45L(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = -Dg(SC,TC,UC,N)*(DBGG(z4)-2*A4_GG(z4))
	VPART45L(21) = T0
	RETURN
	END
	SUBROUTINE VEC_45DT(SC,TC,UC,MF,PT4,VPART45D)
C PARTIE EN DELTA(1-Z2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPART45D(J0MAX)
	DO I = 1,J0MAX
	  VPART45D(I) = 0.D0
	ENDDO
C	j0 = 1 : qi + qk --> qi + qk
      t0 = A(SC,TC,UC,N)*(2*B_QQ()*dlog(pt4)-2*dlog(Mf)*B_QQ()-DCQQ(1.D0
     #))
	VPART45D(1) = T0
C	j0 = 2 : qi + qk --> qi + g
      t0 = -A(SC,TC,UC,N)*DCGQ(1.D0)
	VPART45D(2) = T0
C	j0 = 3 : qi + qbk --> qi + qbk
      t0 = -A(UC,TC,SC,N)*(-2*B_QQ()*dlog(pt4)+2*dlog(Mf)*B_QQ()+DCQQ(1.
     #D0))
	VPART45D(3) = T0
C	j0 = 4 : qi + qbk --> qi + g
      t0 = -A(UC,TC,SC,N)*DCGQ(1.D0)
	VPART45D(4) = T0
C	j0 = 5 : qi + qbi --> qk + qbk
      t0 = -A(UC,SC,TC,N)*(-2*dlog(pt4)*B_QQ()+2*B_QQ()*dlog(Mf)+DCQQ(1.
     #D0))
	VPART45D(5) = T0
C	j0 = 6 : qi + qbi --> qk + g
      t0 = -A(UC,SC,TC,N)*DCGQ(1.D0)
	VPART45D(6) = T0
C	j0 = 7 : qi + g --> qi + qk
      t0 = C(TC,SC,UC,N)*DCQG(1.D0)
	VPART45D(7) = T0
C	j0 = 8 : qi + g --> qi + qbk
      t0 = C(TC,SC,UC,N)*DCQG(1.D0)
	VPART45D(8) = T0
C	j0 = 9 : qi + g --> qk + qbk
      t0 = 0
	VPART45D(9) = T0
C	j0 = 10 : qi + qi --> qi + qi
      t0 = -(DCQQ(1.D0)+2*B_QQ()*dlog(Mf)-2*B_QQ()*dlog(pt4))*(A(SC,TC,U
     #C,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART45D(10) = T0
C	j0 = 11 : qi + qi --> qi + g
      t0 = -DCGQ(1.D0)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))
	VPART45D(11) = T0
C	j0 = 12 : qi + qbi --> qi + qbi
      t0 = -(2*B_QQ()*dlog(Mf)+DCQQ(1.D0)-2*B_QQ()*dlog(pt4))*(A(UC,TC,S
     #C,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART45D(12) = T0
C	j0 = 13 : qi + qbi --> qi + g
      t0 = -DCGQ(1.D0)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))
	VPART45D(13) = T0
C	j0 = 14 : qi + g --> qi + qi
      t0 = C(TC,SC,UC,N)*DCQG(1.D0)
	VPART45D(14) = T0
C	j0 = 15 : qi + g --> qi + qbi
      t0 = C(TC,SC,UC,N)*DCQG(1.D0)
	VPART45D(15) = T0
C	j0 = 16 : qi + qbi --> g + g
      t0 = -C(SC,TC,UC,N)*(-2*B_GG()*dlog(pt4)+2*B_GG()*dlog(Mf)+DCGG(1.
     #D0))
	VPART45D(16) = T0
C	j0 = 17 : qi + g --> qi + g
      t0 = C(TC,SC,UC,N)*(-2*B_GG()*dlog(pt4)+2*B_GG()*dlog(Mf)+DCGG(1.D
     #0))
	VPART45D(17) = T0
C	j0 = 18 : qi + g --> g + g
      t0 = C(UC,SC,TC,N)*DCGQ(1.D0)
	VPART45D(18) = T0
C	j0 = 19 : g + g --> qi +qbi
      t0 = -C(SC,TC,UC,N)*(-2*B_QQ()*dlog(pt4)+2*B_QQ()*dlog(Mf)+DCQQ(1.
     #D0))
	VPART45D(19) = T0
C	j0 = 20 : g + g --> qi + g
      t0 = -C(SC,TC,UC,N)*DCGQ(1.D0)
	VPART45D(20) = T0
C	j0 = 21 : g + g --> g + g
      t0 = Dg(SC,TC,UC,N)*(2*B_GG()*dlog(pt4)-2*B_GG()*dlog(Mf)-DCGG(1.D
     #0))
	VPART45D(21) = T0
	RETURN
	END
C*****************************************************
CC  TERMES QUI VIENNENT DE LA PARTIE VIRTUELLE
C*****************************************************
	SUBROUTINE VEC_VIT(S,SC,TC,UC,MU,PT,PTM,VPARTVI)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
	DIMENSION VPARTVI(J0MAX)
	PI=4.D0*DATAN(1.D0)
	DO I = 1,J0MAX
	  VPARTVI(I) = 0.D0
	ENDDO
	MU2=MU**2
C	j0 = 1 : qi + qk --> qi + qk
C termes provenant des parties collineaires
      teraj = -2*B_QQ()*dlog(ptm**2/S)*A(SC,TC,UC,N)-2*B_QQ()*dlog(pt**2
     #/S)*A(SC,TC,UC,N)+(-CF*(8*dlog(SC/S)-8*dlog(-UC/S)-4*dlog(-TC/S))+
     #N*(4*dlog(SC/S)-2*dlog(-UC/S)-2*dlog(-TC/S)))*A(SC,TC,UC,N)*dlog(p
     #tm**2/S)+2*A4_QQ(1.D0)*A(SC,TC,UC,N)*(dlog(pt**2/S)**2/4-dlog(pt**
     #2/S)*dlog(ptm**2/S)/2)-A4_QQ(1.D0)*A(SC,TC,UC,N)*dlog(ptm**2/S)**2
     #/2
C termes finies d'Ellis et Sexton
      s2 = CF*(-16-2*dlog(-TC/S)**2+dlog(-TC/S)*(6+8*dlog(SC/S)-8*dlog(-
     #UC/S))+(-2*SC**2+2*UC**2)/(SC**2+UC**2)*(Pi**2+dlog(-TC/S)**2-2*dl
     #og(-TC/S)*dlog(SC/S)+dlog(SC/S)**2+(dlog(-TC/S)-dlog(-UC/S))**2)+(
     #2*SC+2*UC)/(SC**2+UC**2)*((SC+UC)*(dlog(-UC/S)-dlog(SC/S))+(UC-SC)
     #*(2*dlog(-TC/S)-dlog(SC/S)-dlog(-UC/S))))
      s3 = N*(85.D0/9.D0+Pi**2+2*dlog(-TC/S)*(dlog(-TC/S)+dlog(-UC/S)-2*
     #dlog(SC/S))+(SC**2-UC**2)/(2*SC**2+2*UC**2)*(Pi**2+2*dlog(-TC/S)**
     #2-4*dlog(-TC/S)*dlog(SC/S)+2*dlog(SC/S)**2+(dlog(-TC/S)-dlog(-UC/S
     #))**2)-SC*TC/(SC**2+UC**2)*(dlog(-TC/S)-dlog(-UC/S))+2*UC*TC/(SC**
     #2+UC**2)*(dlog(-TC/S)-dlog(SC/S))+11.D0/3.D0*dlog(MU2/S)-11.D0/3.D
     #0*dlog(-TC/S))+GTR*(4.D0/3.D0*dlog(-TC/S)-4.D0/3.D0*dlog(MU2/S)-20
     #.D0/9.D0)
      s1 = s2+s3
      s2 = A(SC,TC,UC,N)
      teres = s1*s2
	VPARTVI(1)=(TERAJ+TERES)
C	j0 = 2 : qi + qk --> qi + g
      teraj = 0
      teres = 0
	VPARTVI(2)=(TERAJ+TERES)
C	j0 = 3 : qi + qbk --> qi + qbk
      teraj = -2*B_QQ()*dlog(ptm**2/S)*A(UC,TC,SC,N)-2*B_QQ()*dlog(pt**2
     #/S)*A(UC,TC,SC,N)+(-CF*(8*dlog(-UC/S)-8*dlog(SC/S)-4*dlog(-TC/S))+
     #N*(4*dlog(-UC/S)-2*dlog(SC/S)-2*dlog(-TC/S)))*A(UC,TC,SC,N)*dlog(p
     #tm**2/S)+2*A4_QQ(1.D0)*A(UC,TC,SC,N)*(dlog(pt**2/S)**2/4-dlog(pt**
     #2/S)*dlog(ptm**2/S)/2)-A4_QQ(1.D0)*A(UC,TC,SC,N)*dlog(ptm**2/S)**2
     #/2
      s2 = CF*(-16-2*dlog(-TC/S)**2+dlog(-TC/S)*(6+8*dlog(-UC/S)-8*dlog(
     #SC/S))+(-2*UC**2+2*SC**2)/(UC**2+SC**2)*(Pi**2+(dlog(-TC/S)-dlog(-
     #UC/S))**2+dlog(-TC/S)**2-2*dlog(-TC/S)*dlog(SC/S)+dlog(SC/S)**2)+(
     #2*UC+2*SC)/(UC**2+SC**2)*((UC+SC)*(dlog(SC/S)-dlog(-UC/S))+(SC-UC)
     #*(2*dlog(-TC/S)-dlog(-UC/S)-dlog(SC/S))))
      s3 = N*(85.D0/9.D0+Pi**2+2*dlog(-TC/S)*(dlog(-TC/S)+dlog(SC/S)-2*d
     #log(-UC/S))+(UC**2-SC**2)/(2*UC**2+2*SC**2)*(2*Pi**2+2*(dlog(-TC/S
     #)-dlog(-UC/S))**2+dlog(-TC/S)**2-2*dlog(-TC/S)*dlog(SC/S)+dlog(SC/
     #S)**2)-UC*TC/(UC**2+SC**2)*(dlog(-TC/S)-dlog(SC/S))+2*SC*TC/(UC**2
     #+SC**2)*(dlog(-TC/S)-dlog(-UC/S))+11.D0/3.D0*dlog(MU2/S)-11.D0/3.D
     #0*dlog(-TC/S))+GTR*(4.D0/3.D0*dlog(-TC/S)-4.D0/3.D0*dlog(MU2/S)-20
     #.D0/9.D0)
      s1 = s2+s3
      s2 = A(UC,TC,SC,N)
      teres = s1*s2
	VPARTVI(3)=(TERAJ+TERES)
C	j0 = 4 : qi + qbk --> qi + g
      teraj = 0
      teres = 0
	VPARTVI(4)=(TERAJ+TERES)
C	j0 = 5 : qi + qbi --> qk + qbk
      teraj = -2*B_QQ()*dlog(ptm**2/S)*A(UC,SC,TC,N)-2*B_QQ()*dlog(pt**2
     #/S)*A(UC,SC,TC,N)+(-CF*(8*dlog(-UC/S)-8*dlog(-TC/S)-4*dlog(SC/S))+
     #N*(4*dlog(-UC/S)-2*dlog(-TC/S)-2*dlog(SC/S)))*A(UC,SC,TC,N)*dlog(p
     #tm**2/S)+2*A4_QQ(1.D0)*A(UC,SC,TC,N)*(dlog(pt**2/S)**2/4-dlog(pt**
     #2/S)*dlog(ptm**2/S)/2)-A4_QQ(1.D0)*A(UC,SC,TC,N)*dlog(ptm**2/S)**2
     #/2
      s2 = CF*(-16+2*Pi**2-2*dlog(SC/S)**2+dlog(SC/S)*(6+8*dlog(-UC/S)-8
     #*dlog(-TC/S))+(-2*UC**2+2*TC**2)/(UC**2+TC**2)*(2*dlog(SC/S)**2-2*
     #dlog(SC/S)*dlog(-UC/S)+dlog(-UC/S)**2-2*dlog(SC/S)*dlog(-TC/S)+dlo
     #g(-TC/S)**2)+(2*UC+2*TC)/(UC**2+TC**2)*((UC+TC)*(dlog(-TC/S)-dlog(
     #-UC/S))+(TC-UC)*(2*dlog(SC/S)-dlog(-UC/S)-dlog(-TC/S))))
      s3 = N*(85.D0/9.D0-Pi**2+2*dlog(SC/S)*(dlog(SC/S)+dlog(-TC/S)-2*dl
     #og(-UC/S))+(UC**2-TC**2)/(2*UC**2+2*TC**2)*(3*dlog(SC/S)**2-4*dlog
     #(SC/S)*dlog(-UC/S)+2*dlog(-UC/S)**2-2*dlog(SC/S)*dlog(-TC/S)+dlog(
     #-TC/S)**2)-UC*SC/(UC**2+TC**2)*(dlog(SC/S)-dlog(-TC/S))+2*TC*SC/(U
     #C**2+TC**2)*(dlog(SC/S)-dlog(-UC/S))+11.D0/3.D0*dlog(MU2/S)-11.D0/
     #3.D0*dlog(SC/S))+GTR*(4.D0/3.D0*dlog(SC/S)-4.D0/3.D0*dlog(MU2/S)-2
     #0.D0/9.D0)
      s1 = s2+s3
      s2 = A(UC,SC,TC,N)
      teres = s1*s2
	VPARTVI(5)=(TERAJ+TERES)
C	j0 = 6 : qi + qbi --> qk + g
      teraj = 0
      teres = 0
	VPARTVI(6)=(TERAJ+TERES)
C	j0 = 7 : qi + g --> qi + qk
      teraj = 0
      teres = 0
	VPARTVI(7)=(TERAJ+TERES)
C	j0 = 8 : qi + g --> qi + qbk
      teraj = 0
      teres = 0
	VPARTVI(8)=(TERAJ+TERES)
C	j0 = 9 : qi + g --> qk + qbk
      teraj = 0
      teres = 0
	VPARTVI(9)=(TERAJ+TERES)
C	j0 = 10 : qi + qi --> qi + qi
      s1 = -2*B_QQ()*dlog(ptm**2/S)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC
     #,UC,N))-2*B_QQ()*dlog(pt**2/S)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,T
     #C,UC,N))
      teraj = s1+((-CF*(4*dlog(SC/S)-4*dlog(-TC/S)-4*dlog(-UC/S))+2*N*(2
     #*dlog(SC/S)-dlog(-TC/S)-dlog(-UC/S)))*B(SC,TC,UC,N)+(-CF*(8*dlog(S
     #C/S)-8*dlog(-UC/S)-4*dlog(-TC/S))+N*(4*dlog(SC/S)-2*dlog(-UC/S)-2*
     #dlog(-TC/S)))*A(SC,TC,UC,N)+(-CF*(8*dlog(SC/S)-8*dlog(-TC/S)-4*dlo
     #g(-UC/S))+N*(4*dlog(SC/S)-2*dlog(-UC/S)-2*dlog(-TC/S)))*A(SC,UC,TC
     #,N))*dlog(ptm**2/S)+2*A4_QQ(1.D0)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(S
     #C,TC,UC,N))*(dlog(pt**2/S)**2/4-dlog(pt**2/S)*dlog(ptm**2/S)/2)-A4
     #_QQ(1.D0)*(A(SC,TC,UC,N)+A(SC,UC,TC,N)+B(SC,TC,UC,N))*dlog(ptm**2/
     #S)**2/2
      s1 = (CF*(-16-3.D0/2.D0*(dlog(-TC/S)+dlog(-UC/S))**2+2*dlog(SC/S)*
     #(dlog(-TC/S)+dlog(-UC/S))+2*dlog(-TC/S)+2*dlog(-UC/S)-Pi**2/2)+N*(
     #85.D0/9.D0+5.D0/4.D0*(dlog(-TC/S)+dlog(-UC/S))**2-dlog(-TC/S)*dlog
     #(-UC/S)-2*dlog(SC/S)*(dlog(-TC/S)+dlog(-UC/S))-4.D0/3.D0*dlog(-TC/
     #S)-4.D0/3.D0*dlog(-UC/S)+5.D0/4.D0*Pi**2+11.D0/3.D0*dlog(MU2/S))+G
     #TR*(-20.D0/9.D0+2.D0/3.D0*dlog(-TC/S)+2.D0/3.D0*dlog(-UC/S)-4.D0/3
     #.D0*dlog(MU2/S))+1/N*((Pi**2/2+(dlog(-TC/S)-dlog(-UC/S))**2/2)*UC*
     #TC/SC**2+UC/SC*dlog(-TC/S)+TC/SC*dlog(-UC/S)))*B(SC,TC,UC,N)
      s5 = CF*(-16-2*dlog(-TC/S)**2+dlog(-TC/S)*(6+8*dlog(SC/S)-8*dlog(-
     #UC/S))+(-2*SC**2+2*UC**2)/(SC**2+UC**2)*(Pi**2+dlog(-TC/S)**2-2*dl
     #og(-TC/S)*dlog(SC/S)+dlog(SC/S)**2+(dlog(-TC/S)-dlog(-UC/S))**2)+(
     #2*SC+2*UC)/(SC**2+UC**2)*((SC+UC)*(dlog(-UC/S)-dlog(SC/S))+(UC-SC)
     #*(2*dlog(-TC/S)-dlog(SC/S)-dlog(-UC/S))))
      s6 = N*(85.D0/9.D0+Pi**2+2*dlog(-TC/S)*(dlog(-TC/S)+dlog(-UC/S)-2*
     #dlog(SC/S))+(SC**2-UC**2)/(2*SC**2+2*UC**2)*(Pi**2+2*dlog(-TC/S)**
     #2-4*dlog(-TC/S)*dlog(SC/S)+2*dlog(SC/S)**2+(dlog(-TC/S)-dlog(-UC/S
     #))**2)-SC*TC/(SC**2+UC**2)*(dlog(-TC/S)-dlog(-UC/S))+2*UC*TC/(SC**
     #2+UC**2)*(dlog(-TC/S)-dlog(SC/S))+11.D0/3.D0*dlog(MU2/S)-11.D0/3.D
     #0*dlog(-TC/S))+GTR*(4.D0/3.D0*dlog(-TC/S)-4.D0/3.D0*dlog(MU2/S)-20
     #.D0/9.D0)
      s4 = s5+s6
      s5 = A(SC,TC,UC,N)
      s3 = s4*s5
      s6 = CF*(-16-2*dlog(-UC/S)**2+dlog(-UC/S)*(6+8*dlog(SC/S)-8*dlog(-
     #TC/S))+(-2*SC**2+2*TC**2)/(SC**2+TC**2)*(Pi**2+dlog(-UC/S)**2-2*dl
     #og(-UC/S)*dlog(SC/S)+dlog(SC/S)**2+(dlog(-UC/S)-dlog(-TC/S))**2)+(
     #2*SC+2*TC)/(SC**2+TC**2)*((SC+TC)*(dlog(-TC/S)-dlog(SC/S))+(TC-SC)
     #*(2*dlog(-UC/S)-dlog(SC/S)-dlog(-TC/S))))
      s7 = N*(85.D0/9.D0+Pi**2+2*dlog(-UC/S)*(dlog(-TC/S)+dlog(-UC/S)-2*
     #dlog(SC/S))+(SC**2-TC**2)/(2*SC**2+2*TC**2)*(Pi**2+2*dlog(-UC/S)**
     #2-4*dlog(-UC/S)*dlog(SC/S)+2*dlog(SC/S)**2+(dlog(-UC/S)-dlog(-TC/S
     #))**2)-SC*UC/(SC**2+TC**2)*(dlog(-UC/S)-dlog(-TC/S))+2*TC*UC/(SC**
     #2+TC**2)*(dlog(-UC/S)-dlog(SC/S))+11.D0/3.D0*dlog(MU2/S)-11.D0/3.D
     #0*dlog(-UC/S))+GTR*(4.D0/3.D0*dlog(-UC/S)-4.D0/3.D0*dlog(MU2/S)-20
     #.D0/9.D0)
      s5 = s6+s7
      s6 = A(SC,UC,TC,N)
      s4 = s5*s6
      s2 = s3+s4
      teres = s1+s2
      	VPARTVI(10)=(TERAJ+TERES)
C	j0 = 11 : qi + qi --> qi + g
      teraj = 0
      teres = 0
      	VPARTVI(11)=(TERAJ+TERES)
C	j0 = 12 : qi + qbi --> qi + qbi
      s1 = -2*B_QQ()*dlog(ptm**2/S)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC
     #,SC,N))-2*B_QQ()*dlog(pt**2/S)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,T
     #C,SC,N))
      teraj = s1+((-CF*(4*dlog(-UC/S)-4*dlog(-TC/S)-4*dlog(SC/S))+2*N*(2
     #*dlog(-UC/S)-dlog(-TC/S)-dlog(SC/S)))*B(UC,TC,SC,N)+(-CF*(8*dlog(-
     #UC/S)-8*dlog(SC/S)-4*dlog(-TC/S))+N*(4*dlog(-UC/S)-2*dlog(SC/S)-2*
     #dlog(-TC/S)))*A(UC,TC,SC,N)+(-CF*(8*dlog(-UC/S)-8*dlog(-TC/S)-4*dl
     #og(SC/S))+N*(4*dlog(-UC/S)-2*dlog(SC/S)-2*dlog(-TC/S)))*A(UC,SC,TC
     #,N))*dlog(ptm**2/S)+2*A4_QQ(1.D0)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(U
     #C,TC,SC,N))*(dlog(pt**2/S)**2/4-dlog(pt**2/S)*dlog(ptm**2/S)/2)-A4
     #_QQ(1.D0)*(A(UC,TC,SC,N)+A(UC,SC,TC,N)+B(UC,TC,SC,N))*dlog(ptm**2/
     #S)**2/2
      s1 = (CF*(-16-3.D0/2.D0*dlog(-TC/S)**2-3*dlog(-TC/S)*dlog(SC/S)+Pi
     #**2-3.D0/2.D0*dlog(SC/S)**2+2*dlog(-UC/S)*(dlog(-TC/S)+dlog(SC/S))
     #+2*dlog(-TC/S)+2*dlog(SC/S))+N*(85.D0/9.D0+5.D0/4.D0*dlog(-TC/S)**
     #2+3.D0/2.D0*dlog(-TC/S)*dlog(SC/S)+5.D0/4.D0*dlog(SC/S)**2-2*dlog(
     #-UC/S)*(dlog(-TC/S)+dlog(SC/S))-4.D0/3.D0*dlog(-TC/S)-4.D0/3.D0*dl
     #og(SC/S)+11.D0/3.D0*dlog(MU2/S))+GTR*(-20.D0/9.D0+2.D0/3.D0*dlog(-
     #TC/S)+2.D0/3.D0*dlog(SC/S)-4.D0/3.D0*dlog(MU2/S))+1/N*(-(-dlog(-TC
     #/S)**2/2+dlog(-TC/S)*dlog(SC/S)-dlog(SC/S)**2/2)*SC*TC/UC**2+SC/UC
     #*dlog(-TC/S)+TC/UC*dlog(SC/S)))*B(UC,TC,SC,N)
      s5 = CF*(-16-2*dlog(-TC/S)**2+dlog(-TC/S)*(6+8*dlog(-UC/S)-8*dlog(
     #SC/S))+(-2*UC**2+2*SC**2)/(UC**2+SC**2)*(Pi**2+(dlog(-TC/S)-dlog(-
     #UC/S))**2+dlog(-TC/S)**2-2*dlog(-TC/S)*dlog(SC/S)+dlog(SC/S)**2)+(
     #2*UC+2*SC)/(UC**2+SC**2)*((UC+SC)*(dlog(SC/S)-dlog(-UC/S))+(SC-UC)
     #*(2*dlog(-TC/S)-dlog(-UC/S)-dlog(SC/S))))
      s6 = N*(85.D0/9.D0+Pi**2+2*dlog(-TC/S)*(dlog(-TC/S)+dlog(SC/S)-2*d
     #log(-UC/S))+(UC**2-SC**2)/(2*UC**2+2*SC**2)*(2*Pi**2+2*(dlog(-TC/S
     #)-dlog(-UC/S))**2+dlog(-TC/S)**2-2*dlog(-TC/S)*dlog(SC/S)+dlog(SC/
     #S)**2)-UC*TC/(UC**2+SC**2)*(dlog(-TC/S)-dlog(SC/S))+2*SC*TC/(UC**2
     #+SC**2)*(dlog(-TC/S)-dlog(-UC/S))+11.D0/3.D0*dlog(MU2/S)-11.D0/3.D
     #0*dlog(-TC/S))+GTR*(4.D0/3.D0*dlog(-TC/S)-4.D0/3.D0*dlog(MU2/S)-20
     #.D0/9.D0)
      s4 = s5+s6
      s5 = A(UC,TC,SC,N)
      s3 = s4*s5
      s6 = CF*(-16+2*Pi**2-2*dlog(SC/S)**2+dlog(SC/S)*(6+8*dlog(-UC/S)-8
     #*dlog(-TC/S))+(-2*UC**2+2*TC**2)/(UC**2+TC**2)*(2*dlog(SC/S)**2-2*
     #dlog(SC/S)*dlog(-UC/S)+dlog(-UC/S)**2-2*dlog(-TC/S)*dlog(SC/S)+dlo
     #g(-TC/S)**2)+(2*UC+2*TC)/(UC**2+TC**2)*((UC+TC)*(dlog(-TC/S)-dlog(
     #-UC/S))+(TC-UC)*(2*dlog(SC/S)-dlog(-UC/S)-dlog(-TC/S))))
      s7 = N*(85.D0/9.D0-Pi**2+2*dlog(SC/S)*(dlog(-TC/S)+dlog(SC/S)-2*dl
     #og(-UC/S))+(UC**2-TC**2)/(2*UC**2+2*TC**2)*(3*dlog(SC/S)**2-4*dlog
     #(SC/S)*dlog(-UC/S)+2*dlog(-UC/S)**2-2*dlog(-TC/S)*dlog(SC/S)+dlog(
     #-TC/S)**2)-UC*SC/(UC**2+TC**2)*(dlog(SC/S)-dlog(-TC/S))+2*TC*SC/(U
     #C**2+TC**2)*(dlog(SC/S)-dlog(-UC/S))+11.D0/3.D0*dlog(MU2/S)-11.D0/
     #3.D0*dlog(SC/S))+GTR*(4.D0/3.D0*dlog(SC/S)-4.D0/3.D0*dlog(MU2/S)-2
     #0.D0/9.D0)
      s5 = s6+s7
      s6 = A(UC,SC,TC,N)
      s4 = s5*s6
      s2 = s3+s4
      teres = s1+s2
      	VPARTVI(12)=(TERAJ+TERES)
C	j0 = 13 : qi + qbi --> qi + g
      teraj = 0
      teres = 0
      	VPARTVI(13)=(TERAJ+TERES)
C	j0 = 14 : qi + g --> qi + qi
      teraj = 0
      teres = 0
      	VPARTVI(14)=(TERAJ+TERES)
C	j0 = 15 : qi + g --> qi + qbi
      teraj = 0
      teres = 0
      	VPARTVI(15)=(TERAJ+TERES)
C	j0 = 16 : qi + qbi --> g + g
      s1 = -2*B_QQ()*dlog(ptm**2/S)*C(SC,TC,UC,N)-2*B_GG()*dlog(pt**2/S)
     #*C(SC,TC,UC,N)+(dlog(SC/S)*((2*N**2*(N**2-1)+2*(N**2-1)/N**2)*(TC*
     #*2+UC**2)/TC/UC-4*(N**2-1)**2*(TC**2+UC**2)/SC**2)+4*N**2*(N**2-1)
     #*(dlog(-TC/S)*(UC/TC-2/SC**2*UC**2)+dlog(-UC/S)*(1/UC*TC-2/SC**2*T
     #C**2))-4*(N**2-1)*(UC/TC+1/UC*TC)*(dlog(-TC/S)+dlog(-UC/S)))*dlog(
     #ptm**2/S)
      teraj = s1-32*(SC**2/4+SC*TC/4+TC**2/4)*(N-1)*(N+1)*(dlog(SC/S)*SC
     #*N**4*TC/2+dlog(SC/S)*N**4*SC**2/4+dlog(SC/S)*N**4*TC**2/2-dlog(SC
     #/S)*SC*TC*N**2/2-dlog(SC/S)*TC**2*N**2/2+dlog(SC/S)*SC**2/4+dlog(-
     #TC/S)*N**4*TC**2/2+dlog(-TC/S)*N**4*SC*TC+dlog(-TC/S)*N**4*SC**2/2
     #-dlog(-TC/S)*N**2*SC**2/2+dlog(-UC/S)*N**4*TC**2/2-dlog(-UC/S)*N**
     #2*SC**2/2)/SC**2/TC/N**2/(SC/2+TC/2)+2*A4_GG(1.D0)*C(SC,TC,UC,N)*(
     #dlog(pt**2/S)**2/4-dlog(pt**2/S)*dlog(ptm**2/S)/2)-A4_QQ(1.D0)*C(S
     #C,TC,UC,N)*dlog(ptm**2/S)**2/2
      s1 = (-7*N+7/N+22.D0/3.D0*N*dlog(MU2/S)-8.D0/3.D0*GTR*dlog(MU2/S))
     #*(N**2-1)/N*(1/UC/TC*N**2-1/UC/TC-2*N**2/SC**2)*(TC**2+UC**2)
      s4 = 4*N
      s6 = N**2-1
      s8 = dlog(-TC/S)*dlog(-UC/S)/N*(TC**2+UC**2)/TC/UC/2+(-Pi**2+dlog(
     #SC/S)**2)*(1/N**3*SC**2/TC/UC/4+1/N/8+1/N*UC/TC/4+1/N*TC/UC/4-1/N*
     #TC**2/SC**2/4-1/N*UC**2/SC**2/4-N/SC**2*TC**2/4-N/SC**2*UC**2/4)+d
     #log(SC/S)*(5.D0/8.D0*N-9.D0/8.D0/N-1/N**3-TC/UC*N/2-1/TC*UC*N/2-TC
     #/UC/N**3/2-1/TC*UC/N**3/2-N/SC**2*TC**2/4-N/SC**2*UC**2/4+1/N*TC**
     #2/SC**2/4+1/N*UC**2/SC**2/4)+Pi**2/N/8+3.D0/8.D0*Pi**2/TC*UC/N**3+
     #3.D0/8.D0*Pi**2*TC/UC/N**3+Pi**2/N**3/2+Pi**2/TC*UC*N/8+Pi**2*TC/U
     #C*N/8-Pi**2*N/SC**2*TC**2/2
      s7 = s8-Pi**2*N/SC**2*UC**2/2+N/8-N/SC**2*TC**2/4-N/SC**2*UC**2/4+
     #1/N/8-1/N*TC**2/SC**2/4-1/N*UC**2/SC**2/4+dlog(-TC/S)**2*(N*SC/TC/
     #4-N*UC/SC-N/4+1/N*TC/UC/2-1/N*UC/SC/4+1/TC*UC/N**3/4-1/N**3*SC/UC/
     #2)+dlog(-TC/S)*(N/SC**2*TC**2+N/SC**2*UC**2+3.D0/4.D0*N*TC/SC-5.D0
     #/4.D0/TC*UC*N-N/4-1/N*UC/SC/4-2/N*SC/UC-1/N*SC/TC/2-3.D0/4.D0/N**3
     #*SC/TC-1/N**3/4)+dlog(SC/S)*dlog(-TC/S)*(N/SC**2*TC**2+N/SC**2*UC*
     #*2-1/TC*UC*N/2+1/N*UC/SC/2-1/N*TC/UC+1/N**3*SC/UC-1/TC*UC/N**3/2)
      s5 = s6*s7
      s3 = s4*s5
      s5 = 4*N
      s7 = N**2-1
      s9 = dlog(-TC/S)*dlog(-UC/S)/N*(TC**2+UC**2)/TC/UC/2+(-Pi**2+dlog(
     #SC/S)**2)*(1/N**3*SC**2/TC/UC/4+1/N/8+1/N*UC/TC/4+1/N*TC/UC/4-1/N*
     #TC**2/SC**2/4-1/N*UC**2/SC**2/4-N/SC**2*TC**2/4-N/SC**2*UC**2/4)+d
     #log(SC/S)*(5.D0/8.D0*N-9.D0/8.D0/N-1/N**3-TC/UC*N/2-1/TC*UC*N/2-TC
     #/UC/N**3/2-1/TC*UC/N**3/2-N/SC**2*TC**2/4-N/SC**2*UC**2/4+1/N*TC**
     #2/SC**2/4+1/N*UC**2/SC**2/4)+Pi**2/N/8+3.D0/8.D0*Pi**2/TC*UC/N**3+
     #3.D0/8.D0*Pi**2*TC/UC/N**3+Pi**2/N**3/2+Pi**2/TC*UC*N/8+Pi**2*TC/U
     #C*N/8-Pi**2*N/SC**2*TC**2/2
      s8 = s9-Pi**2*N/SC**2*UC**2/2+N/8-N/SC**2*TC**2/4-N/SC**2*UC**2/4+
     #1/N/8-1/N*TC**2/SC**2/4-1/N*UC**2/SC**2/4+dlog(-UC/S)**2*(N*SC/UC/
     #4-N*TC/SC-N/4+1/N*UC/TC/2-1/N*TC/SC/4+TC/UC/N**3/4-1/N**3*SC/TC/2)
     #+dlog(-UC/S)*(N/SC**2*TC**2+N/SC**2*UC**2+3.D0/4.D0*N*UC/SC-5.D0/4
     #.D0*TC/UC*N-N/4-1/N*TC/SC/4-2/N*SC/TC-1/N*SC/UC/2-3.D0/4.D0/N**3*S
     #C/UC-1/N**3/4)+dlog(SC/S)*dlog(-UC/S)*(N/SC**2*TC**2+N/SC**2*UC**2
     #-TC/UC*N/2+1/N*TC/SC/2-1/N*UC/TC+1/N**3*SC/TC-TC/UC/N**3/2)
      s6 = s7*s8
      s4 = s5*s6
      s2 = s3+s4
      teres = s1+s2
      	VPARTVI(16)=(TERAJ+TERES)
C	j0 = 17 : qi + g --> qi + g
      s1 = (B_QQ()+B_GG())*dlog(ptm**2/S)*C(TC,SC,UC,N)+(B_QQ()+B_GG())*
     #dlog(pt**2/S)*C(TC,SC,UC,N)-(dlog(-TC/S)*((2*N**2*(N**2-1)+2*(N**2
     #-1)/N**2)*(SC**2+UC**2)/SC/UC-4*(N**2-1)**2*(SC**2+UC**2)/TC**2)+4
     #*N**2*(N**2-1)*(dlog(SC/S)*(UC/SC-2/TC**2*UC**2)+dlog(-UC/S)*(1/UC
     #*SC-2/TC**2*SC**2))-4*(N**2-1)*(UC/SC+1/UC*SC)*(dlog(SC/S)+dlog(-U
     #C/S)))*dlog(ptm**2/S)
      teraj = s1+32*(TC**2/4+SC*TC/4+SC**2/4)*(N-1)*(N+1)*(dlog(SC/S)*N*
     #*4*SC**2/2+dlog(SC/S)*N**4*SC*TC+dlog(SC/S)*N**4*TC**2/2-dlog(SC/S
     #)*N**2*TC**2/2+dlog(-TC/S)*N**4*TC**2/4+dlog(-TC/S)*N**4*SC**2/2+d
     #log(-TC/S)*SC*N**4*TC/2-dlog(-TC/S)*N**2*SC**2/2-dlog(-TC/S)*SC*N*
     #*2*TC/2+dlog(-TC/S)*TC**2/4+dlog(-UC/S)*N**4*SC**2/2-dlog(-UC/S)*N
     #**2*TC**2/2)/TC**2/N**2/SC/(SC/2+TC/2)-(A4_QQ(1.D0)+A4_GG(1.D0))*C
     #(TC,SC,UC,N)*(dlog(pt**2/S)**2/4-dlog(pt**2/S)*dlog(ptm**2/S)/2)+(
     #A4_QQ(1.D0)+A4_GG(1.D0))*C(TC,SC,UC,N)*dlog(ptm**2/S)**2/4
      s1 = -(-7*N+7/N+22.D0/3.D0*N*dlog(MU2/S)-8.D0/3.D0*GTR*dlog(MU2/S)
     #)*(N**2-1)/N*(1/UC/SC*N**2-1/UC/SC-2*N**2/TC**2)*(SC**2+UC**2)
      s4 = -4*N
      s6 = N**2-1
      s8 = dlog(SC/S)*dlog(-UC/S)/N*(SC**2+UC**2)/SC/UC/2+dlog(-TC/S)**2
     #*(1/N**3*TC**2/SC/UC/4+1/N/8+1/N*SC/UC/4+1/N*UC/SC/4-1/N*SC**2/TC*
     #*2/4-1/N*UC**2/TC**2/4-N/TC**2*SC**2/4-N/TC**2*UC**2/4)+dlog(-TC/S
     #)*(5.D0/8.D0*N-9.D0/8.D0/N-1/N**3-SC/UC*N/2-1/SC*UC*N/2-SC/UC/N**3
     #/2-1/SC*UC/N**3/2-N/TC**2*SC**2/4-N/TC**2*UC**2/4+1/N*SC**2/TC**2/
     #4+1/N*UC**2/TC**2/4)+Pi**2/N/8+3.D0/8.D0*Pi**2*SC/UC/N**3+3.D0/8.D
     #0*Pi**2/SC*UC/N**3+Pi**2/N**3/2+Pi**2*SC/UC*N/8+Pi**2/SC*UC*N/8-Pi
     #**2*N/TC**2*SC**2/2
      s7 = s8-Pi**2*N/TC**2*UC**2/2+N/8-N/TC**2*SC**2/4-N/TC**2*UC**2/4+
     #1/N/8-1/N*SC**2/TC**2/4-1/N*UC**2/TC**2/4+(-Pi**2+dlog(SC/S)**2)*(
     #N*TC/SC/4-N*UC/TC-N/4+1/N*SC/UC/2-1/N*UC/TC/4+1/SC*UC/N**3/4-1/N**
     #3*TC/UC/2)+dlog(SC/S)*(N/TC**2*SC**2+N/TC**2*UC**2+3.D0/4.D0*N*SC/
     #TC-5.D0/4.D0/SC*UC*N-N/4-1/N*UC/TC/4-2/N*TC/UC-1/N*TC/SC/2-3.D0/4.
     #D0/N**3*TC/SC-1/N**3/4)+dlog(-TC/S)*dlog(SC/S)*(N/TC**2*SC**2+N/TC
     #**2*UC**2-1/SC*UC*N/2+1/N*UC/TC/2-1/N*SC/UC+1/N**3*TC/UC-1/SC*UC/N
     #**3/2)
      s5 = s6*s7
      s3 = s4*s5
      s5 = -4*N
      s7 = N**2-1
      s9 = dlog(SC/S)*dlog(-UC/S)/N*(SC**2+UC**2)/SC/UC/2+dlog(-TC/S)**2
     #*(1/N**3*TC**2/SC/UC/4+1/N/8+1/N*SC/UC/4+1/N*UC/SC/4-1/N*SC**2/TC*
     #*2/4-1/N*UC**2/TC**2/4-N/TC**2*SC**2/4-N/TC**2*UC**2/4)+dlog(-TC/S
     #)*(5.D0/8.D0*N-9.D0/8.D0/N-1/N**3-SC/UC*N/2-1/SC*UC*N/2-SC/UC/N**3
     #/2-1/SC*UC/N**3/2-N/TC**2*SC**2/4-N/TC**2*UC**2/4+1/N*SC**2/TC**2/
     #4+1/N*UC**2/TC**2/4)+Pi**2/N/8+3.D0/8.D0*Pi**2*SC/UC/N**3+3.D0/8.D
     #0*Pi**2/SC*UC/N**3+Pi**2/N**3/2+Pi**2*SC/UC*N/8+Pi**2/SC*UC*N/8-Pi
     #**2*N/TC**2*SC**2/2
      s8 = s9-Pi**2*N/TC**2*UC**2/2+N/8-N/TC**2*SC**2/4-N/TC**2*UC**2/4+
     #1/N/8-1/N*SC**2/TC**2/4-1/N*UC**2/TC**2/4+dlog(-UC/S)**2*(N*TC/UC/
     #4-N*SC/TC-N/4+1/N*UC/SC/2-1/N*SC/TC/4+SC/UC/N**3/4-1/N**3*TC/SC/2)
     #+dlog(-UC/S)*(N/TC**2*SC**2+N/TC**2*UC**2+3.D0/4.D0*N*UC/TC-5.D0/4
     #.D0*SC/UC*N-N/4-1/N*SC/TC/4-2/N*TC/SC-1/N*TC/UC/2-3.D0/4.D0/N**3*T
     #C/UC-1/N**3/4)+dlog(-TC/S)*dlog(-UC/S)*(N/TC**2*SC**2+N/TC**2*UC**
     #2-SC/UC*N/2+1/N*SC/TC/2-1/N*UC/SC+1/N**3*TC/SC-SC/UC/N**3/2)
      s6 = s7*s8
      s4 = s5*s6
      s2 = s3+s4
      teres = s1+s2
      	VPARTVI(17)=(TERAJ+TERES)
C	j0 = 18 : qi + g --> g + g
      teraj = 0
      teres = 0
      	VPARTVI(18)=(TERAJ+TERES)
C	j0 = 19 : g + g --> qi +qbi
      s1 = -2*B_GG()*dlog(ptm**2/S)*C(SC,TC,UC,N)-2*B_QQ()*dlog(pt**2/S)
     #*C(SC,TC,UC,N)+(dlog(SC/S)*((2*N**2*(N**2-1)+2*(N**2-1)/N**2)*(TC*
     #*2+UC**2)/TC/UC-4*(N**2-1)**2*(TC**2+UC**2)/SC**2)+4*N**2*(N**2-1)
     #*(dlog(-TC/S)*(UC/TC-2*UC**2/SC**2)+dlog(-UC/S)*(TC/UC-2*TC**2/SC*
     #*2))-4*(N**2-1)*(UC/TC+TC/UC)*(dlog(-TC/S)+dlog(-UC/S)))*dlog(ptm*
     #*2/S)
      teraj = s1-32*(SC*TC/4+SC**2/4+TC**2/4)*(N-1)*(N+1)*(dlog(SC/S)*N*
     #*4*SC**2/4+dlog(SC/S)*SC*N**4*TC/2+dlog(SC/S)*N**4*TC**2/2-dlog(SC
     #/S)*SC*N**2*TC/2-dlog(SC/S)*TC**2*N**2/2+dlog(SC/S)*SC**2/4+dlog(-
     #TC/S)*N**4*TC**2/2+dlog(-TC/S)*N**4*SC*TC+dlog(-TC/S)*N**4*SC**2/2
     #-dlog(-TC/S)*N**2*SC**2/2+dlog(-UC/S)*N**4*TC**2/2-dlog(-UC/S)*N**
     #2*SC**2/2)/SC**2/TC/N**2/(SC/2+TC/2)+2*A4_QQ(1.D0)*C(SC,TC,UC,N)*(
     #dlog(pt**2/S)**2/4-dlog(pt**2/S)*dlog(ptm**2/S)/2)-A4_GG(1.D0)*C(S
     #C,TC,UC,N)*dlog(ptm**2/S)**2/2
      s1 = (-7*N+7/N+22.D0/3.D0*N*dlog(MU2/S)-8.D0/3.D0*GTR*dlog(MU2/S))
     #*(N**2-1)/N*(1/UC/TC*N**2-1/UC/TC-2*N**2/SC**2)*(TC**2+UC**2)
      s4 = 4*N
      s6 = N**2-1
      s8 = dlog(-TC/S)*dlog(-UC/S)/N*(TC**2+UC**2)/TC/UC/2+(-Pi**2+dlog(
     #SC/S)**2)*(1/N**3*SC**2/TC/UC/4+1/N/8+1/N*TC/UC/4+1/N*UC/TC/4-1/N/
     #SC**2*TC**2/4-1/N/SC**2*UC**2/4-N/SC**2*TC**2/4-N/SC**2*UC**2/4)+d
     #log(SC/S)*(5.D0/8.D0*N-9.D0/8.D0/N-1/N**3-TC/UC*N/2-1/TC*UC*N/2-TC
     #/UC/N**3/2-1/TC*UC/N**3/2-N/SC**2*TC**2/4-N/SC**2*UC**2/4+1/N/SC**
     #2*TC**2/4+1/N/SC**2*UC**2/4)+Pi**2/N/8+3.D0/8.D0*Pi**2*TC/UC/N**3+
     #3.D0/8.D0*Pi**2/TC*UC/N**3+Pi**2/N**3/2+Pi**2*TC/UC*N/8+Pi**2/TC*U
     #C*N/8-Pi**2*N/SC**2*TC**2/2
      s7 = s8-Pi**2*N/SC**2*UC**2/2+N/8-N/SC**2*TC**2/4-N/SC**2*UC**2/4+
     #1/N/8-1/N/SC**2*TC**2/4-1/N/SC**2*UC**2/4+dlog(-TC/S)**2*(N*SC/TC/
     #4-N*UC/SC-N/4+1/N*TC/UC/2-1/N*UC/SC/4+1/TC*UC/N**3/4-1/N**3*SC/UC/
     #2)+dlog(-TC/S)*(N/SC**2*TC**2+N/SC**2*UC**2+3.D0/4.D0*N*TC/SC-5.D0
     #/4.D0/TC*UC*N-N/4-1/N*UC/SC/4-2/N*SC/UC-1/N*SC/TC/2-3.D0/4.D0/N**3
     #*SC/TC-1/N**3/4)+dlog(SC/S)*dlog(-TC/S)*(N/SC**2*TC**2+N/SC**2*UC*
     #*2-1/TC*UC*N/2+1/N*UC/SC/2-1/N*TC/UC+1/N**3*SC/UC-1/TC*UC/N**3/2)
      s5 = s6*s7
      s3 = s4*s5
      s5 = 4*N
      s7 = N**2-1
      s9 = dlog(-TC/S)*dlog(-UC/S)/N*(TC**2+UC**2)/TC/UC/2+(-Pi**2+dlog(
     #SC/S)**2)*(1/N**3*SC**2/TC/UC/4+1/N/8+1/N*TC/UC/4+1/N*UC/TC/4-1/N/
     #SC**2*TC**2/4-1/N/SC**2*UC**2/4-N/SC**2*TC**2/4-N/SC**2*UC**2/4)+d
     #log(SC/S)*(5.D0/8.D0*N-9.D0/8.D0/N-1/N**3-TC/UC*N/2-1/TC*UC*N/2-TC
     #/UC/N**3/2-1/TC*UC/N**3/2-N/SC**2*TC**2/4-N/SC**2*UC**2/4+1/N/SC**
     #2*TC**2/4+1/N/SC**2*UC**2/4)+Pi**2/N/8+3.D0/8.D0*Pi**2*TC/UC/N**3+
     #3.D0/8.D0*Pi**2/TC*UC/N**3+Pi**2/N**3/2+Pi**2*TC/UC*N/8+Pi**2/TC*U
     #C*N/8-Pi**2*N/SC**2*TC**2/2
      s8 = s9-Pi**2*N/SC**2*UC**2/2+N/8-N/SC**2*TC**2/4-N/SC**2*UC**2/4+
     #1/N/8-1/N/SC**2*TC**2/4-1/N/SC**2*UC**2/4+dlog(-UC/S)**2*(N*SC/UC/
     #4-N*TC/SC-N/4+1/N*UC/TC/2-1/N*TC/SC/4+TC/UC/N**3/4-1/N**3*SC/TC/2)
     #+dlog(-UC/S)*(N/SC**2*TC**2+N/SC**2*UC**2+3.D0/4.D0*N*UC/SC-5.D0/4
     #.D0*TC/UC*N-N/4-1/N*TC/SC/4-2/N*SC/TC-1/N*SC/UC/2-3.D0/4.D0/N**3*S
     #C/UC-1/N**3/4)+dlog(SC/S)*dlog(-UC/S)*(N/SC**2*TC**2+N/SC**2*UC**2
     #-TC/UC*N/2+1/N*TC/SC/2-1/N*UC/TC+1/N**3*SC/TC-TC/UC/N**3/2)
      s6 = s7*s8
      s4 = s5*s6
      s2 = s3+s4
      teres = s1+s2
      	VPARTVI(19)=(TERAJ+TERES)
C	j0 = 20 : g + g --> qi + g
      teraj = 0
      teres = 0
      	VPARTVI(20)=(TERAJ+TERES)
C	j0 = 21 : g + g --> g + g
      s1 = -2*B_GG()*dlog(ptm**2/S)*Dg(SC,TC,UC,N)-2*B_GG()*dlog(pt**2/S
     #)*Dg(SC,TC,UC,N)+(16*(N**2-1)*N**3*dlog(SC/S)*(3-2*TC*UC/SC**2+(TC
     #**4+UC**4)/UC**2/TC**2)+16*(N**2-1)*N**3*dlog(-TC/S)*(3-2*UC*SC/TC
     #**2+(UC**4+SC**4)/SC**2/UC**2)+16*(N**2-1)*N**3*dlog(-UC/S)*(3-2*S
     #C*TC/UC**2+(SC**4+TC**4)/SC**2/TC**2))*dlog(ptm**2/S)
      teraj = s1+512*N**3*(SC**2/4+TC*SC/4+TC**2/4)**2*(N-1)*(N+1)*(dlog
     #(SC/S)*SC**2/4+dlog(SC/S)*TC*SC/2+dlog(SC/S)*TC**2/2+dlog(-TC/S)*S
     #C**2/2+dlog(-TC/S)*TC*SC/2+dlog(-TC/S)*TC**2/4+dlog(-UC/S)*SC**2/4
     #+dlog(-UC/S)*TC**2/4)/TC**2/SC**2/(SC/2+TC/2)**2+2*A4_GG(1.D0)*Dg(
     #SC,TC,UC,N)*(dlog(pt**2/S)**2/4-dlog(pt**2/S)*dlog(ptm**2/S)/2)-A4
     #_GG(1.D0)*Dg(SC,TC,UC,N)*dlog(ptm**2/S)**2/2
      s1 = (-1072.D0/9.D0*N+320.D0/9.D0*GTR+16*N*Pi**2+176.D0/3.D0*N*dlo
     #g(MU2/S)-64.D0/3.D0*GTR*dlog(MU2/S))*(N**2-1)*N**2*(3-TC*UC/SC**2-
     #UC*SC/TC**2-SC*TC/UC**2)
      s3 = 4*N**2-4
      s5 = N**2
      s7 = N*((2*TC**2+2*UC**2)/TC/UC*(-Pi**2+dlog(SC/S)**2)+(4*SC*TC/UC
     #**2+4*UC*SC/TC**2-6)*dlog(-TC/S)*dlog(-UC/S)+(4.D0/3.D0*TC*UC/SC**
     #2-14.D0/3.D0*TC/UC-14.D0/3.D0/TC*UC-14-8*TC**2/UC**2-8*UC**2/TC**2
     #)*dlog(SC/S)-1-Pi**2)+GTR*((10.D0/3.D0*TC/UC+10.D0/3.D0/TC*UC+16.D
     #0/3.D0*TC*UC/SC**2-2)*dlog(SC/S)+(-SC**2-UC*TC)/TC/UC*(-Pi**2+dlog
     #(SC/S)**2)+(-2*TC**2-2*UC**2)/TC/UC*dlog(-TC/S)*dlog(-UC/S)+2-Pi**
     #2)+N*(-(-2*UC**2-2*SC**2)/UC/SC*dlog(-TC/S)**2+(4*TC*UC/SC**2+4*SC
     #*TC/UC**2-6)*dlog(-UC/S)*dlog(SC/S)+(4.D0/3.D0*UC*SC/TC**2-14.D0/3
     #.D0*UC/SC-14.D0/3.D0/UC*SC-14-8*UC**2/SC**2-8*SC**2/UC**2)*dlog(-T
     #C/S)-1-Pi**2)
      s6 = s7+GTR*((10.D0/3.D0*UC/SC+10.D0/3.D0/UC*SC+16.D0/3.D0*UC*SC/T
     #C**2-2)*dlog(-TC/S)-(TC**2+UC*SC)/UC/SC*dlog(-TC/S)**2-(2*UC**2+2*
     #SC**2)/UC/SC*dlog(-UC/S)*dlog(SC/S)+2-Pi**2)+N*(-(-2*SC**2-2*TC**2
     #)/SC/TC*dlog(-UC/S)**2+(4*UC*SC/TC**2+4*TC*UC/SC**2-6)*dlog(SC/S)*
     #dlog(-TC/S)+(4.D0/3.D0*SC*TC/UC**2-14.D0/3.D0*SC/TC-14.D0/3.D0/SC*
     #TC-14-8*SC**2/TC**2-8*TC**2/SC**2)*dlog(-UC/S)-1-Pi**2)+GTR*((10.D
     #0/3.D0*SC/TC+10.D0/3.D0/SC*TC+16.D0/3.D0*SC*TC/UC**2-2)*dlog(-UC/S
     #)-(UC**2+TC*SC)/SC/TC*dlog(-UC/S)**2-(2*SC**2+2*TC**2)/SC/TC*dlog(
     #SC/S)*dlog(-TC/S)+2-Pi**2)
      s4 = s5*s6
      s2 = s3*s4
      teres = s1+s2
      	VPARTVI(21)=(TERAJ+TERES)
	RETURN
	END
	SUBROUTINE VEC_VI2T(S,SC,TC,UC,FI,VPARTVI2)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/COUL/N,CF,GTR
	PARAMETER (J0MAX=21)
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
	CALL VEC_H34T(S12,S13,S14,S15,S23,S24,S25,S34,S35,S45
     #	,HH34)
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
c le 28/06/00, correction d'un bug du a l'optimisation du compilateur
C*****************************************************
	SUBROUTINE STRFRAT(X1,IH1,X2,IH2,X3,IH3,X4,IH4,SFD)
	IMPLICIT REAL*8 (A-H,L-Z)
	COMMON/GDF/JNF
	COMMON/SCALE/M,MF,MU
	DIMENSION F1(-6:6),F2(-6:6),D3(-6:6),D4(-6:6)
	PARAMETER (K0MAX=84,J0MAX=21)
	DIMENSION SF(K0MAX),SFD(K0MAX)
	CALL FSTRU(X1,M*M,IH1,F1)
	CALL FSTRU(X2,M*M,IH2,F2)
	CALL DFRAG(X3,MF*MF,IH3,D3)
	CALL DFRAG(X4,MF*MF,IH4,D4)
C initialisation du tableau SF(K0MAX)
	ZERO = 0.D0
	CALL VEC_DINIT(SF,K0MAX,ZERO)
	CALL VEC_DINIT(SFD,K0MAX,ZERO)
C on commence les processus
	DO I=1,JNF-1
	  DO K=I+1,JNF
C	j0 = 1 : qi + qk --> qi + qk
	    SF(1) = F1(I)*F2(K)*D3(I)*D4(K) + 	    
     #              F1(-I)*F2(-K)*D3(-I)*D4(-K) + SF(1)
	    SF(1+21) = F1(K)*F2(I)*D3(I)*D4(K) + 	    
     #                 F1(-K)*F2(-I)*D3(-I)*D4(-K) + SF(1+21)
	    SF(1+42) = F1(I)*F2(K)*D3(K)*D4(I) + 	    
     #                 F1(-I)*F2(-K)*D3(-K)*D4(-I) + SF(1+42)
	    SF(1+63) = F1(K)*F2(I)*D3(K)*D4(I) + 	    
     #                 F1(-K)*F2(-I)*D3(-K)*D4(-I) + SF(1+63)
C	j0 = 3 : qi + qbk --> qi + qbk
	    SF(3) = F1(I)*F2(-K)*D3(I)*D4(-K) + 	    
     #              F1(-I)*F2(K)*D3(-I)*D4(K) + SF(3)
	    SF(3+21) = F1(K)*F2(-I)*D3(I)*D4(-K) + 	    
     #                 F1(-K)*F2(I)*D3(-I)*D4(K) + SF(3+21)
	    SF(3+42) = F1(I)*F2(-K)*D3(K)*D4(-I) + 	    
     #                 F1(-I)*F2(K)*D3(-K)*D4(I) + SF(3+42)
	    SF(3+63) = F1(K)*F2(-I)*D3(K)*D4(-I) + 	    
     #                 F1(-K)*F2(I)*D3(-K)*D4(I) + SF(3+63)
	  ENDDO
	ENDDO 	    
	DO I=1,JNF
	  DO K=1,JNF
	    IF (K.NE.I) THEN
C	j0 = 2 : qi + qk --> qi + g
	      SF(2) = F1(I)*F2(K)*D3(I)*D4(0) + 	    
     #                F1(-I)*F2(-K)*D3(-I)*D4(0) + SF(2)
	      SF(2+21) = F1(K)*F2(I)*D3(I)*D4(0) + 	    
     #                   F1(-K)*F2(-I)*D3(-I)*D4(0) + SF(2+21)
	      SF(2+42) = F1(I)*F2(K)*D3(0)*D4(I) + 	    
     #                   F1(-I)*F2(-K)*D3(0)*D4(-I) + SF(2+42)
	      SF(2+63) = F1(K)*F2(I)*D3(0)*D4(I) + 	    
     #                   F1(-K)*F2(-I)*D3(0)*D4(-I) + SF(2+63)
C	j0 = 4 : qi + qbk --> qi + g
	      SF(4) = F1(I)*F2(-K)*D3(I)*D4(0) + 	    
     #                F1(-I)*F2(K)*D3(-I)*D4(0) + SF(4)
	      SF(4+21) = F1(K)*F2(-I)*D3(I)*D4(0) + 	    
     #                   F1(-K)*F2(I)*D3(-I)*D4(0) + SF(4+21)
	      SF(4+42) = F1(I)*F2(-K)*D3(0)*D4(I) + 	    
     #                   F1(-I)*F2(K)*D3(0)*D4(-I) + SF(4+42)
	      SF(4+63) = F1(K)*F2(-I)*D3(0)*D4(I) + 	    
     #                   F1(-K)*F2(I)*D3(0)*D4(-I) + SF(4+63)
C	j0 = 5 : qi + qbi --> qk + qbk
	      SF(5) = F1(I)*F2(-I)*D3(K)*D4(-K) + SF(5)
	      SF(5+21) = F1(-I)*F2(I)*D3(K)*D4(-K) + SF(5+21)
	      SF(5+42) = F1(I)*F2(-I)*D3(-K)*D4(K) + SF(5+42)
	      SF(5+63) = F1(-I)*F2(I)*D3(-K)*D4(K) + SF(5+63)
C	j0 = 6 : qi + qbi --> qk + g
	      SF(6) = F1(I)*F2(-I)*D3(K)*D4(0) + 	    
     #                F1(-I)*F2(I)*D3(-K)*D4(0) + SF(6)
	      SF(6+21) = F1(-I)*F2(I)*D3(K)*D4(0) + 	    
     #                   F1(I)*F2(-I)*D3(-K)*D4(0) + SF(6+21)
	      SF(6+42) = F1(I)*F2(-I)*D3(0)*D4(K) + 	    
     #                   F1(-I)*F2(I)*D3(0)*D4(-K) + SF(6+42)
	      SF(6+63) = F1(-I)*F2(I)*D3(0)*D4(K) + 	    
     #                   F1(I)*F2(-I)*D3(0)*D4(-K) + SF(6+63)
C	j0 = 7 : qi + g --> qi + qk
	      SF(7) = F1(I)*F2(0)*D3(I)*D4(K) + 	    
     #                F1(-I)*F2(0)*D3(-I)*D4(K) + SF(7)
	      SF(7+21) = F1(0)*F2(I)*D3(I)*D4(K) + 	    
     #                   F1(0)*F2(-I)*D3(-I)*D4(K) + SF(7+21)
	      SF(7+42) = F1(I)*F2(0)*D3(K)*D4(I) + 	    
     #                   F1(-I)*F2(0)*D3(K)*D4(-I) + SF(7+42)
	      SF(7+63) = F1(0)*F2(I)*D3(K)*D4(I) + 	    
     #                   F1(0)*F2(-I)*D3(K)*D4(-I) + SF(7+63)
C	j0 = 8 : qi + g --> qi + qbk
	      SF(8) = F1(I)*F2(0)*D3(I)*D4(-K) + 	    
     #                F1(-I)*F2(0)*D3(-I)*D4(-K) + SF(8)
	      SF(8+21) = F1(0)*F2(I)*D3(I)*D4(-K) + 	    
     #                   F1(0)*F2(-I)*D3(-I)*D4(-K) + SF(8+21)
	      SF(8+42) = F1(I)*F2(0)*D3(-K)*D4(I) + 	    
     #                   F1(-I)*F2(0)*D3(-K)*D4(-I) + SF(8+42)
	      SF(8+63) = F1(0)*F2(I)*D3(-K)*D4(I) + 	    
     #                   F1(0)*F2(-I)*D3(-K)*D4(-I) + SF(8+63)
C	j0 = 9 : qi + g --> qk + qbk
	      SF(9) = F1(I)*F2(0)*D3(K)*D4(-K) + 	    
     #                F1(-I)*F2(0)*D3(K)*D4(-K) + SF(9)
	      SF(9+21) = F1(0)*F2(I)*D3(K)*D4(-K) + 	    
     #                   F1(0)*F2(-I)*D3(K)*D4(-K) + SF(9+21)
	      SF(9+42) = F1(I)*F2(0)*D3(-K)*D4(K) + 	    
     #                   F1(-I)*F2(0)*D3(-K)*D4(K) + SF(9+42)
	      SF(9+63) = F1(0)*F2(I)*D3(-K)*D4(K) + 	    
     #                   F1(0)*F2(-I)*D3(-K)*D4(K) + SF(9+63)
	    ENDIF
	  ENDDO
	ENDDO 	    
	DO I=1,JNF
C	j0 = 10 : qi + qi --> qi + qi
	  SF(10) = (F1(I)*F2(I)*D3(I)*D4(I) + 	    
     #              F1(-I)*F2(-I)*D3(-I)*D4(-I))/4.D0 + SF(10)
	  SF(10+42) = (F1(I)*F2(I)*D3(I)*D4(I) + 	    
     #              F1(-I)*F2(-I)*D3(-I)*D4(-I))/4.D0 + SF(10+42)
C	j0 = 11 : qi + qi --> qi + g
	  SF(11) = (F1(I)*F2(I)*D3(I)*D4(0) + 	    
     #              F1(-I)*F2(-I)*D3(-I)*D4(0))/2.D0 + SF(11)
	  SF(11+42) = (F1(I)*F2(I)*D3(0)*D4(I) + 	    
     #                 F1(-I)*F2(-I)*D3(0)*D4(-I))/2.D0 + SF(11+42)
C	j0 = 12 : qi + qbi --> qi + qbi
	  SF(12) = F1(I)*F2(-I)*D3(I)*D4(-I) + SF(12)
	  SF(12+21) = F1(-I)*F2(I)*D3(-I)*D4(I) + SF(12+21)
	  SF(12+42) = F1(I)*F2(-I)*D3(-I)*D4(I) + SF(12+42)
	  SF(12+63) = F1(-I)*F2(I)*D3(I)*D4(-I) + SF(12+63)
C	j0 = 13 : qi + qbi --> qi + g
	  SF(13) = F1(I)*F2(-I)*D3(I)*D4(0) + 	    
     #             F1(-I)*F2(I)*D3(-I)*D4(0) + SF(13)
	  SF(13+21) = F1(-I)*F2(I)*D3(I)*D4(0) + 	    
     #                F1(I)*F2(-I)*D3(-I)*D4(0) + SF(13+21)
	  SF(13+42) = F1(I)*F2(-I)*D3(0)*D4(I) + 	    
     #                F1(-I)*F2(I)*D3(0)*D4(-I) + SF(13+42)
	  SF(13+63) = F1(-I)*F2(I)*D3(0)*D4(I) + 	    
     #                F1(I)*F2(-I)*D3(0)*D4(-I) + SF(13+63)
C	j0 = 14 : qi + g --> qi + qi
	  SF(14) = (F1(I)*F2(0)*D3(I)*D4(I) + 	    
     #              F1(-I)*F2(0)*D3(-I)*D4(-I))/2.D0 + SF(14)
	  SF(14+21) = (F1(0)*F2(I)*D3(I)*D4(I) + 	    
     #                 F1(0)*F2(-I)*D3(-I)*D4(-I))/2.D0 + SF(14+21)
C	j0 = 15 : qi + g --> qi + qbi
	  SF(15) = F1(I)*F2(0)*D3(I)*D4(-I) + 	    
     #             F1(-I)*F2(0)*D3(I)*D4(-I) + SF(15)
	  SF(15+21) = F1(0)*F2(I)*D3(I)*D4(-I) + 	    
     #                F1(0)*F2(-I)*D3(I)*D4(-I) + SF(15+21)
	  SF(15+42) = F1(I)*F2(0)*D3(-I)*D4(I) + 	    
     #                F1(-I)*F2(0)*D3(-I)*D4(I) + SF(15+42)
	  SF(15+63) = F1(0)*F2(I)*D3(-I)*D4(I) + 	    
     #                F1(0)*F2(-I)*D3(-I)*D4(I) + SF(15+63)
C	j0 = 16 : qi + qbi --> g + g
	  SF(16) = (F1(I)*F2(-I)*D3(0)*D4(0))/2.D0 + SF(16)
	  SF(16+21) = (F1(-I)*F2(I)*D3(0)*D4(0))/2.D0 + SF(16+21)
C	j0 = 17 : qi + g --> qi + g
	  SF(17) = F1(I)*F2(0)*D3(I)*D4(0) + 	    
     #             F1(-I)*F2(0)*D3(-I)*D4(0) + SF(17)
	  SF(17+21) = F1(0)*F2(I)*D3(I)*D4(0) + 	    
     #                F1(0)*F2(-I)*D3(-I)*D4(0) + SF(17+21)
	  SF(17+42) = F1(I)*F2(0)*D3(0)*D4(I) + 	    
     #                F1(-I)*F2(0)*D3(0)*D4(-I) + SF(17+42)
	  SF(17+63) = F1(0)*F2(I)*D3(0)*D4(I) + 	    
     #                F1(0)*F2(-I)*D3(0)*D4(-I) + SF(17+63)
C	j0 = 18 : qi + g --> g + g
	  SF(18) = (F1(I)*F2(0)*D3(0)*D4(0) + 	    
     #              F1(-I)*F2(0)*D3(0)*D4(0))/2.D0 + SF(18)
	  SF(18+21) = (F1(0)*F2(I)*D3(0)*D4(0) + 	    
     #                 F1(0)*F2(-I)*D3(0)*D4(0))/2.D0 + SF(18+21)
C	j0 = 19 : g + g --> qi +qbi
	  SF(19) = (F1(0)*F2(0)*D3(I)*D4(-I))/2.D0 + SF(19)
	  SF(19+42) = (F1(0)*F2(0)*D3(-I)*D4(I))/2.D0 + SF(19+42)
C	j0 = 20 : g + g --> qi + g
	  SF(20) = (F1(0)*F2(0)*D3(I)*D4(0) + 	    
     #              F1(0)*F2(0)*D3(-I)*D4(0))/2.D0 + SF(20)
	  SF(20+42) = (F1(0)*F2(0)*D3(0)*D4(I) + 	    
     #                 F1(0)*F2(0)*D3(0)*D4(-I))/2.D0 + SF(20+42)
	ENDDO 	    
C	j0 = 21 : g + g --> g + g
	SF(21) = F1(0)*F2(0)*D3(0)*D4(0)/4.D0
	SF(21+42) = F1(0)*F2(0)*D3(0)*D4(0)/4.D0
C
     	SF(10+21) = SF(10)
     	SF(10+63) = SF(10+42)
C
     	SF(11+21) = SF(11)
     	SF(11+63) = SF(11+42)
C
     	SF(14+42) = SF(14)
     	SF(14+63) = SF(14+21)
C
     	SF(16+42) = SF(16)
     	SF(16+63) = SF(16+21)
C
     	SF(18+42) = SF(18)
     	SF(18+63) = SF(18+21)
C
     	SF(19+21) = SF(19)
     	SF(19+63) = SF(19+42)
C
     	SF(20+21) = SF(20)
     	SF(20+63) = SF(20+42)
C
     	SF(21+21) = SF(21)
     	SF(21+63) = SF(21+42)
C on divise tout par x1 x2
	XX = 1.D0/(X1*X2)
	CALL VEC_DMULT_CONSTANT(SF,K0MAX,XX,SFD)	    
	RETURN
	END
