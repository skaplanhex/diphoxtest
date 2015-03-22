c
c modification le 21/01/97 pour generer le fi_3 entre 0 et 2 pi
c modification le 22/01/97 pour la separation entre les deux photons
c modification le 14/02/97 pour mettre la coupure sur la masse
c invariante dans les quasi deux donne deux et les deux donne deux
c modification le 7/03/97 pour fixer le nombre total d'evenements au
c depart
c modification le 7/04/97 pour fixer un bug dans f2dim2 et f2dim3
c (modifications aussi dans onefrag.f: subroutine ssv3o et ssv4o)
c modification le 9/06/97 pour l'integration dans f4dim3 et f4dim4
c modification le 3/05/99 pour ajout des termes en r^2*ln(p_{tm})
c modification le 14/05/99 pour introduire le maximum de la masse 
c invariante photon-photon
c dernieres modifs + verification du code par thomas (19/09/99)
c bug dans la generation des 2 -> 3 (4//5) (15/05/00)
c modifications le 31/03/01 pour ameliorer l'efficacite en cas d'isolement
c modification le 22/07/2004 pour la correction d'un bug dans f4dimo3 (fimin
C etait redefinie)
	subroutine onef_sub
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical lo,nlo,gener,integ,intega
        logical lhalf
	logical facteur_symetrie
	character*128 path_bsfile
        integer*4 ntrack,iprov
**         integer*4 ntrack,iprov,maxtrk
**         parameter (maxtrk = 3)
        parameter (irandom=500)
	common/npt/inpt
	common/bparm1/xl(50),xu(50),idim,iwild,ig(50),icall
	common/bparm2/acc1,acc2,itmx1,itmx2
	common/born/iborn
	common/baspring/ispring
** 	common/fixle/iprov,weight
	common/approx/lo,nlo
	common/calcul/integ,intega,gener
 	common/nbevent/inbevent
	common/baseparam/jtmx1,jtmx2,jcall,jcall1,jcall2,ixtry
	common/accuracy/accu
** 	common/sortier/ntrack,we(maxtrk),wpx(maxtrk),
**      #	wpy(maxtrk),wpz(maxtrk)
** 	common/sfragvr/wx3,wx4
	common/randd/rnumb
	common/cheminbs/path_bsfile
	common/long/ilen
	common/symetrie/facteur_symetrie
	common/efficacite/ipositif,inegatif
	dimension wrvec(irandom),wrvec1(irandom),wrvec2(irandom)
	external f4dimo3,f4dimo4,f2dimo1,f2dimo2,f2dimo3,f2dimo3p,
     #	f1dimo,f4dimo4s
c
	pi=datan(1.d0)*4.d0
	ipositif = 0
	inegatif = 0
	half = 0.5d0
c******************************************************************
c       partie 2 --> 3
c******************************************************************
c
cccc
	if (nlo) then
cccc
	iborn = 0
c###
** 	goto 112
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimo3,resf34,sdf34,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimo3,resf34,sdf34,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'threepart1.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
c###
c
c###
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimo4,resf44,sdf44,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimo4,resf44,sdf44,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'threepart2.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
c###
c###
c cette contribution est en general peu importante, 
c pour ne pas perdre du temps, on lui reduit son nombre d'iteration
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
** 	itmx1 = jtmx1
** 	itmx2 = jtmx2
	itmx1 = 5
	itmx2 = 5
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimo4s,resf44s,sdf44s,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
** 	itmx1 = jtmx1
** 	itmx2 = jtmx2
	itmx1 = 5
	itmx2 = 5
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimo4s,resf44s,sdf44s,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'threepart2s.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
c###
c******************************************************************
c       partie 2 --> 2 colineaire
c******************************************************************
c
c###
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo1,resf21,sdf21,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo1,resf21,sdf21,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'qtwopart1.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
c###
c
c###
**  112	continue
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo2,resf22,sdf22,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo2,resf22,sdf22,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'qtwopart2.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
** 	goto 111
c###
c
c###
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo3,resf23,sdf23,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo3,resf23,sdf23,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'qtwopart3.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
c###
c
c###
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo3p,resf23p,sdf23p,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimo3p,resf23p,sdf23p,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'qtwopart3p.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
c###
c
c******************************************************************
c       partie 2 --> 2 virtuelle et born
c******************************************************************
c###
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall2
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	call bases (f1dimo,resf1,sdf1,ctime,it1,it2)
	endif
c###
c
c###
	if (intega) then
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall2
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	call bases (f1dimo,resf1,sdf1,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'twopart1.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c	
	endif
c###
**  111	continue
c
cccc
	endif
cccc
cccc
	if (lo) then
cccc
c
	if (intega.or.integ) then
	iborn = 1
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall2
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	call bases (f1dimo,resfb,sdfb,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'twopart2.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
c
	endif
cccc
	endif
cccc
c******************************************************************
c       calcul de la section efficace totale ordre superieur
c******************************************************************
	if (integ) then
	  resf2 = resf21+resf22+resf23+resf23p
	  sdf2 = sdf21+sdf22+sdf23+sdf23p
c
	  resho = resf34+resf44+resf44s+resf2+resf1
	  sdho = sdf34+sdf44+sdf44s+sdf2+sdf1
	  resb = resfb
	  sdfb = sdfb
	  write (8,*) 'cccccccccccccccccccccc'
	  write (8,*) 'born',resfb,' sdfb',sdfb
	  write (8,*) 'h.o.',resho,' sdfho',sdho
	  write (8,*) 'cccccccccccccccccccccc'
	endif
c******************************************************************
c       calcul du nombre d'evenements a generer pour chacune des
c       parties
c******************************************************************
	if (intega) then
	  if (nlo) then
	    resf2 = resf21+resf22+resf23+resf23p
	    totho = resf34+resf44+resf44s+resf2+resf1
	  else
	    totho = 0.d0
	  endif
	  if (lo) then
	    totb = resfb
	  else
	    totb = 0.d0
	  endif
	  tot = totho + totb
	  if (tot.eq.0.d0) then
	    write (8,*) 'tot=0. job finished'
	    return
	  endif
	  open(28,file=path_bsfile(1:ilen)//'integral.res',
     #	  status='unknown')
	  write(28,*) resf34, resf44, resf44s
	  write(28,*) resf21, resf22, resf23, resf23p
	  write(28,*) resf1,resfb
	  write(28,*) tot
	  close(28)	  
	endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc              generation des evenements                     ccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c******************************************************************
cccc
	if (gener) then
cccc
	open(28,file=path_bsfile(1:ilen)//'integral.res',
     #	status='unknown')
	read(28,*) resf34, resf44, resf44s
	read(28,*) resf21, resf22, resf23, resf23p
	read(28,*) resf1,resfb
	read(28,*) tot
	close(28)
	xnorm = dfloat(inbevent)/tot
c###	  
c   saves for normalization
        write(12,100) inbevent,tot
100     format(1x,i12,d12.5)
cccc
	if (nlo) then
cccc
c###
	iborn = 0
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	inbrealevent = int(resf34*xnorm)
c        
	open(28,file=path_bsfile(1:ilen)//'threepart1.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : 3 partons : 5//3'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 34
        ntrack = 3
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec,irandom)
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb = dble(wrvec(i1+1))
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f4dimo3,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	inbrealevent = int(resf44*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'threepart2.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : 3 partons : 5//4'
	write(6,*)'First Part'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 44
        ntrack = 3
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec,irandom)
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb = dble(wrvec(i1+1))
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f4dimo4,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
** 	itmx1 = jtmx1
** 	itmx2 = jtmx2
	itmx1 = 5
	itmx2 = 5
	icall = jcall1
	idim = 7
	iwild = 7
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	inbrealevent = int(resf44s*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'threepart2s.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : 3 partons : 5//4'
	write(6,*)'Second Part'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 45
        ntrack = 3
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec,irandom)
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb = dble(wrvec(i1+1))
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f4dimo4s,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	inbrealevent = int(resf21*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'qtwopart1.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : quasi 2 partons : col. ini'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 21
        ntrack = 2
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f2dimo1,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	inbrealevent = int(resf22*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'qtwopart2.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : quasi 2 partons : col. fi3'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 22
        ntrack = 2
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f2dimo2,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	inbrealevent = int(resf23*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'qtwopart3.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : quasi 2 partons : col. fi4'
	write(6,*)'first part'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 23
        ntrack = 2
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f2dimo3,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	inbrealevent = int(resf23p*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'qtwopart3p.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : quasi 2 partons : col. fi4'
	write(6,*)'second part'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 33
        ntrack = 2
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f2dimo3p,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall2
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	inbrealevent = int(resf1*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'twopart1.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : 2 partons : virtuelle'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 10
        ntrack = 2
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f1dimo,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
	open(28,file=path_bsfile(1:ilen)//'efficacity_sb.dat',
     #	status='unknown')
          isum = ipositif+inegatif
	  idiff = ipositif-inegatif
	  write(28,*) 'Before the Born term'
	  write(28,*) 'Nb of positif events:',ipositif
	  write(28,*) 'Nb of negatif events:',inegatif
	  write(28,*) 'Asymetry:',dfloat(idiff)/dfloat(isum)
	close(28)
cccc
	endif
cccc
c###
	if (lo) then
c###
	inbrealevent = int(resfb*xnorm)
c	
	iborn = 1
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall2
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
c
	open(28,file=path_bsfile(1:ilen)//'twopart2.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : 2 partons : born'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 11
        ntrack = 2
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec1,irandom)
	    if (facteur_symetrie) then
              call ranlux(wrvec2,irandom)
	    endif
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  if (facteur_symetrie) then
	    rnumb2 = dble(wrvec2(i1+1))
	    lhalf = (rnumb2.le.half)
	  else
	    lhalf = .true.
	  endif
	  call spring(f1dimo,ixtry)
	  call doubletosingleo(rnumb1,lhalf,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
	endif  
c###
	open(28,file=path_bsfile(1:ilen)//'efficacity_ab.dat',
     #	status='unknown')
          isum = ipositif+inegatif
	  idiff = ipositif-inegatif
	  write(28,*) 'After the Born term'
	  write(28,*) 'Nb of positif events:',ipositif
	  write(28,*) 'Nb of negatif events:',inegatif
	  write(28,*) 'Verif:',isum
	  write(28,*) 'Asymetry:',dfloat(idiff)/dfloat(isum)
	close(28)
cccc
	endif  
cccc
	return
	end
c
c*****************************************************
c  partie ou pt5 < ptm
c*****************************************************
c integrale a 4 dimensions
	double precision function f4dimo3(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical isol3,isol3p,isol4,isol4p
	common/hadron/ih1,ih2,ih3,ih4
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/coup/ptm,r
	common/flag/icoup
	common/coupexp/ptm_exp,r_exp
	common/flagexp/icoup_exp
	common/alem/iloopem
	common/scale/m,mf,mu
	common/faconv/hc2
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/randd/rnumb
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension sf0(k0max),camp0(k0max)
	dimension temp0(k0max),temp1(k0max),spcou(k0max)
	dimension p1(4),p2(4),pp3(4)
	dimension p3(4),p4(4),p5(4)
	dimension p1p(4),p2p(4),p4p(4),p5p(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ x3 -------------------------------------------
	x3max = 1.d0
	x3min_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	if (i_flag_iso.eq.1) then
	  x3min_isol = gpt3/(gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  x3min_isol = 1.d0/(1.d0+etmax)
	endif
	x3minp = dmax1(x3min_incl,x3min_isol)
	x3min = dmin1(x3max,x3minp)
	x3 = x3min + (x3max-x3min)*xx(4)
	pt3 = gpt3/x3
c ------------------ pt5 --------------------------------------------
	pt5max = dsqrt(pt3**2+gptsup**2)
	pt5min = ptm
	pt5 = pt5min + (pt5max-pt5min)*xx(5)
c ------------------ fi35 -------------------------------------------
	x11 = pt3/dsqrt(s)*dexp(y3)+dabs(pt3-pt5)/dsqrt(s)*dexp(y4)
	x21 = pt3/dsqrt(s)*dexp(-y3)+dabs(pt3-pt5)/dsqrt(s)*dexp(-y4)
	if (x11.gt.un.or.x21.gt.un) then
	  f4dimo3 = 0.d0
	  return
	endif
	sy5max = dlog(dsqrt(s)/pt5*(1.d0-x11))
	sy5min = -dlog(dsqrt(s)/pt5*(1.d0-x21))
	if (sy5max.lt.sy5min) then
	  f4dimo3 = 0.d0
	  return
	endif
	if (pt5.le.pt3/2.d0) then
	  ffimax = pi
	else
	  ffimax = dacos(-pt3/(2.d0*pt5))
	endif
	gpt4max = dmin1(pt3+pt5,dsqrt(s)/(2.d0*dcosh(y4))
     #	,gptsup)
	xarg = (gpt4max**2-pt3**2-pt5**2)/
     #	(2.d0*pt3*pt5)
	if (xarg.gt.1.d0) then
	  ffimin = 0.d0
	else if (xarg.lt.-1.d0) then
	  f4dimo3 = 0.d0
	  return
	else
	  ffimin = dacos(xarg)
	endif
	if (y3.le.sy5max.and.y3.ge.sy5min) then
	  omegamin = datan(ffimin/(sy5max-y3))
	  omegamax = pi-datan(ffimin/dabs(sy5min-y3))
	  srmax2 = dmax1((sy5max-y3)**2+ffimax**2,
     #	  (sy5min-y3)**2+ffimax**2)
	  rmax2 = dsqrt(srmax2)
	  rmin2 = ffimin
	else if (y3.lt.sy5max.and.y3.lt.sy5min) then
	  omegamin = datan(ffimin/(sy5max-y3))
	  omegamax = datan(ffimax/(sy5min-y3))
	  srmax2 = (sy5max-y3)**2+ffimax**2
	  rmax2 = dsqrt(srmax2)
	  rmin2 = dsqrt(ffimin**2+(sy5min-y3)**2)
	else if (y3.gt.sy5max.and.y3.gt.sy5min) then
	  omegamin = pi-datan(ffimax/dabs(sy5max-y3))
	  omegamax = pi-datan(ffimin/dabs(sy5min-y3))
	  srmax2 = (sy5min-y3)**2+ffimax**2
	  rmax2 = dsqrt(srmax2)
	  rmin2 = dsqrt(ffimin**2+(sy5max-y3)**2)
	endif
	if (omegamin.le.omegamax) then
	  omega = omegamin + (omegamax-omegamin)*xx(6)
	else	  
	   f4dimo3 = 0.d0
	   return
	endif
	xbig = 1.d+20
	xsmall = 1.d-12
	if (omega.lt.xsmall.or.omega.gt.(pi-xsmall)) then
	  rmax1 = xbig
	else
	  rmax1 = pi/dsin(omega)
	endif
	rrmax = dmin1(rmax1,rmax2)
	rrmin = rmin2
	if (rrmin.le.rrmax) then
	  rr = rrmin + (rrmax-rrmin)*xx(7)
	else	  
	   f4dimo3 = 0.d0
	   return
	endif
	fi35 = rr*dsin(omega)
	if (fi35.gt.pi) then
	  f4dimo3 = 0.d0
	  return
	endif
c ------------------ gpt4 -------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	gpt4 = dsqrt(pt3**2+pt5**2+2.d0*pt3*pt5*dcos(fi35))
	if (gpt4.gt.gpt4max.or.gpt4.lt.gpt4min) then
	  f4dimo3 = 0.d0
	  return
	endif
c ------------------ x -------------------------------------------
	x1c = pt3/dsqrt(s)*dexp(y3) + gpt4/dsqrt(s)*dexp(y4)
	x2c = pt3/dsqrt(s)*dexp(-y3) + gpt4/dsqrt(s)*dexp(-y4)
	if (x1c.ge.un.or.x2c.ge.un) then
	  f4dimo3 = 0.d0
	  return
	endif
c ------------------ y5 -------------------------------------------
	y5max = dlog(dsqrt(s)/pt5*(1.d0-x1c))
	y5min = -dlog(dsqrt(s)/pt5*(1.d0-x2c))
	y5 = y3 + rr*dcos(omega)
	if (y5.gt.y5max.or.y5.lt.y5min) then
	  f4dimo3 = 0.d0
	  return
	endif
c ------------------ x -------------------------------------------
	x1 = x1c + pt5/dsqrt(s)*dexp(y5)
	x2 = x2c + pt5/dsqrt(s)*dexp(-y5)
	if (x1.ge.un.or.x2.ge.un) then
	  f4dimo3 = 0.d0
	  return
	endif
c coupure en angle-------------------------------------------------
        r35s = (y3-y5)**2+fi35**2
        rs = r**2
	icorr = 0
        if (r35s.lt.rs) then
           f4dimo3 = 0.d0
           return
**           icorr = 1
        endif
c coupure en isolement-------------------------------------------------
	call isolement(y3,y5,fi35,x3,pt5,gpt3,r_isol,etmax,
     #	isol3)
	cos_fi45 = -(pt3*dcos(fi35)+pt5)/gpt4
	fi45 = zacos(cos_fi45)
	call isolement(y4,y5,fi45,un,pt5,gpt4,r_isol,etmax,
     #	isol4)
** 	if (.not.(isol3)) then
** 	  f4dimo3 = 0.d0
** 	  return
** 	endif
	xisol4 = 1.d0
	if (.not.(isol4).or..not.(isol3)) then
	  if (icorr.eq.0) then
	    f4dimo3 = 0.d0
	    return
	  else if (icorr.eq.1) then
	    xisol4 = 0.d0
	  endif
	endif
c ------------------ fin isolement ------------------------------------
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)*
     #	rr*(rrmax-rrmin)*(omegamax-omegamin)*(x3max-x3min)*
     #	(pt5max-pt5min)
	x1p = (pt3+pt5)/dsqrt(s)*(dexp(y3)+dexp(y4))
	x2p = (pt3+pt5)/dsqrt(s)*(dexp(-y3)+dexp(-y4))
c ------------------------------------------- 
	p1(1) = x1*dsqrt(s)/2.d0
	p1(2) = 0.d0
	p1(3) = 0.d0
	p1(4) = x1*dsqrt(s)/2.d0
c ------------------------------------------- 
	p2(1) = x2*dsqrt(s)/2.d0
	p2(2) = 0.d0
	p2(3) = 0.d0
	p2(4) = -x2*dsqrt(s)/2.d0
c ------------------------------------------- 
	pp3(1) = pt3*dcosh(y3)
	pp3(2) = pt3
	pp3(3) = 0.d0
	pp3(4) = pt3*dsinh(y3)
c ------------------------------------------- 
	p3(1) = gpt3*dcosh(y3)
	p3(2) = gpt3
	p3(3) = 0.d0
	p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	p5(1) = pt5*dcosh(y5)
	p5(2) = pt5*dcos(fi35)
	p5(3) = pt5*dsin(fi35)
	p5(4) = pt5*dsinh(y5)
c ------------------------------------------- 
	p4(1) = p1(1) + p2(1) - pp3(1) - p5(1)
	p4(2) = p1(2) + p2(2) - pp3(2) - p5(2)
	p4(3) = p1(3) + p2(3) - pp3(3) - p5(3)
	p4(4) = p1(4) + p2(4) - pp3(4) - p5(4)
c -------------------------------------------
	p1p(1) = x1p*dsqrt(s)/2.d0
	p1p(2) = 0.d0
	p1p(3) = 0.d0
	p1p(4) = x1p*dsqrt(s)/2.d0
c ------------------------------------------- 
	p2p(1) = x2p*dsqrt(s)/2.d0
	p2p(2) = 0.d0
	p2p(3) = 0.d0
	p2p(4) = -x2p*dsqrt(s)/2.d0
c ------------------------------------------- 
	p5p(1) = pt5*dcosh(y3)
	p5p(2) = pt5
	p5p(3) = 0.d0
	p5p(4) = pt5*dsinh(y3)
c ------------------------------------------- 
	p4p(1) = p1p(1) + p2p(1) - pp3(1) - p5p(1)
	p4p(2) = p1p(2) + p2p(2) - pp3(2) - p5p(2)
	p4p(3) = p1p(3) + p2p(3) - pp3(3) - p5p(3)
	p4p(4) = p1p(4) + p2p(4) - pp3(4) - p5p(4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	mass_cut = 1.d0
	if (iphocut.eq.0) then
c angle azymuthal entre les deux photons plus grand que fimin
	  abspt3 = dsqrt(p3(2)**2+p3(3)**2)
	  abspt4 = dsqrt(p4(2)**2+p4(3)**2)
	  fi34 = zacos((p3(2)*p4(2)+p3(3)*p4(3))/(abspt3*abspt4))
	  if (fi34.le.fimin) then
	    f4dimo3 = 0.d0
	    return
	  endif
	else if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    if (icorr.eq.0) then
	      f4dimo3 = 0.d0
	      return
	    else if (icorr.eq.1) then
	      mass_cut = 0.d0
	    endif
	  endif
	else if (iphocut.eq.2) then
c distance en pseudo-rapidite angle azymuthal entre les 2 photons plus
c grande que dmin
	  dp3p4 = p3(2)*p4(2)+p3(3)*p4(3)
	  dp3p3 = p3(2)*p3(2)+p3(3)*p3(3)
	  dp4p4 = p4(2)*p4(2)+p4(3)*p4(3)
	  fi34 = zacos(dp3p4/dsqrt(dp3p3*dp4p4))
	  dphopho2 = (y3-y4)**2+fi34**2
	  dphopho = dsqrt(dphopho2)	  	  
	  if (dphopho.lt.dmin) then
	    f4dimo3 = 0.d0
	    return
	  endif
	endif	
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	call fspcouo(n,spcou)
c
	alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s^2*alpha/(4*c_i*c_j*pi*s^2)
	cij = alfas(iloop,mu*mu)**2*alpha_em/(pi*s*s)
c
	call strfrao(x1,ih1,x2,ih2,x3,ih3,sf0)
c	call ampo(3,p1,p2,pp3,p4,p5,camp0)
	call ampo3(s,x1,x2,y3,pt3,y5,pt5,fi35,camp0)
	call vec_dmult_vector(camp0,sf0,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f4dimo3 = 0.d0
	do i = j_processo_min,j_processo_max
	  f4dimo3 = temp1(i) + temp1(i+18) + f4dimo3
	enddo
	f4dimo3 = cij*f4dimo3
	if (icorr.eq.1) then
	  mass_cutp = 1.d0
	  if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	    mgg2 = 2.d0*sca(p3,p4p)
	    mgg = dsqrt(mgg2)
	    if (mgg.lt.mmin.or.mgg.gt.mmax) then
	      mass_cutp = 0.d0
	    endif
	  endif	
c coupure en isolement-------------------------------------------------
	  call isolement(y3,y3,0.d0,x3,pt5,gpt3,r_isol,etmax,
     #	  isol3p)
	  call isolement(y4,y3,pi,un,pt5,gpt4,r_isol,etmax,
     #	  isol4p)
	  xisol4p = 1.d0
	  if (.not.(isol4p).or..not.(isol3p)) then
	    xisol4p = 0.d0
	  endif
c ------------------ fin isolement ------------------------------------
	  call choiscale(p3,p4p,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	  alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s^2*alpha/(4*c_i*c_j*pi*s^2)
	  cij = alfas(iloop,mu*mu)**2*alpha_em/(pi*s*s)
c
	  call strfrao(x1p,ih1,x2p,ih2,x3,ih3,sf0)
	  call ampo_corr3(s,x1p,x2p,y3,pt3,pt5,r35s,camp0)
	  call vec_dmult_vector(camp0,sf0,k0max,temp0)
	  call vec_dmult_vector(temp0,spcou,k0max,temp1)
	  f4dimo3p = 0.d0
	  do i = j_processo_min,j_processo_max
	    f4dimo3p = temp1(i) + temp1(i+18) + f4dimo3p
	  enddo
	  f4dimo3p = cij*f4dimo3p
	  f4dimo3 = f4dimo3*mass_cut*xisol4 - 
     #	  f4dimo3p*mass_cutp*xisol4p
	endif
c
	f4dimo3 = f4dimo3 * xjacob * (pt3/x3*pt5) * hc2
	if (ispring.eq.1) then
	  half = 0.5d0
	  if (rnumb.le.half) then
	    xfi12 = 2.d0*pi-zacos(p4(2)/gpt4)
	    xfi13 = fi35
	  else
	    xfi12 = zacos(p4(2)/gpt4)
	    xfi13 = 2.d0*pi-fi35
	  endif
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = pt5
	  xy3 = y5
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f4dimo3
	  f4dimo3 = dabs(f4dimo3)
	endif
c
	return
	end 
c********************************************************************
c integrale a 4 dimensions
	double precision function f4dimo4(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical isol3,isol3p,isol4,isol4p
	common/hadron/ih1,ih2,ih3,ih4
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/coup/ptm,r
	common/flag/icoup
	common/coupexp/ptm_exp,r_exp
	common/flagexp/icoup_exp
	common/alem/iloopem
	common/scale/m,mf,mu
	common/faconv/hc2
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/randd/rnumb
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension sf0(k0max),camp0(k0max)
	dimension temp0(k0max),temp1(k0max),spcou(k0max)
	dimension p1(4),p2(4),pp3(4)
	dimension p3(4),p4(4),p5(4)
	dimension p1p(4),p2p(4),p4p(4),p5p(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ x3 ---------------------------------------------
	x3max = 1.d0
	x3min_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	if (i_flag_iso.eq.1) then
	  x3min_isol = gpt3/(gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  x3min_isol = 1.d0/(1.d0+etmax)
	endif
	x3minp = dmax1(x3min_incl,x3min_isol)
** 	x3minp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3min = dmin1(x3max,x3minp)
	x3 = x3min + (x3max-x3min)*xx(4)
	pt3 = gpt3/x3
c ------------------ pt5 --------------------------------------------
	pt5max = dsqrt(pt3**2+gptsup**2)
	pt5min = ptm
	pt5 = pt5min + (pt5max-pt5min)*xx(5)
c ------------------ fi45 -------------------------------------------
	x11 = pt3/dsqrt(s)*dexp(y3)+dabs(pt3-pt5)/dsqrt(s)*dexp(y4)
	x21 = pt3/dsqrt(s)*dexp(-y3)+dabs(pt3-pt5)/dsqrt(s)*dexp(-y4)
	if (x11.gt.un.or.x21.gt.un) then
	  f4dimo4 = 0.d0
	  return
	endif
	sy5max = dlog(dsqrt(s)/pt5*(1.d0-x11))
	sy5min = -dlog(dsqrt(s)/pt5*(1.d0-x21))
	if (sy5max.lt.sy5min) then
	  f4dimo4 = 0.d0
	  return
	endif
	if (y4.le.sy5max.and.y4.ge.sy5min) then
	  omegamin = 0.d0
	  omegamax = pi
	  srmax2 = dmax1((sy5max-y4)**2+pi**2,(sy5min-y4)**2+pi**2)
	  rmax2 = dsqrt(srmax2)
	else if (y4.lt.sy5max.and.y4.lt.sy5min) then
	  omegamin = 0.d0
	  omegamax = pi/2.d0
	  srmax2 = (sy5max-y4)**2+pi**2
	  rmax2 = dsqrt(srmax2)
	else if (y4.gt.sy5max.and.y4.gt.sy5min) then
	  omegamin = pi/2.d0
	  omegamax = pi
	  srmax2 = (sy5min-y4)**2+pi**2
	  rmax2 = dsqrt(srmax2)
	endif
	omega = omegamin + (omegamax-omegamin)*xx(6)
	xbig = 1.d+20
	xsmall = 1.d-12
	if (omega.lt.xsmall.or.omega.gt.(pi-xsmall)) then
	  rmax1 = xbig
	else
	  rmax1 = pi/dsin(omega)
	endif
	rrmax = dmin1(rmax1,rmax2)
cc	rrmin = r
	rrmin = 0.d0
	rr = rrmin + (rrmax-rrmin)*xx(7)
	fi45 = rr*dsin(omega)
	if (fi45.gt.pi) then
	   f4dimo4 = 0.d0
	   return
	endif
c ------------------ gpt4 ------------------------------------------
	gpt4max = dmin1(dsqrt(s)/(2.d0*dcosh(y4))
     #	,gptsup)
	gpt4min = gptinf
	if (pt5.le.pt3) then
	  gpt4 = -pt5*dcos(fi45)+dsqrt(pt3**2-(pt5*dsin(fi45))**2)
	else if (pt5.gt.pt3) then
	  if ((fi45.ge.pi/2.d0.and.fi45.le.pi).and.
     #	  ((pt3-pt5*dsin(fi45)).ge.0.d0)) then
	    gpt4 = -pt5*dcos(fi45)+dsqrt(pt3**2-(pt5*dsin(fi45))**2)
	  else
	    f4dimo4 = 0.d0
	    return
	  endif
	endif
	if (gpt4.gt.gpt4max.or.gpt4.lt.gpt4min) then
	  f4dimo4 = 0.d0
	  return
	endif
	pt4 = gpt4
c ------------------ x -------------------------------------------
	x1c = pt3/dsqrt(s)*dexp(y3) + gpt4/dsqrt(s)*dexp(y4)
	x2c = pt3/dsqrt(s)*dexp(-y3) + gpt4/dsqrt(s)*dexp(-y4)
	if (x1c.ge.un.or.x2c.ge.un) then
	  f4dimo4 = 0.d0
	  return
	endif
c ------------------ y5 -------------------------------------------
	y5max = dlog(dsqrt(s)/pt5*(1.d0-x1c))
	y5min = -dlog(dsqrt(s)/pt5*(1.d0-x2c))
	y5 = y4 + rr*dcos(omega)
	if (y5.gt.y5max.or.y5.lt.y5min) then
	   f4dimo4 = 0.d0
	   return
	endif
	if (gpt4.lt.0.d0) then
	  write (*,*) 'q1:',pt5,pt3,fi45,gpt4
	endif
c ------------------ x -------------------------------------------
	x1 = x1c + pt5/dsqrt(s)*dexp(y5)
	x2 = x2c + pt5/dsqrt(s)*dexp(-y5)
	if (x1.ge.un.or.x2.ge.un) then
	  f4dimo4 = 0.d0
	  return
	endif
c coupure en angle-------------------------------------------------
        r45s = (y4-y5)**2+fi45**2
        rs = r**2
	icorr = 0
        if (r45s.lt.rs) then
           f4dimo4 = 0.d0
           return
**           icorr = 1
        endif
c coupure en isolement-------------------------------------------------
	cos_fi35 = -(gpt4*dcos(fi45)+pt5)/pt3
	fi35 = zacos(cos_fi35)
	call isolement(y3,y5,fi35,x3,pt5,gpt3,r_isol,etmax,
     #	isol3)
	call isolement(y4,y5,fi45,un,pt5,gpt4,r_isol,etmax,
     #	isol4)
	xisol3 = 1.d0
	if (.not.(isol3).or..not.(isol4)) then
	  if (icorr.eq.0) then
	    f4dimo4 = 0.d0
	    return
	  else if (icorr.eq.1) then
	    xisol3 = 0.d0
	  endif
	endif
c ------------------ fin isolement ------------------------------------
	xjacob = (y3max-y3min)*(y4max-y4min)*(x3max-x3min)*
     #	rr*(rrmax-rrmin)*(omegamax-omegamin)*(gpt3max-gpt3min)*
     #	(pt5max-pt5min)
	x1s = pt3/dsqrt(s)*(dexp(y3)+dexp(y4))
	x2s = pt3/dsqrt(s)*(dexp(-y3)+dexp(-y4))
	gpt4p = (pt3-pt5)
c ------------------------------------------- 
	p1(1) = x1*dsqrt(s)/2.d0
	p1(2) = 0.d0
	p1(3) = 0.d0
	p1(4) = x1*dsqrt(s)/2.d0
c ------------------------------------------- 
	p2(1) = x2*dsqrt(s)/2.d0
	p2(2) = 0.d0
	p2(3) = 0.d0
	p2(4) = -x2*dsqrt(s)/2.d0
c ------------------------------------------- 
	p3(1) = gpt3*dcosh(y3)
	p3(2) = gpt3
	p3(3) = 0.d0
	p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	pp3(1) = pt3*dcosh(y3)
	pp3(2) = pt3
	pp3(3) = 0.d0
	pp3(4) = pt3*dsinh(y3)
c ------------------------------------------- 
c attention fi45 est l'angle entre 4 et 5
	p5(1) = pt5*dcosh(y5)
** 	p5(2) = pt5*dcos(fi45)
** 	p5(3) = pt5*dsin(fi45)
	p5(2) = -pt5*(gpt4*dcos(fi45)+pt5)/pt3
	p5(3) = -pt5*gpt4*dsin(fi45)/pt3
	p5(4) = pt5*dsinh(y5)
c ------------------------------------------- 
	p4(1) = p1(1) + p2(1) - pp3(1) - p5(1)
	p4(2) = p1(2) + p2(2) - pp3(2) - p5(2)
	p4(3) = p1(3) + p2(3) - pp3(3) - p5(3)
	p4(4) = p1(4) + p2(4) - pp3(4) - p5(4)
c -------------------------------------------
	p1p(1) = x1s*dsqrt(s)/2.d0
	p1p(2) = 0.d0
	p1p(3) = 0.d0
	p1p(4) = x1s*dsqrt(s)/2.d0
c ------------------------------------------- 
	p2p(1) = x2s*dsqrt(s)/2.d0
	p2p(2) = 0.d0
	p2p(3) = 0.d0
	p2p(4) = -x2s*dsqrt(s)/2.d0
c ------------------------------------------- 
	p5p(1) = pt5*dcosh(y4)
	p5p(2) = -pt5
	p5p(3) = 0.d0
	p5p(4) = pt5*dsinh(y4)
c ------------------------------------------- 
	p4p(1) = p1p(1) + p2p(1) - pp3(1) - p5p(1)
	p4p(2) = p1p(2) + p2p(2) - pp3(2) - p5p(2)
	p4p(3) = p1p(3) + p2p(3) - pp3(3) - p5p(3)
	p4p(4) = p1p(4) + p2p(4) - pp3(4) - p5p(4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	mass_cut = 1.d0
	if (iphocut.eq.0) then
c angle azymuthal entre les deux photons plus grand que fimin
	  abspt3 = dsqrt(p3(2)**2+p3(3)**2)
	  abspt4 = dsqrt(p4(2)**2+p4(3)**2)
	  fi34 = zacos((p3(2)*p4(2)+p3(3)*p4(3))/(abspt3*abspt4))
	  if (fi34.le.fimin) then
	    f4dimo4 = 0.d0
	    return
	  endif
	else if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    if (icorr.eq.0) then
	      f4dimo4 = 0.d0
	      return
	    else if (icorr.eq.1) then
	      mass_cut = 0.d0
	    endif
	  endif
	else if (iphocut.eq.2) then
c distance en pseudo-rapidite angle azymuthal entre les 2 photons plus
c grande que dmin
	  dp3p4 = p3(2)*p4(2)+p3(3)*p4(3)
	  dp3p3 = p3(2)*p3(2)+p3(3)*p3(3)
	  dp4p4 = p4(2)*p4(2)+p4(3)*p4(3)
	  fi34 = zacos(dp3p4/dsqrt(dp3p3*dp4p4))
	  dphopho2 = (y3-y4)**2+fi34**2
	  dphopho = dsqrt(dphopho2)	  	  
	  if (dphopho.lt.dmin) then
	    f4dimo4 = 0.d0
	    return
	  endif
	endif	
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	call fspcouo(n,spcou)
c
	alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s^2*alpha/(4*c_i*c_j*pi*s^2)
	cij = alfas(iloop,mu*mu)**2*alpha_em/(pi*s*s)
c
	call strfrao(x1,ih1,x2,ih2,x3,ih3,sf0)
cc	call ampo(4,p1,p2,pp3,p4,p5,camp0)
	call ampo4(s,x1,x2,y4,gpt4,y5,pt5,fi45,camp0)
	call vec_dmult_vector(camp0,sf0,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f4dimo4 = 0.d0
	do i = j_processo_min,j_processo_max
	  f4dimo4 = temp1(i) + temp1(i+18) + f4dimo4
	enddo
	f4dimo4 = f4dimo4*cij*gpt4/dsqrt(pt3**2-(pt5*dsin(fi45))**2)
	if (icorr.eq.1) then
	  mass_cutp = 1.d0
	  if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	    mgg2 = 2.d0*sca(p3,p4p)
	    mgg = dsqrt(mgg2)
	    if (mgg.lt.mmin.or.mgg.gt.mmax) then
	      mass_cutp = 0.d0
	    endif
	  endif	
c coupure en isolement-------------------------------------------------
	  call isolement(y3,y4,pi,x3,pt5,gpt3,r_isol,etmax,
     #	  isol3p)
	  call isolement(y4,y4,0.d0,un,pt5,pt3-pt5,r_isol,etmax,
     #	  isol4p)
	  xisol3p = 1.d0
	  if (.not.(isol3p).or..not.(isol4p)) then
	    xisol3p = 0.d0
	  endif
c ------------------ fin isolement ------------------------------------
	  call choiscale(p3,p4p,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	  alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s^2*alpha/(4*c_i*c_j*pi*s^2)
	  cij = alfas(iloop,mu*mu)**2*alpha_em/(pi*s*s)
	  cij = 1.d0/(pi*s*s)
c
	  call strfrao(x1s,ih1,x2s,ih2,x3,ih3,sf0)
	  call ampo_corr4(s,x1s,x2s,y4,pt3-pt5,pt5,r45s,camp0)
	  call vec_dmult_vector(camp0,sf0,k0max,temp0)
	  call vec_dmult_vector(temp0,spcou,k0max,temp1)
	  f4dimo4p = 0.d0
	  do i = j_processo_min,j_processo_max
	    f4dimo4p = temp1(i) + temp1(i+18) + f4dimo4p
	  enddo
	  f4dimo4p = f4dimo4p*cij*(pt3-pt5)/pt3
	  f4dimo4 = f4dimo4*mass_cut*xisol3 - 
     #	  f4dimo4p*mass_cutp*xisol3p
	endif
c
	f4dimo4 = f4dimo4 * xjacob * (pt5/x3) * hc2
     #            *pt3
	if (ispring.eq.1) then
	  half = 0.5d0
	  if (rnumb.le.half) then
	    xfi12 = 2.d0*pi-zacos(p4(2)/gpt4)
	    xfi13 = zacos(-(gpt4*dcos(fi45)+pt5)/pt3)
	  else
	    xfi12 = zacos(p4(2)/gpt4)
	    xfi13 = 2.d0*pi-zacos(-(gpt4*dcos(fi45)+pt5)/pt3)
	  endif
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = pt5
	  xy3 = y5
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f4dimo4
	  f4dimo4 = dabs(f4dimo4)
	endif
c
	return
	end 
c********************************************************************
c integrale a 4 dimensions
	double precision function f4dimo4s(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical isol3,isol3p,isol4,isol4p
	common/hadron/ih1,ih2,ih3,ih4
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/coup/ptm,r
	common/flag/icoup
	common/coupexp/ptm_exp,r_exp
	common/flagexp/icoup_exp
	common/alem/iloopem
	common/scale/m,mf,mu
	common/faconv/hc2
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/randd/rnumb
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension sf0(k0max),camp0(k0max)
	dimension temp0(k0max),temp1(k0max),spcou(k0max)
	dimension p1(4),p2(4),pp3(4)
	dimension p3(4),p4(4),p5(4)
	dimension p1p(4),p2p(4),p4p(4),p5p(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ x3 ---------------------------------------------
	x3max = 1.d0
	x3min_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	if (i_flag_iso.eq.1) then
	  x3min_isol = gpt3/(gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  x3min_isol = 1.d0/(1.d0+etmax)
	endif
	x3minp = dmax1(x3min_incl,x3min_isol)
** 	x3minp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3min = dmin1(x3max,x3minp)
	x3 = x3min + (x3max-x3min)*xx(4)
	pt3 = gpt3/x3
c ------------------ fi45 -------------------------------------------
	fi45min = pi/2.d0
	fi45max = pi
	fi45 = fi45min + (fi45max-fi45min)*xx(6)
c ------------------ pt5 --------------------------------------------
	xbig = 1.d+20
	xsmall = 1.d-12
	if (fi45.lt.xsmall.or.fi45.gt.(pi-xsmall)) then
	  pt5max1 = xbig
	else
	  pt5max1 = pt3/dsin(fi45)
	endif
	pt5max = dmin1(dsqrt(s)/2.d0,pt5max1)
	pt5min = pt3
	pt5 = pt5min + (pt5max-pt5min)*xx(5)
c ------------------ gpt4 ------------------------------------------
	gpt4max = dmin1(dsqrt(s)/(2.d0*dcosh(y4))
     #	,gptsup)
	gpt4min = gptinf
	if ((pt3-pt5*dsin(fi45)).ge.0.d0) then
	  gpt4 = -pt5*dcos(fi45)-dsqrt(pt3**2-(pt5*dsin(fi45))**2)
	else
	  f4dimo4s = 0.d0
	  return
	endif
	if (gpt4.gt.gpt4max.or.gpt4.lt.gpt4min) then
	  f4dimo4s = 0.d0
	  return
	endif
	pt4 = gpt4
c ------------------ x -------------------------------------------
	x1c = pt3/dsqrt(s)*dexp(y3) + gpt4/dsqrt(s)*dexp(y4)
	x2c = pt3/dsqrt(s)*dexp(-y3) + gpt4/dsqrt(s)*dexp(-y4)
	if (x1c.ge.un.or.x2c.ge.un) then
	  f4dimo4s = 0.d0
	  return
	endif
c ------------------ y5 -------------------------------------------
	y5max = dlog(dsqrt(s)/pt5*(1.d0-x1c))
	y5min = -dlog(dsqrt(s)/pt5*(1.d0-x2c))
	y5 = y5min + (y5max-y5min)*xx(7)
c ------------------ x -------------------------------------------
	x1 = x1c + pt5/dsqrt(s)*dexp(y5)
	x2 = x2c + pt5/dsqrt(s)*dexp(-y5)
	if (x1.ge.un.or.x2.ge.un) then
	  f4dimo4s = 0.d0
	  return
	endif
c coupure en angle-------------------------------------------------
        r45s = (y4-y5)**2+fi45**2
        rs = r**2
	icorr = 0
        if (r45s.lt.rs) then
           f4dimo4s = 0.d0
           return
**            icorr = 1
        endif
c coupure en isolement-------------------------------------------------
	cos_fi35 = -(gpt4*dcos(fi45)+pt5)/pt3
	fi35 = zacos(cos_fi35)
	call isolement(y3,y5,fi35,x3,pt5,gpt3,r_isol,etmax,
     #	isol3)
	call isolement(y4,y5,fi45,un,pt5,gpt4,r_isol,etmax,
     #	isol4)
	xisol3 = 1.d0
	if (.not.(isol3).or..not.(isol4)) then
	  if (icorr.eq.0) then
	    f4dimo4s = 0.d0
	    return
	  else if (icorr.eq.1) then
	    xisol3 = 0.d0
	  endif
	endif
c ------------------ fin isolement ------------------------------------
	xjacob = (y3max-y3min)*(y4max-y4min)*(x3max-x3min)*
     #	(y5max-y5min)*(fi45max-fi45min)*(gpt3max-gpt3min)*
     #	(pt5max-pt5min)
	x1s = pt3/dsqrt(s)*(dexp(y3)+dexp(y4))
	x2s = pt3/dsqrt(s)*(dexp(-y3)+dexp(-y4))
	gpt4p = (pt3-pt5)
c ------------------------------------------- 
	p1(1) = x1*dsqrt(s)/2.d0
	p1(2) = 0.d0
	p1(3) = 0.d0
	p1(4) = x1*dsqrt(s)/2.d0
c ------------------------------------------- 
	p2(1) = x2*dsqrt(s)/2.d0
	p2(2) = 0.d0
	p2(3) = 0.d0
	p2(4) = -x2*dsqrt(s)/2.d0
c ------------------------------------------- 
	p3(1) = gpt3*dcosh(y3)
	p3(2) = gpt3
	p3(3) = 0.d0
	p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	pp3(1) = pt3*dcosh(y3)
	pp3(2) = pt3
	pp3(3) = 0.d0
	pp3(4) = pt3*dsinh(y3)
c ------------------------------------------- 
c attention fi45 est l'angle entre 4 et 5
	p5(1) = pt5*dcosh(y5)
** 	p5(2) = pt5*dcos(fi45)
** 	p5(3) = pt5*dsin(fi45)
	p5(2) = -pt5*(gpt4*dcos(fi45)+pt5)/pt3
	p5(3) = -pt5*gpt4*dsin(fi45)/pt3
	p5(4) = pt5*dsinh(y5)
c ------------------------------------------- 
	p4(1) = p1(1) + p2(1) - pp3(1) - p5(1)
	p4(2) = p1(2) + p2(2) - pp3(2) - p5(2)
	p4(3) = p1(3) + p2(3) - pp3(3) - p5(3)
	p4(4) = p1(4) + p2(4) - pp3(4) - p5(4)
c -------------------------------------------
	p1p(1) = x1s*dsqrt(s)/2.d0
	p1p(2) = 0.d0
	p1p(3) = 0.d0
	p1p(4) = x1s*dsqrt(s)/2.d0
c ------------------------------------------- 
	p2p(1) = x2s*dsqrt(s)/2.d0
	p2p(2) = 0.d0
	p2p(3) = 0.d0
	p2p(4) = -x2s*dsqrt(s)/2.d0
c ------------------------------------------- 
	p5p(1) = pt5*dcosh(y4)
	p5p(2) = -pt5
	p5p(3) = 0.d0
	p5p(4) = pt5*dsinh(y4)
c ------------------------------------------- 
	p4p(1) = p1p(1) + p2p(1) - pp3(1) - p5p(1)
	p4p(2) = p1p(2) + p2p(2) - pp3(2) - p5p(2)
	p4p(3) = p1p(3) + p2p(3) - pp3(3) - p5p(3)
	p4p(4) = p1p(4) + p2p(4) - pp3(4) - p5p(4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	mass_cut = 1.d0
	if (iphocut.eq.0) then
c angle azymuthal entre les deux photons plus grand que fimin
	  abspt3 = dsqrt(p3(2)**2+p3(3)**2)
	  abspt4 = dsqrt(p4(2)**2+p4(3)**2)
	  fi34 = zacos((p3(2)*p4(2)+p3(3)*p4(3))/(abspt3*abspt4))
	  if (fi34.le.fimin) then
	    f4dimo4s = 0.d0
	    return
	  endif
	else if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    if (icorr.eq.0) then
	      f4dimo4s = 0.d0
	      return
	    else if (icorr.eq.1) then
	      mass_cut = 0.d0
	    endif
	  endif
	else if (iphocut.eq.2) then
c distance en pseudo-rapidite angle azymuthal entre les 2 photons plus
c grande que dmin
	  dp3p4 = p3(2)*p4(2)+p3(3)*p4(3)
	  dp3p3 = p3(2)*p3(2)+p3(3)*p3(3)
	  dp4p4 = p4(2)*p4(2)+p4(3)*p4(3)
	  fi34 = zacos(dp3p4/dsqrt(dp3p3*dp4p4))
	  dphopho2 = (y3-y4)**2+fi34**2
	  dphopho = dsqrt(dphopho2)	  	  
	  if (dphopho.lt.dmin) then
	    f4dimo4s = 0.d0
	    return
	  endif
	endif	
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	call fspcouo(n,spcou)
c
	alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s^2*alpha/(4*c_i*c_j*pi*s^2)
	cij = alfas(iloop,mu*mu)**2*alpha_em/(pi*s*s)
c
	call strfrao(x1,ih1,x2,ih2,x3,ih3,sf0)
cc	call ampo(4,p1,p2,pp3,p4,p5,camp0)
	call ampo4(s,x1,x2,y4,gpt4,y5,pt5,fi45,camp0)
	call vec_dmult_vector(camp0,sf0,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f4dimo4s = 0.d0
	do i = j_processo_min,j_processo_max
	  f4dimo4s = temp1(i) + temp1(i+18) + f4dimo4s
	enddo
	f4dimo4s = f4dimo4s*cij*gpt4/dsqrt(pt3**2-(pt5*dsin(fi45))**2)
c
	f4dimo4s = f4dimo4s * xjacob * (pt5/x3) * hc2
     #            *pt3
	if (ispring.eq.1) then
	  half = 0.5d0
	  if (rnumb.le.half) then
	    xfi12 = 2.d0*pi-zacos(p4(2)/gpt4)
	    xfi13 = zacos(-(gpt4*dcos(fi45)+pt5)/pt3)
	  else
	    xfi12 = zacos(p4(2)/gpt4)
	    xfi13 = 2.d0*pi-zacos(-(gpt4*dcos(fi45)+pt5)/pt3)
	  endif
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = pt5
	  xy3 = y5
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f4dimo4s
	  f4dimo4s = dabs(f4dimo4s)
	endif
c
	return
	end 
c*****************************************************
c integrale a 2 dimensions
	double precision function f2dimo1(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/alem/iloopem
	common/faconv/hc2
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension sv2(k0max),sy3(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3))
     #	,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ gpt4 ------------------------------------------
c x3 = gpt3/gpt4 et x3min=2*gpt3/dsqrt(s)*dcosh(y3) <= gpt3/gpt4
c d'ou gpt4 <= dsqrt(s)/(2.d0*dcosh(y3))
c pour les coupures d'isolement on a que
c gpt3 <= gpt4 <= gpt3+etmax car
c x3=gpt3/gpt4 >= gpt3/(gpt3+etmax)
	if (i_flag_iso.eq.1) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3*(1.d0+etmax))
	endif
	gpt4max = dmax1(gpt3,gpt4maxp)
	gpt4min = gpt3
	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(4)
c ------------------ u -------------------------------------------
	u2 = xx(5)
	x3 = gpt3/gpt4
	x10 = gpt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = gpt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimo1 = 0.d0
	  return
	endif
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)
     #	*(gpt4max-gpt4min)
	call fspcouo(n,spcou)
c ------------------------------------------- 
	  p3(1) = gpt3*dcosh(y3)
	  p3(2) = gpt3
	  p3(3) = 0.d0
	  p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	  p5(1) = 0.d0
	  p5(2) = 0.d0
	  p5(3) = 0.d0
	  p5(4) = 0.d0
c ------------------------------------------- 
	  p4(1) = gpt4*dcosh(y4)
	  p4(2) = -gpt4
	  p4(3) = 0.d0
	  p4(4) = gpt4*dsinh(y4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimo1 = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alpha_em = alphaem(iloopem,m*m)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s*alpha/(4*c_i*c_j*s**2)
	cbij = alphas*alpha_em*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	call ssy3o(s,gpt3,y3,y4,m,x10,x20,x3,u2,ih1,ih2,ih3,sy3)
	call ssv2o(s,gpt3,y3,y4,x10,x20,x3,u2,ih1,ih2,ih3,sv2)
	call vec_dadd_vector(sy3,sv2,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dim = 0.d0
	do i = j_processo_min,j_processo_max
	  f2dim = temp1(i) + temp1(i+18) + f2dim
	enddo
cc	f2dim = vec_dsum(temp0,k0max)
c
	f2dimo1 = f2dim * xjacob * alpi * cbij * hc2
	if (ispring.eq.1) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f2dimo1
	  f2dimo1 = dabs(f2dimo1)
	endif
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimo2(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/alem/iloopem
	common/faconv/hc2
	common/flagexp/icoup_exp
	common/coupexp/ptm_exp,r_exp
	common/coup/ptmm,r
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension sy3r(k0max),sy3p(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	zero = 0.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3))
     #	,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ gpt4 ------------------------------------------
	if (i_flag_iso.eq.1) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3*(1.d0+etmax))
	endif
	gpt4max = dmax1(gpt3,gpt4maxp)
	gpt4min = gpt3
	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(4)
c ------------------ x3 -------------------------------------------
	x3max = 1.d0
	x3minpp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
c
	x10 = gpt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = gpt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimo2 =0.d0
	  return
	endif
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
	z3max0 = 1.d0
	z3min0 = gpt3/gpt4
	z3min1 = (gpt4-ptm_exp)/gpt4	
	z3sup = (gpt4-ptm)/gpt4
	if (icoup_exp.eq.0) then
	  z3min = z3min0	
	  z3max = z3max0	
	else if (icoup_exp.eq.1) then
	  z3min = dmax1(z3min0,z3min1)	
	  z3max = z3sup	
	endif
c juste pour tester (en principe inutile)
	if (z3min.gt.z3max) then
	  f2dimo2 = 0.d0
	  write (18,*) 'attention z3min > z3max',z3min,z3max
	  return
	endif
	z3 = z3min + (z3max-z3min)*xx(5)
	x3 = gpt3/(z3*gpt4)
c juste pour tester (en principe inutile)
	if (x3.gt.x3max.or.x3.lt.x3minpp) then
	  f2dimo2 = 0.d0
	  write (18,*) 'attention x3 > x3min ou x3 < x3max',
     #	              x3,x3minpp,x3max
	  return
	endif
c juste pour tester (en principe inutile)
	if (i_flag_iso.eq.1) then
	  z3z = gpt3/((gpt3+etmax)*x3)
	else if (i_flag_iso.eq.2) then
	  z3z = 1.d0/((1.d0+etmax)*x3)
	endif
	if (z3.lt.z3z) then
	  f2dimo2 = 0.d0
	  write (18,*) 'attention z3 < z3z',z3,z3z
	  write (18,*) x3*z3, gpt3/(gpt3+etmax), 1.d0/(1.d0+etmax)
	  return
	endif
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)
     #	*(gpt4max-gpt4min)*(z3max-z3min)
	call fspcouo(n,spcou)
c ------------------------------------------- 
	  p3(1) = gpt3*dcosh(y3)
	  p3(2) = gpt3
	  p3(3) = 0.d0
	  p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	  p5(1) = 0.d0
	  p5(2) = 0.d0
	  p5(3) = 0.d0
	  p5(4) = 0.d0
c ------------------------------------------- 
	  p4(1) = gpt4*dcosh(y4)
	  p4(2) = -gpt4
	  p4(3) = 0.d0
	  p4(4) = gpt4*dsinh(y4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimo2 = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
** 	scale_fixe = 100.d0
** 	m = cm*scale_fixe
** 	mu = cmu*scale_fixe
** 	mf = cmf*scale_fixe
	alpha_em = alphaem(iloopem,m*m)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s*alpha/(4*c_i*c_j*s**2)
	cbij = alphas*alpha_em*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	call ssy3po(s,gpt3,gpt4,y3,y4,mf,x10,x20,x3,z3,ih1,ih2,
     #	ih3,sy3p)
	call vec_dinit(sy3p,k0max,zero)	  
	if (z3.le.z3sup) then
	  call ssy3ro(s,gpt4,y3,y4,x10,x20,x3,z3,ih1,ih2,ih3,sy3r)
	else
	  call vec_dinit(sy3r,k0max,zero)	  
	endif
	call vec_dadd_vector(sy3p,sy3r,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dim = 0.d0
	do i = j_processo_min,j_processo_max
	  f2dim = temp1(i) + temp1(i+18) + f2dim
	enddo
cc	f2dim = vec_dsum(temp0,k0max)
c
	f2dimo2 = f2dim * xjacob * alpi * cbij * hc2
	if (ispring.eq.1) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f2dimo2
	  f2dimo2 = dabs(f2dimo2)
	endif
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimo3(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/alem/iloopem
	common/faconv/hc2
	common/coup/ptmm,r
	common/flagexp/icoup_exp
	common/coupexp/ptm_exp,r_exp
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension sy4r(k0max),sy4p(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	zero = 0.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3))
     #	,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ gpt4 ------------------------------------------
	if (i_flag_iso.eq.1) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3*(1.d0+etmax))
	endif
	gpt4max = dmax1(gpt3,gpt4maxp)
	gpt4min = gpt3
	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(4)
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	x10 = gpt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = gpt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimo3 = 0.d0
	  return
	endif
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
c z4 doit etre > a 2*gpt4/dsqrt(s)*dcosh(y3) et 2*gpt4/dsqrt(s)*dcosh(y4)
	x4min = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	x3min = 2.d0*gpt4/dsqrt(s)*dcosh(y3)
	if (i_flag_iso.eq.1) then
	  cut_isol1 = gpt4/(gpt3+etmax)
	  cut_isol2 = gpt4/(gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  cut_isol1 = gpt4/(gpt3*(1.d0+etmax))
	  cut_isol2 = 1.d0/(1.d0+etmax)
	endif
	z4max0 = 1.d0
	z4min0 = dmax1(x4min,x3min,x10,x20,cut_isol1,cut_isol2)
	z4min1 = gpt4/(ptm_exp+gpt4)	
	z4sup = gpt4/(ptm+gpt4)
	if (icoup_exp.eq.0) then
	  z4min = z4min0	
	  z4max = z4max0	
	else if (icoup_exp.eq.1) then
	  z4min = dmax1(z4min0,z4min1)	
	  z4max = z4sup	
	endif
c en principe vrai mais a verif
	if (z4min.gt.z4max) then
	  f2dimo3 = 0.d0
	  write (8,*) 'attention z4 ne verifie pas z4min < z4 < z4max',
     #		z4,z4min,z4max
	  return
	endif
	z4 = z4min + (z4max-z4min)*xx(5)
	x3 = gpt3/gpt4*z4
c en principe vrai mais a verif
	x3minp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	if (x3.gt.un.or.x3.lt.x3minp) then
	  f2dimo3 = 0.d0
	  write (8,*) 'attention x3 ne verifie pas x3min < x3 < x3max',
     #		x3,x3minp,un
	  return
	endif
c
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)
     #	*(gpt4max-gpt4min)*(z4max-z4min)
	call fspcouo(n,spcou)
c ------------------------------------------- 
	  p3(1) = gpt3*dcosh(y3)
	  p3(2) = gpt3
	  p3(3) = 0.d0
	  p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	  p5(1) = 0.d0
	  p5(2) = 0.d0
	  p5(3) = 0.d0
	  p5(4) = 0.d0
c ------------------------------------------- 
	  p4(1) = gpt4*dcosh(y4)
	  p4(2) = -gpt4
	  p4(3) = 0.d0
	  p4(4) = gpt4*dsinh(y4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimo3 = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alpha_em = alphaem(iloopem,m*m)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s*alpha/(4*c_i*c_j*s**2)
	cbij = alphas*alpha_em*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	call ssy4po(s,gpt3,y3,y4,mf,gpt4,x10,x20,x3,z4,ih1,ih2,ih3,sy4p)
	if (z4.le.z4sup) then
	  call ssy4ro(s,gpt3,y3,y4,gpt4,x10,x20,x3,z4,ih1,ih2,ih3,sy4r)
	else
	  call vec_dinit(sy4r,k0max,zero)	  
	endif
	call vec_dadd_vector(sy4p,sy4r,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dim = 0.d0
	do i = j_processo_min,j_processo_max
	  f2dim = temp1(i) + temp1(i+18) + f2dim
	enddo
cc	f2dim = vec_dsum(temp0,k0max)
c
	f2dimo3 = f2dim * xjacob * alpi * cbij * hc2
	if (ispring.eq.1) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f2dimo3
	  f2dimo3 = dabs(f2dimo3)
	endif
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimo3p(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/alem/iloopem
	common/faconv/hc2
	common/coup/ptmm,r
	common/flagexp/icoup_exp
	common/coupexp/ptm_exp,r_exp
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension sy4r(k0max),sy4p(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	zero = 0.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
** c ------------------ gpt3 ------------------------------------------
** 	rs_over_plus = dsqrt(s)/(dexp(y3)+dexp(y4))
** 	rs_over_moins = dsqrt(s)/(dexp(-y3)+dexp(-y4))
** 	gpt3maxpp = dmin1(rs_over_plus,rs_over_moins)
** 	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3))
**      #	,dsqrt(s)/(2.d0*dcosh(y4)),gpt3maxpp)
** 	gpt3max = dmax1(gptinf,gpt3maxp)
** 	gpt3min = gptinf
** 	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
** c ------------------ gpt4 ------------------------------------------
** 	gpt4maxp = gpt3
** 	gpt4max = dmax1(gptinf,gpt4maxp)
** 	gpt4min = gptinf
** 	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(4)
c ------------------ gpt4 ------------------------------------------
	rs_over_plus = dsqrt(s)/(dexp(y3)+dexp(y4))
	rs_over_moins = dsqrt(s)/(dexp(-y3)+dexp(-y4))
	gpt4maxpp = dmin1(rs_over_plus,rs_over_moins)
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3))
     #	,dsqrt(s)/(2.d0*dcosh(y4)),gpt4maxpp)
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(3)
c ------------------ gpt3 ------------------------------------------
	rs_over_plus = dsqrt(s)/(dexp(y3)+dexp(y4))
	rs_over_moins = dsqrt(s)/(dexp(-y3)+dexp(-y4))
	gpt3maxpp = dmin1(rs_over_plus,rs_over_moins)
	if (i_flag_iso.eq.1) then
	  gpt3maxp = dmin1(gptsup
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3maxpp,gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  gpt3maxp = dmin1(gptsup
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3maxpp,gpt4*(1.d0+etmax))
	endif
	gpt3max = dmax1(gpt4,gpt3maxp)
	gpt3min = gpt4
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(4)
	x10 = gpt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = gpt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimo3p = 0.d0
	  return
	endif
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
c z4 doit etre > a 2*gpt4/dsqrt(s)*dcosh(y3) et 2*gpt4/dsqrt(s)*dcosh(y4)
	x4min = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	x3min = 2.d0*gpt4/dsqrt(s)*dcosh(y3)
	if (i_flag_iso.eq.1) then
	  cut_isol1 = gpt4/(gpt3+etmax)
	  cut_isol2 = gpt4/(gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  cut_isol1 = gpt4/(gpt3*(1.d0+etmax))
	  cut_isol2 = 1.d0/(1.d0+etmax)
	endif
	z4max0 = gpt4/gpt3
	z4min0 = dmax1(x4min,x3min,x10,x20,cut_isol1,cut_isol2)
	z4min1 = gpt4/(ptm_exp+gpt4)	
	z4sup = gpt4/(ptm+gpt4)
	un = 1.d0
	if (icoup_exp.eq.0) then
	  z4min = z4min0	
	  z4max = z4max0	
	else if (icoup_exp.eq.1) then
	  z4min = dmax1(z4min0,z4min1)	
	  z4max = z4sup	
	endif
c en principe vrai mais a verif
	if (z4min.gt.z4max) then
	  f2dimo3p = 0.d0
	  write (18,*) 'attention z4 ne verifie pas z4min < z4 < z4max',
     #		z4,z4min,z4max
	  return
	endif
	z4 = z4min + (z4max-z4min)*xx(5)
	x3 = gpt3/gpt4*z4
c en principe vrai mais a verif
	x3minp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	if (x3.gt.un.or.x3.lt.x3minp) then
	  f2dimo3p = 0.d0
	  write (18,*) 'attention x3 ne verifie pas x3min < x3 < x3max',
     #		x3,x3minp,un
	  return
	endif
c
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)
     #	*(gpt4max-gpt4min)*(z4max-z4min)
	call fspcouo(n,spcou)
c ------------------------------------------- 
	  p3(1) = gpt3*dcosh(y3)
	  p3(2) = gpt3
	  p3(3) = 0.d0
	  p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	  p5(1) = 0.d0
	  p5(2) = 0.d0
	  p5(3) = 0.d0
	  p5(4) = 0.d0
c ------------------------------------------- 
	  p4(1) = gpt4*dcosh(y4)
	  p4(2) = -gpt4
	  p4(3) = 0.d0
	  p4(4) = gpt4*dsinh(y4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimo3p = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alpha_em = alphaem(iloopem,m*m)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s*alpha/(4*c_i*c_j*s**2)
	cbij = alphas*alpha_em*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	call ssy4pop(s,gpt3,y3,y4,mf,gpt4,x10,x20,x3,z4,ih1,ih2,
     #	ih3,sy4p)
	if (z4.le.z4sup) then
	  call ssy4ro(s,gpt3,y3,y4,gpt4,x10,x20,x3,z4,ih1,ih2,ih3,sy4r)
	else
	  call vec_dinit(sy4r,k0max,zero)	  
	endif
	call vec_dadd_vector(sy4p,sy4r,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dim = 0.d0
	do i = j_processo_min,j_processo_max
	  f2dim = temp1(i) + temp1(i+18) + f2dim
	enddo
cc	f2dim = vec_dsum(temp0,k0max)
c
	f2dimo3p = f2dim * xjacob * alpi * cbij * hc2
	if (ispring.eq.1) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f2dimo3p
	  f2dimo3p = dabs(f2dimo3p)
	endif
	return
	end 
c*****************************************************
c integrale a 1 dimension
	double precision function f1dimo(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/alem/iloopem
	common/faconv/hc2
	common/born/iborn
cc	common/sortie/p3(4),p4(4),p5(4),resfunc
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/aurenche/iauren
	common/processo/j_processo_min,j_processo_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=36)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	dimension c0sv12(k0max),c0sv3(k0max),c0sv4(k0max),c0fborn(k0max)
	dimension spcou(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3))
     #	,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ gpt4 ------------------------------------------
c x3 = gpt3/gpt4 et x3min=2*gpt3/dsqrt(s)*dcosh(y3) <= gpt3/gpt4
c d'ou gpt4 <= dsqrt(s)/(2.d0*dcosh(y3))
c pour les coupures d'isolement on a que
c gpt3 <= gpt4 <= gpt3+etmax car
c x3=gpt3/gpt4 >= gpt3/(gpt3+etmax)
	if (i_flag_iso.eq.1) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	  ,dsqrt(s)/(2.d0*dcosh(y3)),gpt3*(1.d0+etmax))
	endif
	gpt4max = dmax1(gpt3,gpt4maxp)
	gpt4min = gpt3
	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(4)
c
	x10 = gpt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = gpt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f1dimo = 0.d0
	  return
	endif
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)
     #	*(gpt4max-gpt4min)
	call fspcouo(n,spcou)
	x3 = gpt3/gpt4
c en principe vrai mais a verif
	x3max = 1.d0
	x3minp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3min = dmin1(x3max,x3minp)
	if (x3.gt.x3max.or.x3.lt.x3min) then
	  f1dimo = 0.d0
	  write (18,*) 'attention x3 ne verifie pas x3min < x3 < x3max',
     #		x3,x3min,x3max
	endif
c ------------------------------------------- 
	  p3(1) = gpt3*dcosh(y3)
	  p3(2) = gpt3
	  p3(3) = 0.d0
	  p3(4) = gpt3*dsinh(y3)
c ------------------------------------------- 
	  p5(1) = 0.d0
	  p5(2) = 0.d0
	  p5(3) = 0.d0
	  p5(4) = 0.d0
c ------------------------------------------- 
	  p4(1) = gpt4*dcosh(y4)
	  p4(2) = -gpt4
	  p4(3) = 0.d0
	  p4(4) = gpt4*dsinh(y4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f1dimo = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	if (iauren.eq.0) then
	  call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	else if (iauren.eq.1) then
	  m = gpt3
	  mu = m
	  mf = m
	endif
** 	scale_fixe = 100.d0
** 	m = cm*scale_fixe
** 	mu = cmu*scale_fixe
** 	mf = cmf*scale_fixe
	alpha_em = alphaem(iloopem,m*m)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s*alpha/(4*c_i*c_j*s**2)
	cbij = alphas*alpha_em*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c 
	if (iborn.eq.0) then 
	  call ssv12o(s,gpt3,y3,y4,m,mu,x10,x20,x3,ih1,ih2,ih3,c0sv12)
	  call ssv3o(s,gpt3,gpt4,y3,y4,mf,x10,x20,x3,ih1,ih2,
     #	  ih3,c0sv3)
	  call ssv4o(s,gpt3,gpt4,etmax,y3,y4,mf,x10,x20,x3,ih1,
     #	  ih2,ih3,c0sv4)
	  call vec_dmult_vector(c0sv12,spcou,k0max,temp0)
	  call vec_dmult_vector(c0sv3,spcou,k0max,temp1)
	  call vec_dmult_vector(c0sv4,spcou,k0max,temp2)
	  sv12 = 0.d0
	  sv3 = 0.d0
	  sv4 = 0.d0
	  do i = j_processo_min,j_processo_max
	    sv12 = temp0(i) + temp0(i+18) + sv12
	    sv3 = temp1(i) + temp1(i+18) + sv3
	    sv4 = temp2(i) + temp2(i+18) + sv4
	  enddo
cc	  sv12 = vec_dsum(temp0,k0max)
cc	  sv3 = vec_dsum(temp1,k0max)
cc	  sv4 = vec_dsum(temp2,k0max)
	  fborn = 0.d0
	else if (iborn.eq.1) then
	  call sfborno(s,n,gpt3,y3,y4,x3,ih1,ih2,ih3,c0fborn)
	  call vec_dmult_vector(c0fborn,spcou,k0max,temp3)
	  fborn = 0.d0
	  do i = j_processo_min,j_processo_max
	    fborn = temp3(i) + temp3(i+18) +fborn
	  enddo
cc	  fborn = vec_dsum(temp3,k0max)
	  sv12 = 0.d0
	  sv3 = 0.d0
	  sv4 = 0.d0
	endif
c 
	resb = fborn + (sv12+sv3+sv4)*alpi
	f1dimo = resb * xjacob * cbij * hc2
c
	if (ispring.eq.1) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  xfrag1 = x3
	  xfrag2 = un
	  resfunc = f1dimo
	  f1dimo = dabs(f1dimo)
	endif
c
	return
	end 
c===================================================================
	subroutine ampo3(s,x1,x2,y3,pt3,y5,pt5,fi35,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=36)
	dimension h12(k0max),h13(k0max),h14(k0max)
	dimension h34(k0max),h23(k0max),h24(k0max),cons(k0max)
	dimension camp(k0max)
	s12 =  x1*x2*s/2.d0
	s13 =  x1*pt3*dsqrt(s)/2.d0*dexp(-y3)
	s23 =  x2*pt3*dsqrt(s)/2.d0*dexp(y3)
	s15 =  x1*pt5*dsqrt(s)/2.d0*dexp(-y5)
	s25 =  x2*pt5*dsqrt(s)/2.d0*dexp(y5)
	s35 =  pt3*pt5*coshmcos(y3-y5,fi35)
	s14 =  s12-s13-s15
	s24 =  s12-s23-s25
	s34 =  s13+s23-s35
	s45 =  s15+s25-s35
	e12 = s12/(s15*s25)
	e13 = s13/(s15*s35)
	e23 = s23/(s25*s35)
	e34p = s34/(s15+s25)/s35
c -------------------------------------------
	call sh12o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h12)
	call sh13o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h13)
	call sh23o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h23)
	call sh34o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	call sconso(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,cons)
	do i=1,k0max
	  camp(i) = (0.5d0*h12(i)*e12+h13(i)*e13
     #	             +h23(i)*e23+h34(i)*e34p+0.5d0*cons(i))
	enddo
	return
	end
c===================================================================
	subroutine ampo_corr3(s,x1,x2,y3,pt3,pt5,rs,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=36)
	dimension h12(k0max),h13(k0max),h14(k0max)
	dimension h34(k0max),h23(k0max),h24(k0max),cons(k0max)
	dimension camp(k0max)
	s12 =  x1*x2*s/2.d0
	s13 =  x1*pt3*dsqrt(s)/2.d0*dexp(-y3)
	s23 =  x2*pt3*dsqrt(s)/2.d0*dexp(y3)
	s15 =  x1*pt5*dsqrt(s)/2.d0*dexp(-y3)
	s25 =  x2*pt5*dsqrt(s)/2.d0*dexp(y3)
	s35 =  0.d0
	s14 =  s12-s13-s15
	s24 =  s12-s23-s25
	s34 =  s13+s23-s35
	s45 =  s15+s25-s35
c -------------------------------------------
	call sh13o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h13)
	call sh23o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h23)
	call sh34o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	do i=1,k0max
	  camp(i) = (h13(i)+h23(i)+h34(i))*2.d0/rs/(pt5*pt5)
	enddo
	return
	end
c===================================================================
	subroutine ampo4(s,x1,x2,y4,pt4,y5,pt5,fi45,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=36)
	dimension h12(k0max),h13(k0max),h14(k0max)
	dimension h34(k0max),h23(k0max),h24(k0max),cons(k0max)
	dimension camp(k0max)
	s12 =  x1*x2*s/2.d0
	s14 =  x1*pt4*dsqrt(s)/2.d0*dexp(-y4)
	s24 =  x2*pt4*dsqrt(s)/2.d0*dexp(y4)
	s15 =  x1*pt5*dsqrt(s)/2.d0*dexp(-y5)
	s25 =  x2*pt5*dsqrt(s)/2.d0*dexp(y5)
	s45 =  pt4*pt5*coshmcos(y4-y5,fi45)
	s13 =  s12-s14-s15
	s23 =  s12-s24-s25
	s34 =  s14+s24-s45
	s35 =  s15+s25-s45
	e12 = s12/(s15*s25)
	e13 = s13/(s15*s35)
	e23 = s23/(s25*s35)
	e14 = s14/(s15*s45)
	e24 = s24/(s25*s45)
	e34s = s34/(s15+s25)/s45
c -------------------------------------------
	call sh12o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h12)
	call sh14o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h14)
	call sh24o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h24)
	call sh34o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	call sconso(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,cons)
	do i=1,k0max
	  camp(i) = (0.5d0*h12(i)*e12+h14(i)*e14
     #	             +h24(i)*e24+h34(i)*e34s+0.5d0*cons(i))
	enddo
	return
	end
c===================================================================
	subroutine ampo_corr4(s,x1,x2,y4,pt4,pt5,rs,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=36)
	dimension h12(k0max),h13(k0max),h14(k0max)
	dimension h34(k0max),h23(k0max),h24(k0max),cons(k0max)
	dimension camp(k0max)
	s12 =  x1*x2*s/2.d0
	s14 =  x1*pt4*dsqrt(s)/2.d0*dexp(-y4)
	s24 =  x2*pt4*dsqrt(s)/2.d0*dexp(y4)
	s15 =  x1*pt5*dsqrt(s)/2.d0*dexp(-y4)
	s25 =  x2*pt5*dsqrt(s)/2.d0*dexp(y4)
	s45 =  0.d0
	s13 =  s12-s14-s15
	s23 =  s12-s24-s25
	s34 =  s14+s24-s45
	s35 =  s15+s25-s45
c -------------------------------------------
	call sh14o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h14)
	call sh24o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h24)
	call sh34o(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	do i=1,k0max
	  camp(i) = (h14(i)+h24(i)+h34(i))*2.d0/rs/(pt5*pt5)
	enddo
	return
	end
c********************************************************************
	subroutine count(i)
	implicit real*8(a-h,l-z)
	common/testcount/j(11)
	j(i) = j(i) + 1
	return
	end
c********************************************************************
	subroutine trapnan(i,wx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	wzero = 0.
	if (wx.lt.wzero.and.wx.gt.wzero) then
	  write (*,*) 'i=',' nan=',wx
	endif
	return
	end
c********************************************************************

