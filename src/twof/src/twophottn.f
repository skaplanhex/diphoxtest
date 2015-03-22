c
c version corrige le 6/01/97
c modification le 21/01/97 pour generer le fi_3 entre 0 et 2 pi
c modification le 22/01/97 pour la separation entre les deux photons
c modification le 14/02/97 pour mettre la coupure sur la masse
c invariante dans les quasi deux donne deux et les deux donne deux
c modification le 7/03/97 pour fixer le nombre total d'evenements au
c depart
c modification le 7/04/97 pour fixer un bug dans f2dim2 et f2dim3
c (modifications aussi dans twofrag.f: subroutine ssv3 et ssv4)
c modification le 9/06/97 pour l'integration dans f4dim3 et f4dim4
c modification le 3/05/99 pour ajout des termes en r^2*ln(p_{tm})
c modification le 14/05/99 pour introduire le maximum de la masse 
c invariante photon-photon
c correction d'un bug dans la generation le 11/08/99 dans f2dim2 et f2dim3
c bug dans la generation des 2 -> 3 (4//5) (15/05/00)
c bug dans ssy4pt et ssy4ptp
	subroutine twof_sub
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical lo,nlo,gener,integ,intega
	logical verif
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
**         common/sortier/ntrack,we(maxtrk),wpx(maxtrk),
**      #	wpy(maxtrk),wpz(maxtrk)
** 	common/sfragvr/wx3,wx4
	common/randd/rnumb
	common/cheminbs/path_bsfile
	common/long/ilen
	common/efficacite/ipositif,inegatif
	dimension wrvec(irandom),wrvec1(irandom)
	external f4dimt3,f4dimt4,f2dimt1,f2dimt2,f2dimt2p,f2dimt3,
     #	f2dimt3p,f1dimt
c
	pi=datan(1.d0)*4.d0
	ipositif = 0
	inegatif = 0
c******************************************************************
c       partie 2 --> 3
c******************************************************************
cccc
	if (nlo) then
cccc
	iborn = 0
c###
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 8
	iwild = 8
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimt3,resf34,sdf34,ctime,it1,it2)
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
	idim = 8
	iwild = 8
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimt3,resf34,sdf34,ctime,it1,it2)
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
	idim = 8
	iwild = 8
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimt4,resf44,sdf44,ctime,it1,it2)
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
	idim = 8
	iwild = 8
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimt4,resf44,sdf44,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'threepart2.bs',
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt1,resf21,sdf21,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt1,resf21,sdf21,ctime,it1,it2)
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
	if (integ) then
	ispring = 0
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt2,resf22,sdf22,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt2,resf22,sdf22,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'qtwopart2.bs',
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt2p,resf22p,sdf22p,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt2p,resf22p,sdf22p,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'qtwopart2p.bs',
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt3,resf23,sdf23,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt3,resf23,sdf23,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt3p,resf23p,sdf23p,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimt3p,resf23p,sdf23p,ctime,it1,it2)
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
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	call bases (f1dimt,resf1,sdf1,ctime,it1,it2)
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
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	call bases (f1dimt,resf1,sdf1,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'twopart1.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
	endif
c###
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
	idim = 5
	iwild = 5
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	call bases (f1dimt,resfb,sdfb,ctime,it1,it2)
c on ecrit la grille
	open(28,file=path_bsfile(1:ilen)//'twopart2.bs',
     #	status='unknown',form='unformatted')
	call bswrit(28)
	close(28)
	endif
cccc
	endif
cccc
c******************************************************************
c       calcul de la section efficace totale ordre superieur
c******************************************************************
	if (integ) then
	  resf2 = resf21+resf22+resf23+resf22p+resf23p
	  sdf2 = sdf21+sdf22+sdf23+sdf22p+sdf23p
c
	  resho = resf34+resf44+resf2+resf1
	  sdho = sdf34+sdf44+sdf2+sdf1
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
	    resf2 = resf21+resf22+resf23+resf22p+resf23p
	    totho = resf34+resf44+resf2+resf1
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
	  write(28,*) resf34, resf44
	  write(28,*) resf21, resf22, resf22p, resf23, resf23p
	  write(28,*) resf1, resfb
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
	read(28,*) resf34, resf44
	read(28,*) resf21, resf22, resf22p, resf23, resf23p
	read(28,*) resf1, resfb
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
	idim = 8
	iwild = 8
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
          endif
	  rnumb = dble(wrvec(i1+1))
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f4dimt3,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c
c###
	ispring = 1
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 8
	iwild = 8
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
          endif
	  rnumb = dble(wrvec(i1+1))
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f4dimt4,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
	idim = 6
	iwild = 6
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
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f2dimt1,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
	idim = 6
	iwild = 6
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
	write(6,*)'first part'
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
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f2dimt2,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
	idim = 6
	iwild = 6
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	inbrealevent = int(resf22p*xnorm)
c
	open(28,file=path_bsfile(1:ilen)//'qtwopart2p.bs',
     #	status='unknown',form='unformatted')
	call bsread(28)
	close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : quasi 2 partons : col. fi3'
	write(6,*)'second part'
	write(6,*) inbrealevent,' events'
	write(6,*)'============================================='
	write(6,*)
c	
        iprov = 32
        ntrack = 2
        poid = 1.d0
	do i = 1,inbrealevent
	  i1 = mod(i-1,irandom)
	  if (i1.eq.0) then
            call ranlux(wrvec1,irandom)
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f2dimt2p,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
	idim = 6
	iwild = 6
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
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f2dimt3,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
	idim = 6
	iwild = 6
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
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f2dimt3p,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
	idim = 5
	iwild = 5
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
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f1dimt,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
	idim = 5
	iwild = 5
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
          endif
	  rnumb1 = dble(wrvec1(i1+1))
	  call spring(f1dimt,ixtry)
	  call doubletosinglet(rnumb1,pi,poid,iprov,ntrack)
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
c
	return
	end
c
c*****************************************************
c  partie ou pt5 < ptm
c*****************************************************
c integrale a 4 dimensions
	double precision function f4dimt3(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical isol3,isol3p,isol4,isol4p
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/coup/ptm,r
	common/flag/icoup
	common/coupexp/ptm_exp,r_exp
	common/flagexp/icoup_exp
	common/scale/m,mf,mu
	common/faconv/hc2
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/randd/rnumb
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension sf0(k0max),camp0(k0max)
	dimension temp0(k0max),temp1(k0max),spcou(k0max)
	dimension p1(4),p2(4),pp3(4),pp4(4)
	dimension p3(4),p4(4),p5(4)
	dimension p1p(4),p2p(4),p4p(4),p5p(4),pp4p(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
	expos = 7.d0
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
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(3)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(4)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
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
	x3 = x3min + (x3max-x3min)*xx(5)
	pt3 = gpt3/x3
c ------------------ pt5 --------------------------------------------
cc	pt5max = dsqrt(pt3**2+gptsup**2)
	pt5max = pt3+dsqrt(s)/(2.d0*dcosh(y4))
	pt5min = ptm
	pt5 = pt5min + (pt5max-pt5min)*xx(6)
c ------------------ fi35 -------------------------------------------
	x11 = pt3/dsqrt(s)*dexp(y3)+dabs(pt3-pt5)/dsqrt(s)*dexp(y4)
	x21 = pt3/dsqrt(s)*dexp(-y3)+dabs(pt3-pt5)/dsqrt(s)*dexp(-y4)
	if (x11.gt.un.or.x21.gt.un) then
	  f4dimt3 = 0.d0
	  return
	endif
	sy5max = dlog(dsqrt(s)/pt5*(1.d0-x11))
	sy5min = -dlog(dsqrt(s)/pt5*(1.d0-x21))
	if (sy5max.lt.sy5min) then
	  f4dimt3 = 0.d0
	  return
	endif
	if (y3.le.sy5max.and.y3.ge.sy5min) then
	  omegamin = 0.d0
	  omegamax = pi
	  srmax2 = dmax1((sy5max-y3)**2+pi**2,(sy5min-y3)**2+pi**2)
	  rmax2 = dsqrt(srmax2)
	else if (y3.lt.sy5max.and.y3.lt.sy5min) then
	  omegamin = 0.d0
	  omegamax = pi/2.d0
	  srmax2 = (sy5max-y3)**2+pi**2
	  rmax2 = dsqrt(srmax2)
	else if (y3.gt.sy5max.and.y3.gt.sy5min) then
	  omegamin = pi/2.d0
	  omegamax = pi
	  srmax2 = (sy5min-y3)**2+pi**2
	  rmax2 = dsqrt(srmax2)
	endif
	omega = omegamin + (omegamax-omegamin)*xx(7)
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
	rr = rrmin + (rrmax-rrmin)*xx(8)
	fi35 = rr*dsin(omega)
	if (fi35.gt.pi) then
	   f4dimt3 = 0.d0
	   return
	endif
c ------------------ pt4 -------------------------------------------
	pt4 = dsqrt(pt3**2+pt5**2+2.d0*pt3*pt5*dcos(fi35))
	x4 = gpt4/pt4
	x4minp = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	x4min = dmin1(un,x4minp)
	if (x4.lt.x4min.or.x4.gt.un) then
	  f4dimt3 = 0.d0
	  return
	endif
c ------------------ x -------------------------------------------
	x1c = pt3/dsqrt(s)*dexp(y3) + pt4/dsqrt(s)*dexp(y4)
	x2c = pt3/dsqrt(s)*dexp(-y3) + pt4/dsqrt(s)*dexp(-y4)
	if (x1c.ge.un.or.x2c.ge.un) then
	  f4dimt3 = 0.d0
	  return
	endif
c ------------------ y5 -------------------------------------------
	y5max = dlog(dsqrt(s)/pt5*(1.d0-x1c))
	y5min = -dlog(dsqrt(s)/pt5*(1.d0-x2c))
	y5 = y3 + rr*dcos(omega)
	if (y5.gt.y5max.or.y5.lt.y5min) then
	   f4dimt3 = 0.d0
	   return
	endif
c ------------------ x -------------------------------------------
	x1 = x1c + pt5/dsqrt(s)*dexp(y5)
	x2 = x2c + pt5/dsqrt(s)*dexp(-y5)
	if (x1.ge.un.or.x2.ge.un) then
	  f4dimt3 = 0.d0
	  return
	endif
c coupure en angle-------------------------------------------------
        r35s = (y3-y5)**2+fi35**2
        rs = r**2
	icorr = 0
        if (r35s.lt.rs) then
           f4dimt3 = 0.d0
           return
**           icorr = 1
        endif
c coupure en isolement-------------------------------------------------
	call isolement(y3,y5,fi35,x3,pt5,gpt3,r_isol,etmax,
     #	isol3)
	cos_fi45 = -(pt3*dcos(fi35)+pt5)/pt4
	fi45 = zacos(cos_fi45)
	call isolement(y4,y5,fi45,x4,pt5,gpt4,r_isol,etmax,
     #	isol4)
** 	if (.not.(isol3)) then
** 	  f4dimt3 = 0.d0
** 	  return
** 	endif
	xisol4 = 1.d0
	if (.not.(isol4).or..not.(isol3)) then
	  if (icorr.eq.0) then
	    f4dimt3 = 0.d0
	    return
	  else if (icorr.eq.1) then
	    xisol4 = 0.d0
	  endif
	endif
c ------------------ fin isolement ------------------------------------
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu3max-uu3min)*
     #	rr*(rrmax-rrmin)*(omegamax-omegamin)*(x3max-x3min)*
     #	(pt5max-pt5min)*(uu4max-uu4min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
	pt4p = (pt3+pt5)
	x1p = pt4p/dsqrt(s)*(dexp(y3)+dexp(y4))
	x2p = pt4p/dsqrt(s)*(dexp(-y3)+dexp(-y4))
	x4p = gpt4/pt4p
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
	pp4(1) = p1(1) + p2(1) - pp3(1) - p5(1)
	pp4(2) = p1(2) + p2(2) - pp3(2) - p5(2)
	pp4(3) = p1(3) + p2(3) - pp3(3) - p5(3)
	pp4(4) = p1(4) + p2(4) - pp3(4) - p5(4)
c ------------------------------------------- 
	p4(1) = x4*pp4(1)
	p4(2) = x4*pp4(2)
	p4(3) = x4*pp4(3)
	p4(4) = x4*pp4(4)
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
	pp4p(1) = p1p(1) + p2p(1) - pp3(1) - p5p(1)
	pp4p(2) = p1p(2) + p2p(2) - pp3(2) - p5p(2)
	pp4p(3) = p1p(3) + p2p(3) - pp3(3) - p5p(3)
	pp4p(4) = p1p(4) + p2p(4) - pp3(4) - p5p(4)
c ------------------------------------------- 
	p4p(1) = x4p*pp4p(1)
	p4p(2) = x4p*pp4p(2)
	p4p(3) = x4p*pp4p(3)
	p4p(4) = x4p*pp4p(4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	mass_cut = 1.d0
	if (iphocut.eq.0) then
c angle azymuthal entre les deux photons plus grand que fimin
	  abspt3 = dsqrt(p3(2)**2+p3(3)**2)
	  abspt4 = dsqrt(p4(2)**2+p4(3)**2)
	  fi34 = zacos((p3(2)*p4(2)+p3(3)*p4(3))/(abspt3*abspt4))
	  if (fi34.le.fimin) then
	    f4dimt3 = 0.d0
	    return
	  endif
	else if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    if (icorr.eq.0) then
	      f4dimt3 = 0.d0
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
	    f4dimt3 = 0.d0
	    return
	  endif
	endif	
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	call fspcout(n,spcou)
c
c c_{ij} = alpha_s^3/(4*c_i*c_j*pi*s^2)
	cij = alfas(iloop,mu*mu)**3/(pi*s*s)
c
	call strfrat(x1,ih1,x2,ih2,x3,ih3,x4,ih4,sf0)
cc	call ampt(3,p1,p2,pp3,pp4,p5,camp0)
	call ampt3(s,x1,x2,y3,pt3,y5,pt5,fi35,camp0)
	call vec_dmult_vector(camp0,sf0,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f4dimt3 = 0.d0
	do i = j_processt_min,j_processt_max
	  f4dimt3 = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	  temp1(i+63) + f4dimt3
	enddo
	f4dimt3 = cij*f4dimt3/pt4
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
	  call isolement(y4,y3,pi,x4p,pt5,gpt4,r_isol,etmax,
     #	  isol4p)
	  xisol4p = 1.d0
	  if (.not.(isol4p).or..not.(isol3p)) then
	    xisol4p = 0.d0
	  endif
c ------------------ fin isolement ------------------------------------
	  call choiscale(p3,p4p,cm,cmu,cmf,ichoi_scale,m,mu,mf)
c c_{ij} = alpha_s^3/(4*c_i*c_j*pi*s^2)
	  cij = alfas(iloop,mu*mu)**3/(pi*s*s)
c
	  call strfrat(x1p,ih1,x2p,ih2,x3,ih3,x4p,ih4,sf0)
	  call ampt_corr3(s,x1p,x2p,y3,pt3,pt5,r35s,camp0)
	  call vec_dmult_vector(camp0,sf0,k0max,temp0)
	  call vec_dmult_vector(temp0,spcou,k0max,temp1)
	  f4dimt3p = 0.d0
	  do i = j_processt_min,j_processt_max
	    f4dimt3p = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	    temp1(i+63) + f4dimt3p
	  enddo
	  f4dimt3p = cij*f4dimt3p/pt4p
	  f4dimt3 = f4dimt3*mass_cut*xisol4 - 
     #	  f4dimt3p*mass_cutp*xisol4p
	endif
c
	f4dimt3 = f4dimt3 * xjacob * (pt3/x3*pt5) * hc2
	if (facteur_symetrie) then
	  f4dimt3 = f4dimt3/2.d0
	endif	
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
	  xfrag2 = x4
	  resfunc = f4dimt3
	  f4dimt3 = dabs(f4dimt3)
	endif
c
	return
	end 
c********************************************************************
c integrale a 4 dimensions
	double precision function f4dimt4(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical isol3,isol3p,isol4,isol4p
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/coup/ptm,r
	common/flag/icoup
	common/coupexp/ptm_exp,r_exp
	common/flagexp/icoup_exp
	common/scale/m,mf,mu
	common/faconv/hc2
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/randd/rnumb
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension sf0(k0max),camp0(k0max)
	dimension temp0(k0max),temp1(k0max),spcou(k0max)
	dimension p1(4),p2(4),pp3(4),pp4(4)
	dimension p3(4),p4(4),p5(4)
	dimension p1p(4),p2p(4),p3p(4),pp3p(4),p5p(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
	expos = 7.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(3)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(4)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
c ------------------ x4 -------------------------------------------
	x4max = 1.d0
	x4min_incl = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	if (i_flag_iso.eq.1) then
	  x4min_isol = gpt4/(gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  x4min_isol = 1.d0/(1.d0+etmax)
	endif
	x4minp = dmax1(x4min_incl,x4min_isol)
	x4min = dmin1(x4max,x4minp)
	x4 = x4min + (x4max-x4min)*xx(5)
	pt4 = gpt4/x4
c ------------------ pt5 --------------------------------------------
cc	pt5max = dsqrt(pt4**2+gptsup**2)
	pt5max = pt4+dsqrt(s)/(2.d0*dcosh(y3))
	pt5min = ptm
	pt5 = pt5min + (pt5max-pt5min)*xx(6)
c ------------------ fi45 -------------------------------------------
	x11 = pt4/dsqrt(s)*dexp(y4)+dabs(pt4-pt5)/dsqrt(s)*dexp(y3)
	x21 = pt4/dsqrt(s)*dexp(-y4)+dabs(pt4-pt5)/dsqrt(s)*dexp(-y3)
	if (x11.gt.un.or.x21.gt.un) then
	  f4dimt4 = 0.d0
	  return
	endif
	sy5max = dlog(dsqrt(s)/pt5*(1.d0-x11))
	sy5min = -dlog(dsqrt(s)/pt5*(1.d0-x21))
	if (sy5max.lt.sy5min) then
	  f4dimt4 = 0.d0
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
	omega = omegamin + (omegamax-omegamin)*xx(7)
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
	rr = rrmin + (rrmax-rrmin)*xx(8)
	fi45 = rr*dsin(omega)
	if (fi45.gt.pi) then
	   f4dimt4 = 0.d0
	   return
	endif
c ------------------ pt3 -------------------------------------------
	pt3 = dsqrt(pt4**2+pt5**2+2.d0*pt5*pt4*dcos(fi45))
	x3 = gpt3/pt3
	x3minp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3min = dmin1(un,x3minp)
	if (x3.lt.x3min.or.x3.gt.un) then
	  f4dimt4 = 0.d0
	  return
	endif
c ------------------ x -------------------------------------------
	x1c = pt3/dsqrt(s)*dexp(y3) + pt4/dsqrt(s)*dexp(y4)
	x2c = pt3/dsqrt(s)*dexp(-y3) + pt4/dsqrt(s)*dexp(-y4)
	if (x1c.ge.un.or.x2c.ge.un) then
	  f4dimt4 = 0.d0
	  return
	endif
c ------------------ y5 -------------------------------------------
	y5max = dlog(dsqrt(s)/pt5*(1.d0-x1c))
	y5min = -dlog(dsqrt(s)/pt5*(1.d0-x2c))
	y5 = y4 + rr*dcos(omega)
	if (y5.gt.y5max.or.y5.lt.y5min) then
	   f4dimt4 = 0.d0
	   return
	endif
c ------------------ x -------------------------------------------
	x1 = x1c + pt5/dsqrt(s)*dexp(y5)
	x2 = x2c + pt5/dsqrt(s)*dexp(-y5)
	if (x1.ge.un.or.x2.ge.un) then
	  f4dimt4 = 0.d0
	  return
	endif
c coupure en angle-------------------------------------------------
        r45s = (y4-y5)**2+fi45**2
        rs = r**2
	icorr = 0
        if (r45s.lt.rs) then
           f4dimt4 = 0.d0
           return
**            icorr = 1
        endif
c coupure en isolement-------------------------------------------------
	call isolement(y4,y5,fi45,x4,pt5,gpt4,r_isol,etmax,
     #	isol4)
	cos_fi35 = -(pt4*dcos(fi45)+pt5)/pt3
	fi35 = zacos(cos_fi35)
	call isolement(y3,y5,fi35,x3,pt5,gpt3,r_isol,etmax,
     #	isol3)
** 	if (.not.(isol4)) then
** 	  f4dimt4 = 0.d0
** 	  return
** 	endif
	xisol3 = 1.d0
	if (.not.(isol3).or..not.(isol4)) then
	  if (icorr.eq.0) then
	    f4dimt4 = 0.d0
	    return
	  else if (icorr.eq.1) then
	    xisol3 = 0.d0
	  endif
	endif
c ------------------ fin isolement ------------------------------------
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu4max-uu4min)*
     #	rr*(rrmax-rrmin)*(omegamax-omegamin)*(uu3max-uu3min)*
     #	(pt5max-pt5min)*(x4max-x4min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
	pt3p = (pt4+pt5)
	x1s = pt3p/dsqrt(s)*(dexp(y3)+dexp(y4))
	x2s = pt3p/dsqrt(s)*(dexp(-y3)+dexp(-y4))
	x3p = gpt3/pt3p
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
	p4(1) = gpt4*dcosh(y4)
	p4(2) = gpt4
	p4(3) = 0.d0
	p4(4) = gpt4*dsinh(y4)
c ------------------------------------------- 
	pp4(1) = pt4*dcosh(y4)
	pp4(2) = pt4
	pp4(3) = 0.d0
	pp4(4) = pt4*dsinh(y4)
c ------------------------------------------- 
	p5(1) = pt5*dcosh(y5)
	p5(2) = pt5*dcos(fi45)
	p5(3) = pt5*dsin(fi45)
	p5(4) = pt5*dsinh(y5)
c ------------------------------------------- 
	pp3(1) = p1(1) + p2(1) - pp4(1) - p5(1)
	pp3(2) = p1(2) + p2(2) - pp4(2) - p5(2)
	pp3(3) = p1(3) + p2(3) - pp4(3) - p5(3)
	pp3(4) = p1(4) + p2(4) - pp4(4) - p5(4)
c -------------------------------------------
	p3(1) = x3*pp3(1)
	p3(2) = x3*pp3(2)
	p3(3) = x3*pp3(3)
	p3(4) = x3*pp3(4)
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
	p5p(2) = pt5
	p5p(3) = 0.d0
	p5p(4) = pt5*dsinh(y4)
c ------------------------------------------- 
	pp3p(1) = p1p(1) + p2p(1) - pp4(1) - p5p(1)
	pp3p(2) = p1p(2) + p2p(2) - pp4(2) - p5p(2)
	pp3p(3) = p1p(3) + p2p(3) - pp4(3) - p5p(3)
	pp3p(4) = p1p(4) + p2p(4) - pp4(4) - p5p(4)
c -------------------------------------------
	p3p(1) = x3p*pp3p(1)
	p3p(2) = x3p*pp3p(2)
	p3p(3) = x3p*pp3p(3)
	p3p(4) = x3p*pp3p(4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	mass_cut = 1.d0
	if (iphocut.eq.0) then
c angle azymuthal entre les deux photons plus grand que fimin
	  abspt3 = dsqrt(p3(2)**2+p3(3)**2)
	  abspt4 = dsqrt(p4(2)**2+p4(3)**2)
	  fi34 = zacos((p3(2)*p4(2)+p3(3)*p4(3))/(abspt3*abspt4))
	  if (fi34.le.fimin) then
	    f4dimt4 = 0.d0
	    return
	  endif
	else if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    if (icorr.eq.0) then
	      f4dimt4 = 0.d0
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
	    f4dimt4 = 0.d0
	    return
	  endif
	endif	
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	call fspcout(n,spcou)
c
c c_{ij} = alpha_s^3/(4*c_i*c_j*pi*s^2)
	cij = alfas(iloop,mu*mu)**3/(pi*s*s)
c
	call strfrat(x1,ih1,x2,ih2,x3,ih3,x4,ih4,sf0)
cc	call ampt(4,p1,p2,pp3,pp4,p5,camp0)
	call ampt4(s,x1,x2,y4,pt4,y5,pt5,fi45,camp0)
	call vec_dmult_vector(camp0,sf0,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f4dimt4 = 0.d0
	do i = j_processt_min,j_processt_max
	  f4dimt4 = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	  temp1(i+63) + f4dimt4
	enddo
	f4dimt4 = cij*f4dimt4/pt3
	if (icorr.eq.1) then
	  mass_cutp = 1.d0
	  if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	    mgg2 = 2.d0*sca(p3p,p4)
	    mgg = dsqrt(mgg2)
	    if (mgg.lt.mmin.or.mgg.gt.mmax) then
	      mass_cutp = 0.d0
	    endif
	  endif	
c coupure en isolement-------------------------------------------------
	  call isolement(y4,y4,0.d0,x4,pt5,gpt4,r_isol,etmax,
     #	  isol4p)
	  call isolement(y3,y4,pi,x3p,pt5,gpt3,r_isol,etmax,
     #	  isol3p)
	  xisol3p = 1.d0
	  if (.not.(isol3p).or..not.(isol4p)) then
	    xisol3p = 0.d0
	  endif
c ------------------ fin isolement ------------------------------------
	  call choiscale(p3p,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
c
c c_{ij} = alpha_s^3/(4*c_i*c_j*pi*s^2)
	  cij = alfas(iloop,mu*mu)**3/(pi*s*s)
c
	  call strfrat(x1s,ih1,x2s,ih2,x3p,ih3,x4,ih4,sf0)
	  call ampt_corr4(s,x1s,x2s,y4,pt4,pt5,r45s,camp0)
	  call vec_dmult_vector(camp0,sf0,k0max,temp0)
	  call vec_dmult_vector(temp0,spcou,k0max,temp1)
	  f4dimt4p = 0.d0
	  do i = j_processt_min,j_processt_max
	    f4dimt4p = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	    temp1(i+63) + f4dimt4p
	  enddo
	  f4dimt4p = cij*f4dimt4p/pt3p
	  f4dimt4 = f4dimt4*mass_cut*xisol3 - 
     #	  f4dimt4p*mass_cutp*xisol3p
	endif
c
	f4dimt4 = f4dimt4 * xjacob * (pt4/x4*pt5) * hc2
	if (facteur_symetrie) then
	  f4dimt4 = f4dimt4/2.d0
	endif	
	if (ispring.eq.1) then
	  half = 0.5d0
	  if (rnumb.le.half) then
	    xfi12 = 2.d0*pi-zacos(p3(2)/gpt3)
	    xfi13 = fi45
	  else
	    xfi12 = zacos(p3(2)/gpt3)
	    xfi13 = 2.d0*pi-fi45
	  endif
	  xpt1 = gpt4
	  xy1 = y4
	  xpt2 = gpt3
	  xy2 = y3
	  xpt3 = pt5
	  xy3 = y5
	  xfrag1 = x4
	  xfrag2 = x3
	  resfunc = f4dimt4
	  f4dimt4 = dabs(f4dimt4)
	endif
c
	return
	end 
c*****************************************************
c integrale a 2 dimensions
	double precision function f2dimt1(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/faconv/hc2
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension sv2(k0max),sy3(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	expos = 7.d0
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
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(3)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(4)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
c ------------------ x3 -------------------------------------------
c pour l'isolement, on demande que x3 >= gpt3/(gpt3+etmax) et
c x4 >= gpt4/(gpt4+etmax) mais a cause de la conservation de
c l'impulsion transverse x4 = gpt4/gpt3*x3 d'ou
c x3 doit etre >= gpt3/(gpt4+etmax)
	x3max = dmin1(1.d0,gpt3/gpt4)
	x3min1_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3min2_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y4)
	if (i_flag_iso.eq.1) then
	  x3min1_isol = gpt3/(gpt3+etmax)
	  x3min2_isol = gpt3/(gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  x3min1_isol = 1.d0/(1.d0+etmax)
	  x3min2_isol = gpt3/(gpt4*(1.d0+etmax))
	endif
	x3minp_incl = dmax1(x3min1_incl,x3min2_incl)
	x3minp_isol = dmax1(x3min1_isol,x3min2_isol)
	x3minp = dmax1(x3minp_incl,x3minp_isol)
	x3min = dmin1(x3max,x3minp)
	x3 = x3min + (x3max-x3min)*xx(5)
	pt3 = gpt3/x3
c ------------------ u -------------------------------------------
	u2 = xx(6)
	x4 = gpt4/gpt3*x3
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu3max-uu3min)
     #	*(uu4max-uu4min)*(x3max-x3min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
	call fspcout(n,spcou)
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
	    f2dimt1 = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s^2/(4*c_i*c_j*s**2)
	cbij = alphas**2*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	call ssy3t(s,gpt3,y3,y4,m,x3,x4,u2,ih1,ih2,ih3,ih4,sy3)
	call ssv2t(s,gpt3,y3,y4,x3,x4,u2,ih1,ih2,ih3,ih4,sv2)
	call vec_dadd_vector(sy3,sv2,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dim = 0.d0
	do i = j_processt_min,j_processt_max
	  f2dim = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	  temp1(i+63) + f2dim
	enddo
cc	f2dim = vec_dsum(temp0,k0max)
c
	f2dimt1 = f2dim * xjacob * alpi * cbij * hc2
	if (facteur_symetrie) then
	  f2dimt1 = f2dimt1/2.d0
	endif	
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
	  xfrag2 = x4
	  resfunc = f2dimt1
	  f2dimt1 = dabs(f2dimt1)
	endif
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimt2(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/faconv/hc2
	common/flagexp/icoup_exp
	common/coupexp/ptm_exp,r_exp
	common/coup/ptmm,r
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension sy3r(k0max),sy3p(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max)
	dimension temp3(k0max),temp4(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	zero = 0.d0
	un = 1.d0
	expos = 7.d0
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
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)),
     #	dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(3)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(4)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
	call fspcout(n,spcou)
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
c ------------------ x3 -------------------------------------------
	x3minpp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3max = dmin1(gpt3/gpt4,un)
	x3min = dmin1(x3max,x3minpp)
	x3 = x3min + (x3max-x3min)*xx(5)
	pt3 = gpt3/x3
	x10 = pt3*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = pt3*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimt2 = 0.d0
	  return
	endif
c ------------- z3 ----------------------------------
	z3sup = pt3/(pt3+ptm)
	if (i_flag_iso.eq.1) then
	  z3_isol = dmax1(gpt3/((gpt3+etmax)*x3),gpt3/((gpt4+etmax)*x3))
	else if (i_flag_iso.eq.2) then
	  z3_isol = dmax1(1.d0/((1.d0+etmax)*x3),
     #              gpt3/((gpt4*(1.d0+etmax))*x3))
	endif
	if (z3_isol.ge.un) then
	  f2dimt2 = 0.d0
	  return
	endif
	z3min1 = dmax1(x3minpp/x3,x10,x20,z3_isol)
	if (icoup_exp.eq.0) then
	  z3min = dmin1(z3min1,1.d0)	
	  z3max = un	
	else if (icoup_exp.eq.1) then
	  z3min = dmin1(z3min1,1.d0)	
	  z3max = z3sup	
	endif
	z3 = z3min + (z3max-z3min)*xx(6)
	x4 = gpt4/gpt3*x3*z3
c
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu3max-uu3min)
     #	*(uu4max-uu4min)*(x3max-x3min)*(z3max-z3min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimt2 = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s^2/(4*c_i*c_j*s**2)
	cbij = alphas**2*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	if (z3.le.z3sup) then
	  call ssy3rt(s,gpt3,y3,y4,x3,x4,z3,ih1,ih2,ih3,ih4,sy3r)
	else
	  call vec_dinit(sy3r,k0max,zero)	  
	endif
	call ssy3pt(s,gpt3,gpt4,y3,y4,mf,x3,x4,z3,ih1,ih2,ih3,
     #	ih4,sy3p)
	call vec_dadd_vector(sy3p,sy3r,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
c
	f2dim = 0.d0
	do i = j_processt_min,j_processt_max
	  f2dim = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	  temp1(i+63) + f2dim
	enddo
cc	f2dim = vec_dsum(temp4,k0max)
c
	f2dimt2 = f2dim * xjacob * alpi * cbij * hc2
	if (facteur_symetrie) then
	  f2dimt2 = f2dimt2/2.d0
	endif	
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
	  xfrag2 = x4
	  resfunc = f2dimt2
	  f2dimt2 = dabs(f2dimt2)
	endif
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimt2p(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/faconv/hc2
	common/flagexp/icoup_exp
	common/coupexp/ptm_exp,r_exp
	common/coup/ptmm,r
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension sy3r(k0max),sy3p(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max)
	dimension temp3(k0max),temp4(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	zero = 0.d0
	un = 1.d0
	expos = 7.d0
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
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)),
     #	dsqrt(s)/(2.d0*dcosh(y3)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(3)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)),gpt4)
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(4)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
	call fspcout(n,spcou)
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
c ------------------ x3 -------------------------------------------
	x3minpp = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3minp = dmax1(x3minpp,gpt3/gpt4)
	x3max = 1.d0
	x3min = dmin1(x3max,x3minp)
	x3 = x3min + (x3max-x3min)*xx(5)
	if (x3min.eq.x3max) then
	  f2dimt2p = 0.d0
	  return
	endif
	pt3 = gpt3/x3
	x10 = pt3*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = pt3*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimt2p = 0.d0
	  return
	endif
c ------------- z3 ----------------------------------
	z3sup = pt3/(pt3+ptm)
	if (i_flag_iso.eq.1) then
	  z3_isol = dmax1(gpt3/((gpt3+etmax)*x3),gpt3/((gpt4+etmax)*x3))
	else if (i_flag_iso.eq.2) then
	  z3_isol = dmax1(1.d0/((1.d0+etmax)*x3),
     #              gpt3/((gpt4*(1.d0+etmax))*x3))
	endif
	if (z3_isol.ge.un) then
	  f2dimt2p = 0.d0
	  return
	endif
	z3min1 = dmax1(x3minpp/x3,x10,x20,z3_isol)
	z3maxp = gpt3/gpt4/x3
	if (icoup_exp.eq.0) then
	  z3min = dmin1(z3min1,z3maxp)	
	  z3max = z3maxp	
	else if (icoup_exp.eq.1) then
	  z3min = dmin1(z3min1,1.d0)	
	  z3max = z3sup	
	endif
	z3 = z3min + (z3max-z3min)*xx(6)
	x4 = gpt4/gpt3*x3*z3
c
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu3max-uu3min)
     #	*(uu4max-uu4min)*(x3max-x3min)*(z3max-z3min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimt2p = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s^2/(4*c_i*c_j*s**2)
	cbij = alphas**2*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	if (z3.le.z3sup) then
	  call ssy3rt(s,gpt3,y3,y4,x3,x4,z3,ih1,ih2,ih3,ih4,sy3r)
	else
	  call vec_dinit(sy3r,k0max,zero)	  
	endif
	call ssy3ptp(s,gpt3,gpt4,y3,y4,mf,x3,x4,z3,ih1,ih2,ih3,
     #	ih4,sy3p)
	call vec_dadd_vector(sy3p,sy3r,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
c
	f2dim = 0.d0
	do i = j_processt_min,j_processt_max
	  f2dim = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	  temp1(i+63) + f2dim
	enddo
cc	f2dim = vec_dsum(temp4,k0max)
c
	f2dimt2p = f2dim * xjacob * alpi * cbij * hc2
	if (facteur_symetrie) then
	  f2dimt2p = f2dimt2p/2.d0
	endif	
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
	  xfrag2 = x4
	  resfunc = f2dimt2p
	  f2dimt2p = dabs(f2dimt2p)
	endif
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimt3(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/faconv/hc2
	common/flagexp/icoup_exp
	common/coupexp/ptm_exp,r_exp
	common/coup/ptmm,r
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension sy4r(k0max),sy4p(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max)
	dimension temp3(k0max),temp4(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	zero = 0.d0
	un = 1.d0
	expos = 7.d0
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
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)),
     #	dsqrt(s)/(2.d0*dcosh(y4)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(3)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(4)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
	call fspcout(n,spcou)
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
c ------------------ x4 -------------------------------------------
	x4minpp = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	x4max = dmin1(gpt4/gpt3,un)
	x4min = dmin1(x4max,x4minpp)
	x4 = x4min + (x4max-x4min)*xx(5)
	pt4 = gpt4/x4
	x10 = pt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = pt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimt3 = 0.d0
	  return
	endif
c ------------- z4 ----------------------------------
	z4sup = pt4/(pt4+ptm)
	if (i_flag_iso.eq.1) then
	  z4_isol = dmax1(gpt4/((gpt4+etmax)*x4),gpt4/((gpt3+etmax)*x4))
	else if (i_flag_iso.eq.2) then
	  z4_isol = dmax1(1.d0/((1.d0+etmax)*x4),
     #              gpt4/((gpt3*(1.d0+etmax))*x4))
	endif
	if (z4_isol.ge.un) then
	  f2dimt3 = 0.d0
	  return
	endif
	z4min1 = dmax1(x4minpp/x4,x10,x20,z4_isol)
	if (icoup_exp.eq.0) then
	  z4min = dmin1(z4min1,1.d0)	
	  z4max = un	
	else if (icoup_exp.eq.1) then
	  z4min = dmin1(z4min1,1.d0)	
	  z4max = z4sup	
	endif
	z4 = z4min + (z4max-z4min)*xx(6)
	x3 = gpt3/gpt4*x4*z4
c
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu3max-uu3min)
     #	*(uu4max-uu4min)*(x4max-x4min)*(z4max-z4min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimt3 = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s^2/(4*c_i*c_j*s**2)
	cbij = alphas**2*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	if (z4.le.z4sup) then
	  call ssy4rt(s,gpt4,y3,y4,x3,x4,z4,ih1,ih2,ih3,ih4,sy4r)
	else
	  call vec_dinit(sy4r,k0max,zero)	  
	endif
	call ssy4pt(s,gpt3,gpt4,y3,y4,mf,x3,x4,z4,ih1,ih2,ih3,
     #	ih4,sy4p)
	call vec_dadd_vector(sy4p,sy4r,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
c
	f2dim = 0.d0
	do i = j_processt_min,j_processt_max
	  f2dim = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	  temp1(i+63) + f2dim
	enddo
cc	f2dim = vec_dsum(temp4,k0max)
c
	f2dimt3 = f2dim * xjacob * alpi * cbij * hc2
	if (facteur_symetrie) then
	  f2dimt3 = f2dimt3/2.d0
	endif	
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
	  xfrag2 = x4
	  resfunc = f2dimt3
	  f2dimt3 = dabs(f2dimt3)
	endif
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimt3p(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/faconv/hc2
	common/flagexp/icoup_exp
	common/coupexp/ptm_exp,r_exp
	common/coup/ptmm,r
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension sy4r(k0max),sy4p(k0max),spcou(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max)
	dimension temp3(k0max),temp4(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	zero = 0.d0
	un = 1.d0
	expos = 7.d0
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
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)),
     #	dsqrt(s)/(2.d0*dcosh(y4)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(3)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)),gpt3)
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(4)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
	call fspcout(n,spcou)
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
c ------------------ x4 -------------------------------------------
	x4minpp = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	x4minp = dmax1(x4minpp,gpt4/gpt3)
	x4max = 1.d0
	x4min = dmin1(x4max,x4minp)
	x4 = x4min + (x4max-x4min)*xx(5)
	if (x4min.eq.x4max) then
	  f2dimt3p = 0.d0
	  return
	endif
	pt4 = gpt4/x4
	x10 = pt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = pt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimt3p = 0.d0
	  return
	endif
c ------------- z4 ----------------------------------
	z4sup = pt4/(pt4+ptm)
	if (i_flag_iso.eq.1) then
	  z4_isol = dmax1(gpt4/((gpt4+etmax)*x4),gpt4/((gpt3+etmax)*x4))
	else if (i_flag_iso.eq.2) then
	  z4_isol = dmax1(1.d0/((1.d0+etmax)*x4),
     #              gpt4/((gpt3*(1.d0+etmax))*x4))
	endif
	if (z4_isol.ge.un) then
	  f2dimt3p = 0.d0
	  return
	endif
	z4min1 = dmax1(x4minpp/x4,x10,x20,z4_isol)
	z4maxp = gpt4/gpt3/x4
	if (icoup_exp.eq.0) then
	  z4min = dmin1(z4min1,z4maxp)	
	  z4max = z4maxp	
	else if (icoup_exp.eq.1) then
	  z4min = dmin1(z4min1,1.d0)	
	  z4max = z4sup	
	endif
	z4 = z4min + (z4max-z4min)*xx(6)
	x3 = gpt3/gpt4*x4*z4
c
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu3max-uu3min)
     #	*(uu4max-uu4min)*(x4max-x4min)*(z4max-z4min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
c -------------coupure entre les 2 photons-------------------------
	if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    f2dimt3p = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s^2/(4*c_i*c_j*s**2)
	cbij = alphas**2*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c
	if (z4.le.z4sup) then
	  call ssy4rt(s,gpt4,y3,y4,x3,x4,z4,ih1,ih2,ih3,ih4,sy4r)
	else
	  call vec_dinit(sy4r,k0max,zero)	  
	endif
	call ssy4ptp(s,gpt3,gpt4,y3,y4,mf,x3,x4,z4,ih1,ih2,ih3,
     #	ih4,sy4p)
	call vec_dadd_vector(sy4p,sy4r,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
c
	f2dim = 0.d0
	do i = j_processt_min,j_processt_max
	  f2dim = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	  temp1(i+63) + f2dim
	enddo
cc	f2dim = vec_dsum(temp4,k0max)
c
	f2dimt3p = f2dim * xjacob * alpi * cbij * hc2
	if (facteur_symetrie) then
	  f2dimt3p = f2dimt3p/2.d0
	endif	
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
	  xfrag2 = x4
	  resfunc = f2dimt3p
	  f2dimt3p = dabs(f2dimt3p)
	endif
	return
	end 
c********************************************************************
c integrale a 1 dimension
	double precision function f1dimt(xx)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical facteur_symetrie
	common/hadron/ih1,ih2,ih3,ih4
	common/scale/m,mf,mu
	common/coul/n,cf,gtr
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/faconv/hc2
	common/born/iborn
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/sfragv/xfrag1,xfrag2
	common/baspring/ispring
	common/nboucle/iloop
	common/processt/j_processt_min,j_processt_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/symetrie/facteur_symetrie
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=84)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	dimension c0sv12(k0max),c0sv3(k0max),c0sv4(k0max),c0fborn(k0max)
	dimension spcou(k0max)
	dimension p3(4),p4(4),p5(4)
	pi=datan(1.d0)*4.d0
	un = 1.d0
	expos = 7.d0
c ------------------ y3 -------------------------------------------
	y3max = ysup
	y3min = yinf
	y3 = y3min + (y3max-y3min)*xx(1)
c ------------------ y4 -------------------------------------------
	y4max = ysup
	y4min = yinf
	y4 = y4min + (y4max-y4min)*xx(2)
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	uu4max = dexp(-expos*dlog(gpt4min))
	uu4min = dexp(-expos*dlog(gpt4max))
	uu4 = uu4min + (uu4max-uu4min)*xx(3)
	gpt4 = dexp(-dlog(uu4)/expos)
	gpt4_jac = dexp(-(expos+1.d0)*dlog(uu4)/expos)/expos
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	uu3max = dexp(-expos*dlog(gpt3min))
	uu3min = dexp(-expos*dlog(gpt3max))
	uu3 = uu3min + (uu3max-uu3min)*xx(4)
	gpt3 = dexp(-dlog(uu3)/expos)
	gpt3_jac = dexp(-(expos+1.d0)*dlog(uu3)/expos)/expos
c ------------------ x3 -------------------------------------------
c pour l'isolement, on demande que x3 >= gpt3/(gpt3+etmax) et
c x4 >= gpt4/(gpt4+etmax) mais a cause de la conservation de
c l'impulsion transverse x4 = gpt4/gpt3*x3 d'ou
c x3 doit etre >= gpt3/(gpt4+etmax)
	x3max = dmin1(1.d0,gpt3/gpt4)
	x3min1_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	x3min2_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y4)
	if (i_flag_iso.eq.1) then
	  x3min1_isol = gpt3/(gpt3+etmax)
	  x3min2_isol = gpt3/(gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  x3min1_isol = 1.d0/(1.d0+etmax)
	  x3min2_isol = gpt3/(gpt4*(1.d0+etmax))
	endif
	x3minp_incl = dmax1(x3min1_incl,x3min2_incl)
	x3minp_isol = dmax1(x3min1_isol,x3min2_isol)
	x3minp = dmax1(x3minp_incl,x3minp_isol)
	x3min = dmin1(x3max,x3minp)
	x3 = x3min + (x3max-x3min)*xx(5)
	pt3 = gpt3/x3
	x4 = gpt4/gpt3*x3
c
	x10 = pt3*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = pt3*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f1dimt = 0.d0
	  return
	endif
c
	xjacob = (y3max-y3min)*(y4max-y4min)*(uu3max-uu3min)
     #	*(uu4max-uu4min)*(x3max-x3min)
	xjacob = xjacob * gpt3_jac * gpt4_jac
	call fspcout(n,spcou)
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
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
	    f1dimt = 0.d0
	    return
	  endif
	endif
c set the scales-------------------------------------------------
	call choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	alphas = alfas(iloop,mu*mu)
c c^b_{ij} = 2*pi*alpha_s^2/(4*c_i*c_j*s**2)
	cbij = alphas**2*2.d0*pi/(s*s)
	alpi = alphas/(2.d0*pi)
c 
	if (i_flag_iso.eq.1) then
	  z3_isol = dmax1(gpt3/((gpt3+etmax)*x3),gpt3/((gpt4+etmax)*x3))
	  z4_isol = dmax1(gpt4/((gpt4+etmax)*x4),gpt4/((gpt3+etmax)*x4))
	else if (i_flag_iso.eq.2) then
	  z3_isol = dmax1(1.d0/((1.d0+etmax)*x3),
     #              gpt3/((gpt4*(1.d0+etmax))*x3))
	  z4_isol = dmax1(1.d0/((1.d0+etmax)*x4),
     #              gpt4/((gpt3*(1.d0+etmax))*x4))
	endif
	if ((z3_isol.ge.un).or.(z4_isol.ge.un)) then
	  f1dimt = 0.d0
	  return
	endif
c 
	if (iborn.eq.0) then 
	  call ssv12t(s,gpt3,y3,y4,m,mu,x3,x4,x10,x20,
     #	  ih1,ih2,ih3,ih4,c0sv12)
	  call ssv3t(s,gpt3,y3,mf,x3,x4,x10,x20,z3_isol,ih1,ih2,
     #	  ih3,ih4,c0sv3)
	  call ssv4t(s,gpt4,y4,mf,x3,x4,x10,x20,z4_isol,ih1,ih2,
     #	  ih3,ih4,c0sv4)
	  call vec_dmult_vector(c0sv12,spcou,k0max,temp0)
	  call vec_dmult_vector(c0sv3,spcou,k0max,temp1)
	  call vec_dmult_vector(c0sv4,spcou,k0max,temp2)
	  sv12 = 0.d0
	  sv3 = 0.d0
	  sv4 = 0.d0
	  do i = j_processt_min,j_processt_max
	    sv12 = temp0(i) + temp0(i+21) + temp0(i+42) + 
     #	    temp0(i+63) + sv12
	    sv3 = temp1(i) + temp1(i+21) + temp1(i+42) + 
     #	    temp1(i+63) + sv3
	    sv4 = temp2(i) + temp2(i+21) + temp2(i+42) + 
     #	    temp2(i+63) + sv4
	  enddo
cc	  sv12 = vec_dsum(temp0,k0max)
cc	  sv3 = vec_dsum(temp1,k0max)
cc	  sv4 = vec_dsum(temp2,k0max)
	  fborn = 0.d0
	else if (iborn.eq.1) then
	  call sfbornt(s,n,gpt3,y3,y4,x3,x4,ih1,ih2,ih3,ih4,c0fborn)
	  call vec_dmult_vector(c0fborn,spcou,k0max,temp3)
	  fborn = 0.d0
	  do i = j_processt_min,j_processt_max
	    fborn = temp3(i) + temp3(i+21) + temp3(i+42) + 
     #	    temp3(i+63) + fborn
	  enddo
cc	  fborn = vec_dsum(temp3,k0max)
	  sv12 = 0.d0
	  sv3 = 0.d0
	  sv4 = 0.d0
	endif
c 
	resb = fborn + (sv12+sv3+sv4)*alpi
	f1dimt = resb * xjacob * cbij * hc2
	if (facteur_symetrie) then
	  f1dimt = f1dimt/2.d0
	endif	
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
	  xfrag2 = x4
	  resfunc = f1dimt
	  f1dimt = dabs(f1dimt)
	endif
c
	return
	end 
c===================================================================
	subroutine ampt3(s,x1,x2,y3,pt3,y5,pt5,fi35,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=84)
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
	call sh12t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h12)
	call sh13t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h13)
	call sh23t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h23)
	call sh34t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	call sconst(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,cons)
	do i=1,k0max
	  camp(i) = (0.5d0*h12(i)*e12+h13(i)*e13
     #	             +h23(i)*e23+h34(i)*e34p+0.5d0*cons(i))
	enddo
	return
	end
c===================================================================
	subroutine ampt_corr3(s,x1,x2,y3,pt3,pt5,rs,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=84)
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
	call sh13t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h13)
	call sh23t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h23)
	call sh34t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	do i=1,k0max
	  camp(i) = (h13(i)+h23(i)+h34(i))*2.d0/rs/(pt5*pt5)
	enddo
	return
	end
c===================================================================
	subroutine ampt4(s,x1,x2,y4,pt4,y5,pt5,fi45,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=84)
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
	call sh12t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h12)
	call sh14t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h14)
	call sh24t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h24)
	call sh34t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	call sconst(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,cons)
	do i=1,k0max
	  camp(i) = (0.5d0*h12(i)*e12+h14(i)*e14
     #	             +h24(i)*e24+h34(i)*e34s+0.5d0*cons(i))
	enddo
	return
	end
c===================================================================
	subroutine ampt_corr4(s,x1,x2,y4,pt4,pt5,rs,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=84)
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
	call sh14t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h14)
	call sh24t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h24)
	call sh34t(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	do i=1,k0max
	  camp(i) = (h14(i)+h24(i)+h34(i))*2.d0/rs/(pt5*pt5)
	enddo
	return
	end
c********************************************************************
