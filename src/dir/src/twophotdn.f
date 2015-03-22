c
c modification le 21/01/97 pour generer le fi_3 entre 0 et 2 pi
c modification le 22/01/97 pour la separation entre les deux photons
c modification le 14/02/97 pour mettre la coupure sur la masse
c invariante dans les quasi deux donne deux et les deux donne deux
c modification le 7/03/97 pour fixer le nombre total d'evenements au
c depart
c modification le 7/04/97 pour fixer un bug dans f2dim2 et f2dim3
c (modifications aussi dans dir.f: subroutine ssv3 et ssv4)
c modification le 9/06/97 pour l'integration dans f4dim3 et f4dim4
c modification le 3/05/99 pour ajout des termes en r^2*ln(p_{tm})
c modification le 14/05/99 pour introduire le maximum de la masse 
c invariante photon-photon
	subroutine dir_sub
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical lo,nlo,gener,integ,intega
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
	common/processd/j_processd_min,j_processd_max
**         common/sortier/ntrack,we(maxtrk),wpx(maxtrk),
**      #	wpy(maxtrk),wpz(maxtrk)
	common/randd/rnumb
	common/cheminbs/path_bsfile
	common/long/ilen
	common/efficacite/ipositif,inegatif
	dimension wrvec(irandom),wrvec1(irandom)
	external f4dimd3,f4dimd4,f2dimd1,f2dimd2,f2dimd3,f1dimd
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
	idim = 6
	iwild = 6
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimd3,resf34,sdf34,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimd3,resf34,sdf34,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimd4,resf44,sdf44,ctime,it1,it2)
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
	idim = 6
	iwild = 6
	do i=1,idim
	  xl(i) = 0.d0
	  xu(i) = 1.d0
	enddo
	call bases (f4dimd4,resf44,sdf44,ctime,it1,it2)
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
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimd1,resf21,sdf21,ctime,it1,it2)
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
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	call bases (f2dimd1,resf21,sdf21,ctime,it1,it2)
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
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	if (j_processd_max.eq.2) then
	  call bases (f2dimd2,resf22,sdf22,ctime,it1,it2)
	else
	  resf22 = 0.d0
	  sdf22 = 0.d0
	endif
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
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	if (j_processd_max.eq.2) then
	  call bases (f2dimd2,resf22,sdf22,ctime,it1,it2)
c on ecrit la grille
	  open(28,file=path_bsfile(1:ilen)//'qtwopart2.bs',
     #	  status='unknown',form='unformatted')
	  call bswrit(28)
	  close(28)
c	
	else
	  resf22 = 0.d0
	  sdf22 = 0.d0
	endif
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
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	if (j_processd_max.eq.2) then
	  call bases (f2dimd3,resf23,sdf23,ctime,it1,it2)
	else
	  resf23 = 0.d0
	  sdf23 = 0.d0
	endif
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
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	if (j_processd_max.eq.2) then
	  call bases (f2dimd3,resf23,sdf23,ctime,it1,it2)
c on ecrit la grille
	  open(28,file=path_bsfile(1:ilen)//'qtwopart3.bs',
     #	  status='unknown',form='unformatted')
	  call bswrit(28)
	  close(28)
c	
	else
	  resf23 = 0.d0
	  sdf23 = 0.d0
	endif
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
	idim = 3
	iwild = 3
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	if (j_processd_min.eq.1) then
	  call bases (f1dimd,resf1,sdf1,ctime,it1,it2)
	else
	  resf1 = 0.d0
	  sdf1 = 0.d0
	endif
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
	idim = 3
	iwild = 3
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	if (j_processd_min.eq.1) then
	  call bases (f1dimd,resf1,sdf1,ctime,it1,it2)
c on ecrit la grille
	  open(28,file=path_bsfile(1:ilen)//'twopart1.bs',
     #	  status='unknown',form='unformatted')
	  call bswrit(28)
	  close(28)
c	
	else
	  resf1 = 0.d0
	  sdf1 = 0.d0
	endif
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
	idim = 3
	iwild = 3
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	call bases (f1dimd,resfb,sdfb,ctime,it1,it2)
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
	  resf2 = resf21+resf22+resf23
	  sdf2 = sdf21+sdf22+sdf23
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
	    resf2 = resf21+resf22+resf23
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
	  write(28,*) resf21, resf22, resf23
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
	read(28,*) resf34, resf44
	read(28,*) resf21, resf22, resf23
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
	ispring = 2
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 6
	iwild = 6
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
	  call spring(f4dimd3,ixtry)
	  call doubletosingled(rnumb1,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c
c###
	ispring = 2
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall1
	idim = 6
	iwild = 6
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
	  call spring(f4dimd4,ixtry)
	  call doubletosingled(rnumb1,pi,poid,iprov,ntrack)
	  call gfill
	 enddo
c###
c
c###
	ispring = 2
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 4
	iwild = 4
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
	  call spring(f2dimd1,ixtry)
	  call doubletosingled(rnumb1,pi,poid,iprov,ntrack)
	  call gfill
	enddo
c###
c
c###
	ispring = 2
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	if (j_processd_max.eq.2) then
	  inbrealevent = int(resf22*xnorm)
c
	  open(28,file=path_bsfile(1:ilen)//'qtwopart2.bs',
     #	  status='unknown',form='unformatted')
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
            endif
	    rnumb1 = dble(wrvec1(i1+1))
	    call spring(f2dimd2,ixtry)
	    call doubletosingled(rnumb1,pi,poid,iprov,ntrack)
	    call gfill
	  enddo
	endif
c###
c
c###
	ispring = 2
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall
	idim = 4
	iwild = 4
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
	enddo
	if (j_processd_max.eq.2) then
	  inbrealevent = int(resf23*xnorm)
c
	  open(28,file=path_bsfile(1:ilen)//'qtwopart3.bs',
     #	  status='unknown',form='unformatted')
	  call bsread(28)
	  close(28)
c	
	write(6,*)
	write(6,*)'============================================='
	write(6,*)'event generation : quasi 2 partons : col. fi4'
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
	    call spring(f2dimd3,ixtry)
	    call doubletosingled(rnumb1,pi,poid,iprov,ntrack)
	    call gfill
	  enddo
	endif
c###
c
c###
	ispring = 2
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall2
	idim = 3
	iwild = 3
	do i=1,idim
          xl(i) = 0.d0
          xu(i) = 1.d0
        enddo
	if (j_processd_min.eq.1) then
	  inbrealevent = int(resf1*xnorm)
c
	  open(28,file=path_bsfile(1:ilen)//'twopart1.bs',
     #	  status='unknown',form='unformatted')
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
	    call spring(f1dimd,ixtry)
	    call doubletosingled(rnumb1,pi,poid,iprov,ntrack)
	    call gfill
	  enddo
	endif
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
	ispring = 2
	call bsinit
	acc1 = accu
	acc2 = accu
	itmx1 = jtmx1
	itmx2 = jtmx2
	icall = jcall2
	idim = 3
	iwild = 3
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
	  call spring(f1dimd,ixtry)
	  call doubletosingled(rnumb1,pi,poid,iprov,ntrack)
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
	double precision function f4dimd3(xx)
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
	common/baspring/ispring
	common/nboucle/iloop
	common/aurenche/iauren
	common/processd/j_processd_min,j_processd_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/randd/rnumb
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
	dimension xx(20)
	parameter (k0max=4)
	dimension sf0(k0max),camp0(k0max)
	dimension temp0(k0max),temp1(k0max),spcou(k0max)
	dimension p1(4),p2(4),p3(4),p4(4),p5(4)
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
c ------------------ pt5 --------------------------------------------
	pt5max = dsqrt(gpt3**2+gptsup**2)
	pt5min = ptm
	pt5 = pt5min + (pt5max-pt5min)*xx(4)
c ------------------ fi35 -------------------------------------------
	x11 = gpt3/dsqrt(s)*dexp(y3)+dabs(gpt3-pt5)/dsqrt(s)*dexp(y4)
	x21 = gpt3/dsqrt(s)*dexp(-y3)+dabs(gpt3-pt5)/dsqrt(s)*dexp(-y4)
	if (x11.gt.un.or.x21.gt.un) then
	  f4dimd3 = 0.d0
	  return
	endif
	sy5max = dlog(dsqrt(s)/pt5*(1.d0-x11))
	sy5min = -dlog(dsqrt(s)/pt5*(1.d0-x21))
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
	omega = omegamin + (omegamax-omegamin)*xx(5)
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
	rr = rrmin + (rrmax-rrmin)*xx(6)
	fi35 = rr*dsin(omega)
	if (fi35.gt.pi) then
	   f4dimd3 = 0.d0
	   return
	endif
c ------------------ gpt4 -------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	gpt4 = dsqrt(gpt3**2+pt5**2+2.d0*pt5*gpt3*dcos(fi35))
	if (gpt4.lt.gpt4min.or.gpt4.gt.gpt4max) then
	   f4dimd3 = 0.d0
	   return
	endif
c ------------------ x -------------------------------------------
	x1c = gpt3/dsqrt(s)*dexp(y3) + gpt4/dsqrt(s)*dexp(y4)
	x2c = gpt3/dsqrt(s)*dexp(-y3) + gpt4/dsqrt(s)*dexp(-y4)
	if (x1c.ge.un.or.x2c.ge.un) then
	  f4dimd3 = 0.d0
	  return
	endif
	sup_pt5 = dsqrt(s*(1.d0-x1c)*(1.d0-x2c))
	if (pt5.gt.sup_pt5) then
	  f4dimd3 = 0.d0
	  return
	endif	
c ------------------ y5 -------------------------------------------
	y5max = dlog(dsqrt(s)/pt5*(1.d0-x1c))
	y5min = -dlog(dsqrt(s)/pt5*(1.d0-x2c))
	y5 = y3 + rr*dcos(omega)
	if (y5.gt.y5max.or.y5.lt.y5min) then
	   f4dimd3 = 0.d0
	   return
	endif
c ------------------ x -------------------------------------------
	x1 = x1c + pt5/dsqrt(s)*dexp(y5)
	x2 = x2c + pt5/dsqrt(s)*dexp(-y5)
	if (x1.ge.un.or.x2.ge.un) then
	  f4dimd3 = 0.d0
	  return
	endif
c coupure en angle-------------------------------------------------
        r35s = (y3-y5)**2+fi35**2
        rs = r**2
	icorr = 0
        if (r35s.lt.rs) then
           f4dimd3 = 0.d0
           return
**            icorr = 1
        endif
c coupure en isolement-------------------------------------------------
	call isolement(y3,y5,fi35,un,pt5,gpt3,r_isol,etmax,
     #	isol3)
	cos_fi45 = -(gpt3*dcos(fi35)+pt5)/gpt4
	fi45 = zacos(cos_fi45)
	call isolement(y4,y5,fi45,un,pt5,gpt4,r_isol,etmax,
     #	isol4)
	xisol4 = 1.d0
	if (.not.(isol4).or..not.(isol3)) then
	  if (icorr.eq.0) then
	    f4dimd3 = 0.d0
	    return
	  else if (icorr.eq.1) then
	    xisol4 = 0.d0
	  endif
	endif
c ------------------ fin isolement ------------------------------------
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)*
     #	rr*(rrmax-rrmin)*(omegamax-omegamin)*
     #	(pt5max-pt5min)
	x1p = (gpt3+pt5)/dsqrt(s)*(dexp(y3)+dexp(y4))
	x2p = (gpt3+pt5)/dsqrt(s)*(dexp(-y3)+dexp(-y4))
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
	p5(1) = pt5*dcosh(y5)
	p5(2) = pt5*dcos(fi35)
	p5(3) = pt5*dsin(fi35)
	p5(4) = pt5*dsinh(y5)
c ------------------------------------------- 
	p4(1) = p1(1) + p2(1) - p3(1) - p5(1)
	p4(2) = p1(2) + p2(2) - p3(2) - p5(2)
	p4(3) = p1(3) + p2(3) - p3(3) - p5(3)
	p4(4) = p1(4) + p2(4) - p3(4) - p5(4)
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
	p4p(1) = p1p(1) + p2p(1) - p3(1) - p5p(1)
	p4p(2) = p1p(2) + p2p(2) - p3(2) - p5p(2)
	p4p(3) = p1p(3) + p2p(3) - p3(3) - p5p(3)
	p4p(4) = p1p(4) + p2p(4) - p3(4) - p5p(4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	mass_cut = 1.d0
	if (iphocut.eq.0) then
c angle azymuthal entre les deux photons plus grand que fimin
	  abspt3 = dsqrt(p3(2)**2+p3(3)**2)
	  abspt4 = dsqrt(p4(2)**2+p4(3)**2)
	  fi34 = zacos((p3(2)*p4(2)+p3(3)*p4(3))/(abspt3*abspt4))
	  if (fi34.le.fimin) then
	    f4dimd3 = 0.d0
	    return
	  endif
	else if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    if (icorr.eq.0) then
	      f4dimd3 = 0.d0
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
	    f4dimd3 = 0.d0
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
	call fspcoud(n,spcou)
c	
	alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s*alpha^2/(4*c_i*c_j*pi*s^2)
	cij = alfas(iloop,mu*mu)*alpha_em**2/(pi*s*s)
c
	call strfrad(x1,ih1,x2,ih2,sf0)
	call ampd3(s,x1,x2,y3,gpt3,y5,pt5,fi35,camp0)
	call vec_dmult_vector(camp0,sf0,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f4dimd3 = 0.d0
	do i = j_processd_min,j_processd_max
	  f4dimd3 = temp1(i) + temp1(2+i) + f4dimd3
	enddo
	f4dimd3 = cij*f4dimd3
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
	  call isolement(y3,y3,0.d0,un,pt5,gpt3,r_isol,etmax,
     #	  isol3p)
	  call isolement(y4,y3,pi,un,pt5,gpt3+pt5,r_isol,etmax,
     #	  isol4p)
	  xisol4p = 1.d0
	  if (.not.(isol4p).or..not.(isol3p)) then
	    xisol4p = 0.d0
	  endif
c ------------------ fin isolement ------------------------------------
	  if (iauren.eq.0) then
	    call choiscale(p3,p4p,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	  else if (iauren.eq.1) then
	    m = gpt3
	    mu = m
	    mf = m
	  endif
	  alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s*alpha^2/(4*c_i*c_j*pi*s^2)
	  cij = alfas(iloop,mu*mu)*alpha_em**2/(pi*s*s)
c
	  call strfrad(x1p,ih1,x2p,ih2,sf0)
	  call ampd_corr3(s,x1p,x2p,y3,gpt3,pt5,r35s,camp0)
	  call vec_dmult_vector(camp0,sf0,k0max,temp0)
	  call vec_dmult_vector(temp0,spcou,k0max,temp1)
	  f4dimd3p = 0.d0
	  do i = j_processd_min,j_processd_max
	    f4dimd3p = temp1(i) + temp1(2+i) + f4dimd3p
	  enddo
	  f4dimd3p = cij*f4dimd3p
	  f4dimd3 = f4dimd3*mass_cut*xisol4 - 
     #	  f4dimd3p*mass_cutp*xisol4p
	endif
c
	f4dimd3 = f4dimd3 * xjacob * (gpt3*pt5) * hc2/2.d0
c
	if (ispring.eq.1) then
	  f4dimd3 = dabs(f4dimd3)
	else if (ispring.eq.2) then
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
	  resfunc = f4dimd3
	  f4dimd3 = dabs(f4dimd3)
	endif
c
	return
	end 
c*****************************************************
c integrale a 4 dimensions
	double precision function f4dimd4(xx)
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
	common/baspring/ispring
	common/nboucle/iloop
	common/aurenche/iauren
	common/processd/j_processd_min,j_processd_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/randd/rnumb
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
	dimension xx(20)
	parameter (k0max=4)
	dimension sf0(k0max),camp0(k0max)
	dimension temp0(k0max),temp1(k0max),spcou(k0max)
	dimension p1(4),p2(4),p3(4),p4(4),p5(4)
	dimension p1p(4),p2p(4),p3p(4),p5p(4)
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
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(3)
c ------------------ pt5 --------------------------------------------
	pt5max = dsqrt(gpt4**2+gptsup**2)
	pt5min = ptm
	pt5 = pt5min + (pt5max-pt5min)*xx(4)
c ------------------ fi45 -------------------------------------------
	x11 = gpt4/dsqrt(s)*dexp(y4)+dabs(gpt4-pt5)/dsqrt(s)*dexp(y3)
	x21 = gpt4/dsqrt(s)*dexp(-y4)+dabs(gpt4-pt5)/dsqrt(s)*dexp(-y3)
	if (x11.gt.un.or.x21.gt.un) then
	  f4dimd4 = 0.d0
	  return
	endif
	sy5max = dlog(dsqrt(s)/pt5*(1.d0-x11))
	sy5min = -dlog(dsqrt(s)/pt5*(1.d0-x21))
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
	omega = omegamin + (omegamax-omegamin)*xx(5)
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
	rr = rrmin + (rrmax-rrmin)*xx(6)
	fi45 = rr*dsin(omega)
	if (fi45.gt.pi) then
	   f4dimd4 = 0.d0
	   return
	endif
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = dsqrt(gpt4**2+pt5**2+2.d0*pt5*gpt4*dcos(fi45))
	if (gpt3.lt.gpt3min.or.gpt3.gt.gpt3max) then
	   f4dimd4 = 0.d0
	   return
	endif
c ------------------ x -------------------------------------------
	x1c = gpt3/dsqrt(s)*dexp(y3) + gpt4/dsqrt(s)*dexp(y4)
	x2c = gpt3/dsqrt(s)*dexp(-y3) + gpt4/dsqrt(s)*dexp(-y4)
	if (x1c.ge.un.or.x2c.ge.un) then
	  f4dimd4 = 0.d0
	  return
	endif
	sup_pt5 = dsqrt(s*(1.d0-x1c)*(1.d0-x2c))
	if (pt5.gt.sup_pt5) then
	  f4dimd4 = 0.d0
	  return
	endif	
c ------------------ y5 -------------------------------------------
	y5max = dlog(dsqrt(s)/pt5*(1.d0-x1c))
	y5min = -dlog(dsqrt(s)/pt5*(1.d0-x2c))
	y5 = y4 + rr*dcos(omega)
	if (y5.gt.y5max.or.y5.lt.y5min) then
	   f4dimd4 = 0.d0
	   return
	endif
c ------------------ x -------------------------------------------
	x1 = x1c + pt5/dsqrt(s)*dexp(y5)
	x2 = x2c + pt5/dsqrt(s)*dexp(-y5)
	if (x1.ge.un.or.x2.ge.un) then
	  f4dimd4 = 0.d0
	  return
	endif
c coupure en angle-------------------------------------------------
        r45s = (y4-y5)**2+fi45**2
        rs = r**2
	icorr = 0
        if (r45s.lt.rs) then
           f4dimd4 = 0.d0
           return
**            icorr = 1
        endif
c coupure en isolement-------------------------------------------------
	cos_fi35 = -(gpt4*dcos(fi45)+pt5)/gpt3
	fi35 = zacos(cos_fi35)
	call isolement(y3,y5,fi35,un,pt5,gpt3,r_isol,etmax,
     #	isol3)
	call isolement(y4,y5,fi45,un,pt5,gpt4,r_isol,etmax,
     #	isol4)
	xisol3 = 1.d0
	if (.not.(isol3).or..not.(isol4)) then
	  if (icorr.eq.0) then
	    f4dimd4 = 0.d0
	    return
	  else if (icorr.eq.1) then
	    xisol3 = 0.d0
	  endif
	endif
c ------------------ fin isolement ------------------------------------
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt4max-gpt4min)*
     #	rr*(rrmax-rrmin)*(omegamax-omegamin)*
     #	(pt5max-pt5min)
	x1s = (gpt4+pt5)/dsqrt(s)*(dexp(y3)+dexp(y4))
	x2s = (gpt4+pt5)/dsqrt(s)*(dexp(-y3)+dexp(-y4))
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
	p5(1) = pt5*dcosh(y5)
	p5(2) = pt5*dcos(fi45)
	p5(3) = pt5*dsin(fi45)
	p5(4) = pt5*dsinh(y5)
c ------------------------------------------- 
	p3(1) = p1(1) + p2(1) - p4(1) - p5(1)
	p3(2) = p1(2) + p2(2) - p4(2) - p5(2)
	p3(3) = p1(3) + p2(3) - p4(3) - p5(3)
	p3(4) = p1(4) + p2(4) - p4(4) - p5(4)
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
	p3p(1) = p1p(1) + p2p(1) - p4(1) - p5p(1)
	p3p(2) = p1p(2) + p2p(2) - p4(2) - p5p(2)
	p3p(3) = p1p(3) + p2p(3) - p4(3) - p5p(3)
	p3p(4) = p1p(4) + p2p(4) - p4(4) - p5p(4)
c -------------------------------------------
c -------------coupure entre les 2 photons-------------------------
	mass_cut = 1.d0
	if (iphocut.eq.0) then
c angle azymuthal entre les deux photons plus grand que fimin
	  abspt3 = dsqrt(p3(2)**2+p3(3)**2)
	  abspt4 = dsqrt(p4(2)**2+p4(3)**2)
	  fi34 = zacos((p3(2)*p4(2)+p3(3)*p4(3))/(abspt3*abspt4))
	  if (fi34.le.fimin) then
	    f4dimd4 = 0.d0
	    return
	  endif
	else if (iphocut.eq.1) then
c masse invariante entre les deux photons plus grande que mmin
c et plus petite que mmax
	  mgg2 = 2.d0*sca(p3,p4)
	  mgg = dsqrt(mgg2)
	  if (mgg.lt.mmin.or.mgg.gt.mmax) then
	    if (icorr.eq.0) then
	      f4dimd4 = 0.d0
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
	    f4dimd4 = 0.d0
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
	call fspcoud(n,spcou)
c
	alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s*alpha^2/(4*c_i*c_j*pi*s^2)
	cij = alfas(iloop,mu*mu)*alpha_em**2/(pi*s*s)
c
	call strfrad(x1,ih1,x2,ih2,sf0)
	call ampd4(s,x1,x2,y4,gpt4,y5,pt5,fi45,camp0)
	call vec_dmult_vector(camp0,sf0,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f4dimd4 = 0.d0
	do i = j_processd_min,j_processd_max
	  f4dimd4 = temp1(i) + temp1(2+i) + f4dimd4
	enddo
	f4dimd4 = f4dimd4*cij
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
	  call isolement(y3,y4,pi,un,pt5,gpt3+pt5,r_isol,etmax,
     #	  isol3p)
	  call isolement(y4,y4,0.d0,un,pt5,gpt4,r_isol,etmax,
     #	  isol4p)
	  xisol3p = 1.d0
	  if (.not.(isol3p).or..not.(isol4p)) then
	    xisol3p = 0.d0
	  endif
c ------------------ fin isolement ------------------------------------
	  if (iauren.eq.0) then
	    call choiscale(p3p,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	  else if (iauren.eq.1) then
	    m = gpt4+pt5
	    mu = m
	    mf = m
	  endif
	  alpha_em = alphaem(iloopem,m*m)
c c_{ij} = alpha_s*alpha^2/(4*c_i*c_j*pi*s^2)
	  cij = alfas(iloop,mu*mu)*alpha_em**2/(pi*s*s)
c
	  call strfrad(x1s,ih1,x2s,ih2,sf0)
	  call ampd_corr4(s,x1s,x2s,y4,gpt4,pt5,r45s,camp0)
	  call vec_dmult_vector(camp0,sf0,k0max,temp0)
	  call vec_dmult_vector(temp0,spcou,k0max,temp1)
	  f4dimd4p = 0.d0
	  do i = j_processd_min,j_processd_max
	    f4dimd4p = temp1(i) + temp1(2+i) + f4dimd4p
	  enddo
	  f4dimd4p = cij*f4dimd4p
	  f4dimd4 = f4dimd4*mass_cut*xisol3 - 
     #	  f4dimd4p*mass_cutp*xisol3p
	endif
c
	f4dimd4 = f4dimd4 * xjacob * (gpt4*pt5) * hc2/2.d0
c
	if (ispring.eq.1) then
	  f4dimd4 = dabs(f4dimd4)
	else if (ispring.eq.2) then
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
	  resfunc = f4dimd4
	  f4dimd4 = dabs(f4dimd4)
	endif
c
	return
	end 
c********************************************************************
c integrale a 2 dimensions
	double precision function f2dimd1(xx)
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
	common/baspring/ispring
	common/nboucle/iloop
	common/aurenche/iauren
	common/processd/j_processd_min,j_processd_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	dimension xx(20)
	parameter (k0max=4)
	dimension sv2(k0max),sy3(k0max)
	dimension temp0(k0max),temp1(k0max)
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
c ------------------ gpt3 ------------------------------------------
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3))
     #	,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
c ------------------ u -------------------------------------------
	u1 = xx(4)
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)
	call fspcoud(n,spcou)
	gpt4 = gpt3
	x10 = gpt3*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = gpt3*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimd1 = 0.d0
	  return
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
	    f2dimd1 = 0.d0
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
	alpha_em = alphaem(iloopem,m*m)
c c^b_{ij} = 2*pi*alpha^2/(4*c_i*c_j*s**2)
	cbij = alpha_em**2*2.d0*pi/(s*s)
	alpi = alfas(iloop,mu*mu)/(2.d0*pi)
c
	call ssy3(s,gpt3,y3,y4,x10,x20,m,u1,ih1,ih2,sy3)
	call ssv2(s,gpt3,y3,y4,x10,x20,u1,ih1,ih2,sv2)
	call vec_dadd_vector(sy3,sv2,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dimd1 = 0.d0
	do i = j_processd_min,j_processd_max
	  f2dimd1 = temp1(i) + temp1(2+i) + f2dimd1
	enddo
ccc	f2dimd1 = vec_dsum(temp1,k0max)
c
	f2dimd1 = f2dimd1 * xjacob * alpi * cbij * hc2/2.d0
c
	if (ispring.eq.1) then
	  f2dimd1 = dabs(f2dimd1)
	else if (ispring.eq.2) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  resfunc = f2dimd1
	  f2dimd1 = dabs(f2dimd1)
	endif
c
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimd2(xx)
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
	common/baspring/ispring
	common/nboucle/iloop
	common/aurenche/iauren
	common/processd/j_processd_min,j_processd_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=4)
	dimension sy3r(k0max),sy3p(k0max)
	dimension temp0(k0max),temp1(k0max)
	dimension spcou(k0max)
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
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	pt3 = gpt3
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	x3min_incl = 2.d0*gpt3/dsqrt(s)*dcosh(y3)
	if (i_flag_iso.eq.1) then
	  x3min_isol = gpt3/(gpt3+etmax)
	else if (i_flag_iso.eq.2) then
	  x3min_isol = 1.d0/(1.d0+etmax)
	endif
	x3min = dmax1(x3min_incl,x3min_isol)
	x10 = pt3*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = pt3*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimd2 =0.d0
	  return
	endif
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
	z3max0 = 1.d0
	z3min0 = dmax1(x3min,gpt3/gpt4maxp,x10,x20)
	z3min1 = pt3/(ptm_exp+pt3)	
	z3sup = pt3/(ptm+pt3)
	un = 1.d0
	if (icoup_exp.eq.0) then
	  z3min = z3min0	
	  z3max = z3max0	
	else if (icoup_exp.eq.1) then
	  z3min = dmax1(z3min0,z3min1)	
	  z3max = z3sup	
	endif
c juste pour tester (en principe inutile)
	if (z3min.gt.z3max) then
	  f2dimd2 = 0.d0
	  write (18,*) 'attention z3min > z3max',z3min,z3max
	  return
	endif
	z3 = z3min + (z3max-z3min)*xx(4)
c
	gpt4 = gpt3/z3
c juste pour tester (en principe inutile)
	if (gpt4.gt.gpt4max.or.gpt4.lt.gpt4min) then
	  f2dimd2 = 0.d0
	  write (18,*) 'attention gpt4 > gpt4min ou gpt4 < gpt4max',
     #	              gpt4,gpt4min,gpt4max
	  return
	endif
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)*
     #  (z3max-z3min)	
	call fspcoud(n,spcou)
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
	    f2dimd2 = 0.d0
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
	alpha_em = alphaem(iloopem,m*m)
c c^b_{ij} = 2*pi*alpha^2/(4*c_i*c_j*s**2)
	cbij = alpha_em**2*2.d0*pi/(s*s)
	alpi = alfas(iloop,mu*mu)/(2.d0*pi)
c
	call ssy3p(s,gpt3,x10,x20,z3,y3,y4,mf,ih1,ih2,sy3p)
	if (z3.le.z3sup) then
	  call ssy3r(s,gpt3,x10,x20,z3,r,y3,y4,ih1,ih2,sy3r)
	else
	  call vec_dinit(sy3r,k0max,zero)	  
	endif
	call vec_dadd_vector(sy3r,sy3p,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dimd2 = 0.d0
	do i = j_processd_min,j_processd_max
	  f2dimd2 = temp1(i) + temp1(2+i) + f2dimd2
	enddo
ccc	f2dimd2 = vec_dsum(temp1,k0max)
	f2dimd2 = f2dimd2 * xjacob * alpi * cbij * hc2/2.d0
c
	if (ispring.eq.1) then
	  f2dimd2 = dabs(f2dimd2)
	else if (ispring.eq.2) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  resfunc = f2dimd2
	  f2dimd2 = dabs(f2dimd2)
	endif
c
	return
	end 
c integrale a 2 dimensions
	double precision function f2dimd3(xx)
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
	common/baspring/ispring
	common/nboucle/iloop
	common/aurenche/iauren
	common/processd/j_processd_min,j_processd_max
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=4)
	dimension sy4r(k0max),sy4p(k0max)
	dimension temp0(k0max),temp1(k0max)
	dimension spcou(k0max)
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
c ------------------ gpt4 ------------------------------------------
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4))
     #	,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt4max = dmax1(gptinf,gpt4maxp)
	gpt4min = gptinf
	gpt4 = gpt4min + (gpt4max-gpt4min)*xx(3)
	gpt3maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt3min = gptinf
	pt4 = gpt4
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	x4min_incl = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	if (i_flag_iso.eq.1) then
	  x4min_isol = gpt4/(gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  x4min_isol = 1.d0/(1.d0+etmax)
	endif
	x4min = dmax1(x4min_incl,x4min_isol)
	x10 = pt4*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = pt4*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f2dimd3 = 0.d0
	  return
	endif
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c
	z4max0 = 1.d0
	z4min0 = dmax1(x4min,gpt4/gpt3maxp,x10,x20)
	z4min1 = pt4/(ptm_exp+pt4)	
	z4sup = pt4/(ptm+pt4)
	if (icoup_exp.eq.0) then
	  z4min = z4min0
	  z4max = z4max0
	else if (icoup_exp.eq.1) then
	  z4min = dmax1(z4min0,z4min1)
	  z4max = z4sup
	endif
c juste pour tester (en principe inutile)
	if (z4min.gt.z4max) then
	  f2dimd3 = 0.d0
	  write (18,*) 'attention z4min > z4max',z4min,z4max
	  return
	endif
	z4 = z4min + (z4max-z4min)*xx(4)
c
	gpt3 = gpt4/z4
c juste pour tester (en principe inutile)
	if (gpt3.gt.gpt3max.or.gpt3.lt.gpt3min) then
	  f2dimd3 = 0.d0
	  write (18,*) 'attention gpt3 > gpt3min ou gpt3 < gpt3max',
     #	              gpt3,gpt3min,gpt3max
	  return
	endif
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt4max-gpt4min)
     #  *(z4max-z4min)
	call fspcoud(n,spcou)
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
	    f2dimd3 = 0.d0
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
	alpha_em = alphaem(iloopem,m*m)
c c^b_{ij} = 2*pi*alpha^2/(4*c_i*c_j*s**2)
	cbij = alpha_em**2*2.d0*pi/(s*s)
	alpi = alfas(iloop,mu*mu)/(2.d0*pi)
c
	call ssy4p(s,gpt4,x10,x20,z4,y3,y4,mf,ih1,ih2,sy4p)
	if (z4.le.z4sup) then
	  call ssy4r(s,gpt4,x10,x20,z4,r,y3,y4,ih1,ih2,sy4r)
	else
	  call vec_dinit(sy4r,k0max,zero)	  
	endif
	call vec_dadd_vector(sy4r,sy4p,k0max,temp0)
	call vec_dmult_vector(temp0,spcou,k0max,temp1)
	f2dimd3 = 0.d0
	do i = j_processd_min,j_processd_max
	  f2dimd3 = temp1(i) + temp1(2+i) + f2dimd3
	enddo
ccc	f2dimd3 = vec_dsum(temp1,k0max)
c
	f2dimd3 = f2dimd3 * xjacob * alpi * cbij * hc2/2.d0
	if (ispring.eq.1) then
	  f2dimd3 = dabs(f2dimd3)
	else if (ispring.eq.2) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  resfunc = f2dimd3
	  f2dimd3 = dabs(f2dimd3)
	endif
c
	return
	end 
c*****************************************************
c integrale a 1 dimension
	double precision function f1dimd(xx)
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
	common/baspring/ispring
	common/nboucle/iloop
	common/aurenche/iauren
	common/processd/j_processd_min,j_processd_max
	common/box/ibox
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/coup_histo/r_isol,etmax
	common/flag_iso/i_flag_iso
	dimension xx(20)
	parameter (k0max=4)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	dimension c0sv12(k0max),c0sv3(k0max),c0sv4(k0max),c0fborn(k0max)
	dimension spcou(k0max)
	dimension p3(4),p4(4),p5(4)
	pi = datan(1.d0)*4.d0
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
	gpt3sup = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y3)))
	gpt4sup = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
	gpt4maxp = dmin1(gptsup,dsqrt(s)/(2.d0*dcosh(y4)))
c	gpt3max = dmax1(gptinf,gpt3maxp)
	gpt4max = dmax1(gptinf,gpt4maxp)
c	gpt3min = gptinf
	gpt4min = gptinf
	gpt3 = gpt3min + (gpt3max-gpt3min)*xx(3)
	gpt4 = gpt3
c
	x10 = gpt3*(dexp(y3)+dexp(y4))/dsqrt(s)
	x20 = gpt3*(dexp(-y3)+dexp(-y4))/dsqrt(s)
	if (x10.gt.un.or.x20.gt.un) then
	  f1dimd = 0.d0
	  return
	endif
	xjacob = (y3max-y3min)*(y4max-y4min)*(gpt3max-gpt3min)
	call fspcoud(n,spcou)
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
	    f1dimd = 0.d0
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
	alpha_em = alphaem(iloopem,m*m)
c c^b_{ij} = 2*pi*alpha^2/(4*c_i*c_j*s**2)
	cbij = alpha_em**2*2.d0*pi/(s*s)
	alphas = alfas(iloop,mu*mu)
	alpi = alphas/(2.d0*pi)
c
	if (iborn.eq.0) then 
	  if (i_flag_iso.eq.1) then
	    x3min_isol = gpt3/(gpt3+etmax)
	    x4min_isol = gpt4/(gpt4+etmax)
	  else if (i_flag_iso.eq.2) then
	    x3min_isol = 1.d0/(1.d0+etmax)
	    x4min_isol = 1.d0/(1.d0+etmax)
	  endif
	  call ssv12(s,gpt3,gpt3,y3,y4,x10,x20,m,mu,ih1,ih2,c0sv12)
	  call ssv3(s,gpt3,gpt4sup,x3min_isol,y3,y4,x10,x20,mf,ih1,ih2,
     #              c0sv3)
	  call ssv4(s,gpt4,gpt3sup,x4min_isol,y3,y4,x10,x20,mf,ih1,ih2,
     #              c0sv4)
	  call vec_dmult_vector(c0sv12,spcou,k0max,temp0)
	  call vec_dmult_vector(c0sv3,spcou,k0max,temp1)
	  call vec_dmult_vector(c0sv4,spcou,k0max,temp2)
	  sv12 = 0.d0
	  sv3 = 0.d0
	  sv4 = 0.d0
	  do i = j_processd_min,j_processd_max
	    sv12 = temp0(i) + temp0(2+i) + sv12 
	    sv3 = temp1(i) + temp1(2+i) + sv3
	    sv4 = temp2(i) + temp2(2+i) + sv4
	  enddo
cc	  sv12 = vec_dsum(temp0,k0max)
cc	  sv3 = vec_dsum(temp1,k0max)
cc	  sv4 = vec_dsum(temp2,k0max)
	  born = 0.d0
	else if (iborn.eq.1) then
	  box = 0.d0
	  fborn = 0.d0
	  if (ibox.eq.0.or.ibox.eq.2) then
	    call sfborn(s,n,gpt3,y3,y4,ih1,ih2,c0fborn)
	    call vec_dmult_vector(c0fborn,spcou,k0max,temp3)
	    do i = j_processd_min,j_processd_max
	      fborn = temp3(i) + temp3(2+i) +fborn
	    enddo
cc	    fborn = vec_dsum(temp3,k0max)
	  endif
	  if (ibox.eq.1.or.ibox.eq.2) then
	    call glgl_gaga(s,n,gpt3,y3,y4,ih1,ih2,boite)
	    box =  boite*alphas**2/(16.d0*pi**2)
	  endif
	  sv12 = 0.d0
	  sv3 = 0.d0
	  sv4 = 0.d0
	  born = fborn + box
	endif
c 
	resb = born + (sv12+sv3+sv4)*alpi
	f1dimd = resb * xjacob * cbij * hc2/2.d0
c
	if (ispring.eq.1) then
	  f1dimd = dabs(f1dimd)
	else if (ispring.eq.2) then
	  xfi12 = pi
	  xfi13 = 0.d0
	  xpt1 = gpt3
	  xy1 = y3
	  xpt2 = gpt4
	  xy2 = y4
	  xpt3 = 0.d0
	  xy3 = 0.d0
	  resfunc = f1dimd
	  f1dimd = dabs(f1dimd)
	endif
	return
	end 
c********************************************************************
c===================================================================
	subroutine ampd3(s,x1,x2,y3,pt3,y5,pt5,fi35,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=4)
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
	call sh12d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h12)
	call sh13d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h13)
	call sh23d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h23)
	call sh34d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	call sconsd(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,cons)
	do i=1,k0max
	  camp(i) = (0.5d0*h12(i)*e12+h13(i)*e13
     #	             +h23(i)*e23+h34(i)*e34p+0.5d0*cons(i))
	enddo
	return
	end
c===================================================================
	subroutine ampd_corr3(s,x1,x2,y3,pt3,pt5,rs,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=4)
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
	call sh13d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h13)
	call sh23d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h23)
	call sh34d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	do i=1,k0max
	  camp(i) = (h13(i)+h23(i)+h34(i))*2.d0/rs/(pt5*pt5)
	enddo
	return
	end
c===================================================================
	subroutine ampd4(s,x1,x2,y4,pt4,y5,pt5,fi45,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=4)
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
	call sh12d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h12)
	call sh14d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h14)
	call sh24d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h24)
	call sh34d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	call sconsd(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,cons)
	do i=1,k0max
	  camp(i) = (0.5d0*h12(i)*e12+h14(i)*e14
     #	             +h24(i)*e24+h34(i)*e34s+0.5d0*cons(i))
	enddo
	return
	end
c===================================================================
	subroutine ampd_corr4(s,x1,x2,y4,pt4,pt5,rs,camp)
	implicit real*8(a-h,l-z)
	parameter (k0max=4)
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
	call sh14d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h14)
	call sh24d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h24)
	call sh34d(s12,s13,s14,s15,s23,s24,s25,s34,s35,s45,h34)
	do i=1,k0max
	  camp(i) = (h14(i)+h24(i)+h34(i))*2.d0/rs/(pt5*pt5)
	enddo
	return
	end
c********************************************************************
cc	subroutine trapnan(i,x)
cc	implicit real*8(a-h,k-z)
cc	zero = 0.d0
cc	if (.not.(x.gt.zero).and..not.(x.le.zero)) then
cc	  write (19,*) 'nan',i,x
cc	  call flu(19)
cc	endif
cc	return
cc	end
c********************************************************************
cdebug	subroutine xcount
cdebug	common/nbcoup/icutte
cdebug	icutte= icutte+1
cdebug	return
cdebug	end
c********************************************************************
