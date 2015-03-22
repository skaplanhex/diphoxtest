c
	program main
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	logical lo,nlo,gener,integ,intega
	logical dir,onef,twof
	logical distrib,resom
	logical pt_pair,pt_photon
	character*128 path_bsfile,path_rzfile
        integer*4 ntrack,iprov,maxtrk
        parameter (maxtrk = 3)
	common/npt/inpt
	common/fixle/iprov,weight
	common/sortier/ntrack,we(maxtrk),wpx(maxtrk),
     #	wpy(maxtrk),wpz(maxtrk)
	common/sfragvr/wx3,wx4
	common/approx/lo,nlo
	common/calcul/integ,intega,gener
	common/nbevent/inbevent
        parameter(iwpawc = 3000000)
        common/pawc/hmemor(iwpawc)
	common/quest/iquest(100)
	common/pclass/dir,onef,twof
	common/cheminbs/path_bsfile
	common/long/ilen
	common/cheminrz/path_rzfile
	common/longrz/ilenrz
	common/calchoi/distrib,resom
	common/diffpt/pt_pair,pt_photon
	common/qtpara1/ibinqt
	common/qtpara2/wbinqt(100)
	common/ptpara1/ibinpt
	common/ptpara2/wbinpt(100)
	character*128 name_experiment
	common/name/name_experiment
	common/lname/ilname
	call param
c###
	if (gener) then
        iwpawt = 150000*7*inbevent/1000000
        if(iwpawt.gt.iwpawc) then
          write(6,*)' you must have iwpawc at least ',iwpawt 
          stop  
        endif
	call hlimit(iwpawc)
	endif
c###
c###
	if (distrib.or.resom) then
	  call hlimit(iwpawc)
	endif
c###
	open(unit=8,file=path_bsfile(1:ilen)//'sigmaint.res',
     #	status='unknown')
c
c******************************************************************
c###
	if (gener) then
	iquest(10) = 65000
	call hcdir('//pawc',' ')
	call hropen(50,'ggd',
     #	path_rzfile(1:ilenrz)//'.rz',
     #	'nq',8192,istat)
	if (istat.ne.0) then
	  write (8,*) 'error in hropen',istat
	  stop
	endif
	call hbnt(10,'gamma-gamma events','d')
        call hldir('//pawc','t')
        call hldir('//gg','t')
	call hbname(10,'fixblo',iprov,'iprov, weight')
	call hbname(10,'varblok2',ntrack,'ntrack[1,3],'//
     +            'we(ntrack), wpx(ntrack), wpy(ntrack), wpz(ntrack)')
     	if (onef.or.twof) then
	  call hbname(10,'fragblo',wx3,'x3, x4')
	endif
	call hbset('bsize',8192,ierr)
	if (ierr.ne.0) then
	  write (8,*) 'error in hbset',ierr
	  stop
	endif
	
	write(6,*)' before rzldir deb'	
	call rzldir(' ',' ')
	endif	
c###
c###
	if (distrib.or.resom) then
	  call hropen(1,'resommation',
     #	  path_rzfile(1:ilenrz)//'.dat','n',1024,istat)
	  if (name_experiment(1:ilname).eq.'d0') then
c nb. de bins et bornes pour les histogrammes en pt de la paire
	    ibinqt = 14
	    wbinqt(1) = 0.
	    wbinqt(2) = 2.
	    wbinqt(3) = 4.
	    wbinqt(4) = 6.
	    wbinqt(5) = 8.
	    wbinqt(6) = 10.
	    wbinqt(7) = 12.
	    wbinqt(8) = 14.
	    wbinqt(9) = 16.
	    wbinqt(10) = 20.
	    wbinqt(11) = 24.
	    wbinqt(12) = 32.
	    wbinqt(13) = 40.
	    wbinqt(14) = 56.
	    wbinqt(15) = 80.
c nb. de bins et bornes pour les histogrammes en pt
c nb. de bins et bornes pour les histogrammes en pt
	    ibinpt = 9
	    wbinpt(1) = 16.
	    wbinpt(2) = 18.
	    wbinpt(3) = 20.
	    wbinpt(4) = 22.
	    wbinpt(5) = 24.
	    wbinpt(6) = 28.
	    wbinpt(7) = 32.
	    wbinpt(8) = 40.
	    wbinpt(9) = 56.
	    wbinpt(10) = 80.
	  elseif (name_experiment(1:ilname).eq.'lhc') then
c nb. de bins et bornes pour les histogrammes en pt de la paire
	    ibinqt = 20
	    qtmin = 0.d0
	    qtmax = 100.d0
	    qtlarg = (qtmax-qtmin)/dfloat(ibinqt)
	    do i=1,ibinqt+1
	      binqt = dfloat(i-1)*qtlarg+qtmin
	      wbinqt(i) = sngl(binqt)
	    enddo
c nb. de bins et bornes pour les histogrammes en pt
	    ibinpt = 10
	    ptpmin = 40.d0
	    ptpmax = 140.d0
	    ptplarg = (ptpmax-ptpmin)/dfloat(ibinpt)
	    do i=1,ibinpt+1
	      binpt = dfloat(i-1)*ptplarg+ptpmin
	      wbinpt(i) = sngl(binpt)
	    enddo
	  elseif (name_experiment(1:ilname).eq.'wa70') then
c nb. de bins et bornes pour les histogrammes en qt
	    ibinqt = 12
	    qtmin = 0.d0
	    qtmax = 3.84d0
	    qtlarg = (qtmax-qtmin)/dfloat(ibinqt)
	    do i=1,ibinqt+1
	      binqt = dfloat(i-1)*qtlarg+qtmin
	      wbinqt(i) = sngl(binqt)
	    enddo
c nb. de bins et bornes pour les histogrammes en pt
	    ibinpt = 8
	    ptpmin = 3.d0
	    ptpmax = 7.d0
	    ptplarg = (ptpmax-ptpmin)/dfloat(ibinpt)
	    do i=1,ibinpt+1
	      binpt = dfloat(i-1)*ptplarg+ptpmin
	      wbinpt(i) = sngl(binpt)
	    enddo
	  endif
	  if (distrib) then
	    if (pt_pair) then
	      call hbookb(91,'dsigmares/dqt',ibinqt,wbinqt,0.)
	    endif
	    if (pt_photon) then
	      call hbookb(81,'dsigmares/dpt',ibinpt,wbinpt,0.)
	    endif
	  endif
	  if (resom) then
	    if (pt_pair) then
	      call hbookb(92,'dsigmares/dqt',ibinqt,wbinqt,0.)
	    endif
	    if (pt_photon) then
	      call hbookb(82,'dsigmaress/dpt',ibinpt,wbinpt,0.)
	    endif
	  endif
	endif
c###
c******************************************************************
c        
	open(unit=12,
     #	file=path_rzfile(1:ilenrz)//'.res'
     #	,status='unknown')
cccc
	if ((dir.and.onef).or.(dir.and.twof).or.(onef.and.twof)) then
	  write(6,*) 'stop...too many options selected'
	  write(8,*) 'stop...too many options selected'
	  stop
	endif
cccc
	if (dir) then
	  call dir_sub
	endif
	if (onef) then
	  call onef_sub
	endif
	if (twof) then
	  call twof_sub
	endif
cccc	
        close(unit=12)
c	
c
c###
	if (gener) then
	write(6,*)' before hcdir //gg'
	call hcdir('//ggd',' ')
	call hldir('//ggd','t')
	write(6,*)' before rzldir end'
	call rzldir(' ',' ')
	write(6,*)' before hrout'
	call hrout(10,icycle,' ')
        call hldir(' ',' ')
        call hldir(' ','a')
	write(6,*)' before hrend'
	call hrend('ggd')
	endif  
c###
	if (distrib) then
	  if (pt_pair) then
	    call hrout(91,icycle,' ')
	  endif
	  if (pt_photon) then
	    call hrout(81,icycle,' ')
	  endif
	endif
	if (resom) then
	  if (pt_pair) then
	    call hrout(92,icycle,' ')
	  endif
	  if (pt_photon) then
	    call hrout(82,icycle,' ')
	  endif
	endif
	if (distrib.or.resom) then
	  call hrend('resommation')
	endif
	close (8)
	end
