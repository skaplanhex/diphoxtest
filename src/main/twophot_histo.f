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
	integer*4 nbhisto
        parameter (maxtrk = 3)
	common/npt/inpt
	common/fixle/iprov,weight
	common/sortier/ntrack,we(maxtrk),wpx(maxtrk),
     #	wpy(maxtrk),wpz(maxtrk)
	common/sfragvr/wx3,wx4
	common/approx/lo,nlo
	common/calcul/integ,intega,gener
	common/nbevent/inbevent
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
c	
c	
	call param
c###
	open(unit=8,file=path_bsfile(1:ilen)//'sigmaint.res',
     #	status='unknown')
c
c******************************************************************
c###
	if (gener) then
	  call histo_init(inbevent)
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
c normalisation
	  open(unit=30,file=path_rzfile(1:ilenrz)//'.res',
     #	  status='unknown')
	  read(30,100) inumb_event,wfirst_int
100       format(1x,i12,e12.5)
** 	write (*,*) 'coucou1',inumb_event,wfirst_int
c
	  call histo_close(inbevent,wfirst_int)
	endif  
c###
	close (8)
	end
