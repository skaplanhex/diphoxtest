	subroutine histo_init(inbevent)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	character*128 path_rzfile
	parameter(iwpawc = 4200000)
	common/pawc/hmemor(iwpawc)
	common/cheminrz/path_rzfile
	common/longrz/ilenrz
	parameter (inbhisto_equi=4)
	common/init_close_equi/ibin(inbhisto_equi),wxmin(inbhisto_equi),
     #	wxmax(inbhisto_equi)
	iwpawt = 150000*7*inbevent/1000000
	if(iwpawt.gt.iwpawc) then
	  write(6,*)' you must have iwpawc at least ',iwpawt
	stop
	endif
	call hlimit(iwpawc)
	call hropen(1,'fixed_order',path_rzfile(1:ilenrz)//'.dat'
     #	,'n',1024,istat)
	ibin(1) = 40
	wxmin(1) = 25.
	wxmax(1) = 185.
c
	ibin(2) = 40
	wxmin(2) = 25.
	wxmax(2) = 185.
c
	ibin(3) = 20
	wxmin(3) = 0.
	wxmax(3) = 100.
c
	ibin(4) = 20
	wxmin(4) = 0.
	wxmax(4) = 100.
c
	call hbook1(20,'dsigma_lo/dpt_gamma',ibin(1),wxmin(1),wxmax(1),0.)
	call hbook1(920,'dummy',ibin(1),wxmin(1),wxmax(1),0.)
	call hbarx(20)
	call hbarx(920)
c
	call hbook1(21,'dsigma_nlo/dpt_gamma',ibin(2),wxmin(2),wxmax(2),0.)
	call hbook1(921,'dummy',ibin(2),wxmin(2),wxmax(2),0.)
	call hbarx(21)
	call hbarx(921)
c
	call hbook1(22,'dsigma_lo/dqt',ibin(3),wxmin(3),wxmax(3),0.)
	call hbook1(922,'dummy',ibin(3),wxmin(3),wxmax(3),0.)
	call hbarx(22)
	call hbarx(922)
c
	call hbook1(23,'dsigma_nlo/dqt',ibin(4),wxmin(4),wxmax(4),0.)
	call hbook1(923,'dummy',ibin(4),wxmin(4),wxmax(4),0.)
	call hbarx(23)
	call hbarx(923)
c
	return
	end
