	subroutine histo_close(inbevent,wfirst_int)
	implicit real*8 (a-h,l-v,x-z)
	implicit real*4 (w)
	parameter (inbhisto_equi=4)
	common/init_close_equi/ibin(inbhisto_equi),wxmin(inbhisto_equi),
     #	wxmax(inbhisto_equi)
	dimension wbin_size(inbhisto_equi),wnorma(inbhisto_equi)
	do i=1,inbhisto_equi
	  wbin_size(i) = (wxmax(i)-wxmin(i))/float(ibin(i))
	  wnorma(i) = wfirst_int/(float(inbevent)*wbin_size(i))
	enddo
	call hopera(20,'+',920,20,wnorma(1),1.)
	call hopera(21,'+',921,21,wnorma(2),1.)
	call hopera(22,'+',922,22,wnorma(3),1.)
	call hopera(23,'+',923,23,wnorma(4),1.)
	call hrout(20,ICYCLE,' ')
	call hrout(21,ICYCLE,' ')
	call hrout(22,ICYCLE,' ')
	call hrout(23,ICYCLE,' ')
	call hrend('fixed_order')
	return
	end
