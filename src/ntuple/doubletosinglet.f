	subroutine doubletosinglet(rnumb1,pi,poid,iiprov,intrack)
	implicit real*8 (a-h,l-v,x-z)
        implicit real*4 (w)
        integer*4 ntrack,iprov,maxtrk
        parameter (maxtrk = 3)
	common/sortie/xpt1,xy1,xpt2,xy2,xfi12,xpt3,xy3,xfi13,resfunc
	common/fixle/iprov,weight
        common/sortier/ntrack,we(maxtrk),wpx(maxtrk),
     #	wpy(maxtrk),wpz(maxtrk)
	common/sfragvr/wx3,wx4
	common/sfragv/xfrag1,xfrag2
	common/efficacite/ipositif,inegatif
	iprov = iiprov
	ntrack = intrack
	wzero = 0.
	dpi = 2.d0*pi     
	xfi = dpi*rnumb1
	xfi1 = xfi
	xfi2 = dmod(xfi12+xfi,dpi)
	xfi3 = dmod(xfi13+xfi,dpi)
	do i = 1,ntrack
	  if(i.eq.1) then
	    if (iiprov.eq.44) then
	      we(i) = sngl(xpt2*dcosh(xy2))
	      wpx(i) = sngl(xpt2*dcos(xfi2))
	      wpy(i) = sngl(xpt2*dsin(xfi2))
	      wpz(i) = sngl(xpt2*dsinh(xy2))
	    else
	      we(i) = sngl(xpt1*dcosh(xy1))
	      wpx(i) = sngl(xpt1*dcos(xfi1))
	      wpy(i) = sngl(xpt1*dsin(xfi1))
	      wpz(i) = sngl(xpt1*dsinh(xy1))
	    endif
	  elseif(i.eq.2)then
	    if (iiprov.eq.44) then
	      we(i) = sngl(xpt1*dcosh(xy1))
	      wpx(i) = sngl(xpt1*dcos(xfi1))
	      wpy(i) = sngl(xpt1*dsin(xfi1))
	      wpz(i) = sngl(xpt1*dsinh(xy1))
	    else
	      we(i) = sngl(xpt2*dcosh(xy2))
	      wpx(i) = sngl(xpt2*dcos(xfi2))
	      wpy(i) = sngl(xpt2*dsin(xfi2))
	      wpz(i) = sngl(xpt2*dsinh(xy2))
	    endif
	  elseif(i.eq.3)then
	    we(i) = sngl(xpt3*dcosh(xy3))
	    wpx(i) = sngl(xpt3*dcos(xfi3))
	    wpy(i) = sngl(xpt3*dsin(xfi3))
	    wpz(i) = sngl(xpt3*dsinh(xy3))
	  endif	  
	enddo
	if (iiprov.eq.44) then
	  wx3 = sngl(xfrag2)
	  wx4 = sngl(xfrag1)
	else
	  wx3 = sngl(xfrag1)
	  wx4 = sngl(xfrag2)
	endif
	weight = sngl( dsign(poid,resfunc) )
	if (weight.lt.wzero) then
	  inegatif = inegatif+1
	else
	  ipositif = ipositif+1
	endif
	return
	end

c********************************************************************
