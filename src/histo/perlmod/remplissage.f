	subroutine remplissage(iprov,xmgg,pt3,pt4,qt,fi34,
     #	y3,y4,yboost,ystar,ygg,p_out3,p_out4,z3_trig,z4_trig,
     #	cos_thetas,wx3,wx4,weight)
	implicit real*4 (a-h,l-z)
	  if (iprov.eq.11) then
	  call hfill (20,xmgg,0.,weight)
	  call hfill (22,qt,0.,weight)
	  call hfill (24,fi34,0.,weight)
	  endif
	  call hfill (21,xmgg,0.,weight)
	  call hfill (23,qt,0.,weight)
	  call hfill (25,fi34,0.,weight)
	return
	end
