c********************************************************************
	real function wacos(x)
	implicit real*4 (a-h,l-z)
	un = 1.
	unp = 1. + 1.e-4
	zero = 0.
	pi=atan(1.)*4.
	if (abs(x).le.un) then
	  wacos = acos(x)
	else if (abs(x).gt.un.and.abs(x).lt.unp) then
	  if (x.lt.zero) then
	    wacos = pi
	  else
	    wacos = 0.
	  endif
	else
	  write (*,*) 'alerte dans wacos argument > 1:',x
	  wacos = 0.
	endif
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
C correction d'un bug le 29/07/99, pt5 n'etait pas definie pour le twof
C correction d'un bug le 12/05/00, mauvaise definition de l'energie depose 
c pour le onef
	subroutine gfill
	implicit real*8 (a-h,l-v,x-z)
        implicit real*4 (w)
	logical ident
	integer maxtrk,ntrack
        parameter (maxtrk = 3)
	logical rapidite,pete,isol,masscut
	logical dir,onef,twof
        common/sortier/xe(maxtrk),xpx(maxtrk),
     #	xpy(maxtrk),xpz(maxtrk)
	common/fixle/iprov,weight,ntrack
	common/pclass/dir,onef,twof
	common/sfragvr/xx3,xx4
** 	call further_param(pt1min,pt2min,ymax,ymin,etmin,
**      #	r,xmasmax,xmasmin)
	call further_param(pt1min,pt2min,y1max,y1min,
     #	y2max,y2min,xmasmax,xmasmin,ident)
c dans cette version on ne permet pas de mettre un critere d'isolement
c ici
	r = 0.4
	etmin = 7000.
	isol = .true.
	xinfty = 1.d+8
	zero = 0.d0
	un = 1.d0	
	if (dir) then
	  xx3 = 1.
	  xx4 = 1.
	endif 
	xmgg2 = 2.*
     #	(xe(1)*xe(2)-xpx(1)*xpx(2)-xpy(1)*xpy(2)-xpz(1)*xpz(2))
	xmgg = sqrt(xmgg2)
	pt3 = sqrt(xpx(1)**2+xpy(1)**2)
	pt4 = sqrt(xpx(2)**2+xpy(2)**2)
	y3 = log((xe(1)+xpz(1))/(xe(1)-xpz(1)))*0.5
	y4 = log((xe(2)+xpz(2))/(xe(2)-xpz(2)))*0.5
	ystar = (y3-y4)/2.
	yboost = (y3+y4)/2.
	cos_thetas = abs(tanh(ystar))
c ygg est la rapidite du systeme photon-photon (pion-pion)
	ygg = log((xe(1)+xe(2)+xpz(1)+xpz(2))/
     #	(xe(1)+xe(2)-xpz(1)-xpz(2)))*0.5
	if (ntrack.eq.3) then
	  pt5 = sqrt(xpx(3)**2+xpy(3)**2)
	  xnumy5 = (xe(3)+xpz(3))
	  xdeny5 = (xe(3)-xpz(3))
	  if (xnumy5.eq.zero.or.xdeny5.eq.zero) then
	    if (xpz(3).le.zero) then
	      y5 = -xinfty
	    else
	      y5 = xinfty
	    endif
	  else if (xnumy5.gt.zero.and.xdeny5.gt.zero) then
	    y5 = log((xe(3)+xpz(3))/(xe(3)-xpz(3)))*0.5
	  endif
	else if (iprov.eq.22.or.iprov.eq.32) then
	  y5 = y3
	  if (dir) then
	    pt5 = pt4-pt3
	  else if (onef) then
	    if (xx4.eq.1.) then
	      pt5 = pt4-pt3/xx3
	    else if (xx3.eq.1.) then
	      pt5 = pt3-pt4/xx4
	    endif
	  else if (twof) then
	    pt5 = pt4/xx4-pt3/xx3
	  endif
	else if (iprov.eq.23.or.iprov.eq.33) then
	  y5 = y4
	  if (dir) then
	    pt5 = pt3-pt4
	  else if (onef) then
	    if (xx4.eq.un) then
	      pt5 = pt3/xx3-pt4
	    else if (xx3.eq.un) then
	      pt5 = pt4/xx4-pt3
	    endif
	  else if (twof) then
	    pt5 = pt3/xx3-pt4/xx4
	  endif
	else
	  pt5 =0.
	  y5 = 0.
	endif
	cfi34 = (xpx(1)*xpx(2)+xpy(1)*xpy(2))/(pt3*pt4)
	if (ntrack.eq.3) then
	  cfi35 = (xpx(1)*xpx(3)+xpy(1)*xpy(3))/(pt3*pt5)
	  cfi45 = (xpx(3)*xpx(2)+xpy(3)*xpy(2))/(pt5*pt4)
	  fi35 = zacos(cfi35)
	  fi45 = zacos(cfi45)
	endif
	pi = atan(un)*4.
	if (abs(cfi34).gt.1.0001) then
	  write (6,*) 'alerte',cfi34
	else if (abs(cfi34).gt.un) then
	  cfi34 = sign(un,cfi34)
	else
	  cfi34 = cfi34
	endif
	fi34 = zacos( cfi34 )
	qt = sqrt( (xpx(1)+xpx(2))**2+(xpy(1)+xpy(2))**2 )
	p_out3 = (xpx(1)*xpy(2)-xpx(2)*xpy(1))/pt3
	p_out4 = (xpx(2)*xpy(1)-xpx(1)*xpy(2))/pt4
	z3_trig = -pt4/pt3*cfi34
	z4_trig = -pt3/pt4*cfi34
c
	if (ntrack.eq.3) then
	  dist35 = sqrt((y3-y5)**2+(fi35)**2)
	  dist45 = sqrt((y4-y5)**2+(fi45)**2)
	  if ((dist35.gt.r)
     #	  .and.(dist45.gt.r)) then
	    edep3 = (1.-xx3)/xx3*pt3
	    edep4 = (1.-xx4)/xx4*pt4
** 	    isol = (edep3.le.etmin).and.(edep4.le.etmin)
	  else if ((dist35.le.r)
     #	  .and.(dist45.gt.r)) then
	    edep3 = pt5 + (1.-xx3)/xx3*pt3
	    edep4 = (1.-xx4)/xx4*pt4
** 	    isol = (edep3.le.etmin).and.(edep4.le.etmin)
	  else if ((dist35.gt.r)
     #	  .and.(dist45.le.r)) then
	    edep3 = (1.-xx3)/xx3*pt3
	    edep4 = pt5 + (1.-xx4)/xx4*pt4
** 	    isol = (edep3.le.etmin).and.(edep4.le.etmin)
	  endif
	else if (iprov.eq.22.or.iprov.eq.32) then
	  if (dir) then
	    edep3 = pt5 
	    edep4 = 0.
	  else if (onef) then
	    if (xx4.eq.un) then
	      edep3 = pt5 + (1.-xx3)/xx3*pt3
	      edep4 = 0.
	    else if (xx3.eq.un) then
	      edep4 = pt5 + (1.-xx4)/xx4*pt4
	      edep3 = 0.
	    endif
	  else if (twof) then
	    edep3 = pt5 + (1.-xx3)/xx3*pt3
	    edep4 = (1.-xx4)/xx4*pt4
	  endif
** 	  isol = (edep3.le.etmin).and.(edep4.le.etmin)
	else if (iprov.eq.23.or.iprov.eq.33) then
	  if (dir) then
	    edep4 = pt5 
	    edep3 = 0.
	  else if (onef) then
	    if (xx4.eq.un) then
	      edep4 = pt5
	      edep3 = (1.-xx3)/xx3*pt3
	    else if (xx3.eq.un) then
	      edep3 = pt5
	      edep4 = (1.-xx4)/xx4*pt4
	    endif
	  else if (twof) then
	    edep4 = pt5 + (1.-xx4)/xx4*pt4
	    edep3 = (1.-xx3)/xx3*pt3
	  endif
** 	  isol = (edep3.le.etmin).and.(edep4.le.etmin)
	else
	  dist35 = 0.
	  dist45 = 0.
	  edep3 = (1.-xx3)/xx3*pt3
	  edep4 = (1.-xx4)/xx4*pt4
** 	  isol = (edep3.le.etmin).and.(edep4.le.etmin)
	endif
	x4min = pt4/(pt4+etmin)
	x3min = pt3/(pt3+etmin)
	if (ident) then	
	  if (pt3.ge.pt4) then
	    pt1 = pt3
	    pt2 = pt4
	    y1 = y3
	    y2 = y4
	  else
	    pt1 = pt4
	    pt2 = pt3
	    y1 = y4
	    y2 = y3
	  endif
	else
	  pt1 = pt3
	  pt2 = pt4
	  y1 = y3
	  y2 = y4
	endif
	edep = max(edep3,edep4)
c
	masscut = (xmgg.le.xmasmax.and.xmgg.ge.xmasmin)
c	
	rapidite = (y1.le.y1max.and.y1.ge.y1min).and.
     #	(y2.le.y2max.and.y2.ge.y2min)
	pete = pt1.ge.pt1min.and.pt2.ge.pt2min
c on convertit les variables
	wxmgg = sngl(xmgg)
	wpt3 = sngl(pt3)
	wpt4 = sngl(pt4)
	wqt = sngl(qt)
	wfi34 = sngl(fi34)
	wy3 = sngl(y3)
	wy4 = sngl(y4)
	wyboost = sngl(yboost)
	wystar = sngl(ystar)
	wygg = sngl(ygg)
	wp_out3 = sngl(p_out3)
	wp_out4 = sngl(p_out4)
	wz3_trig = sngl(z3_trig)
	wz4_trig = sngl(z4_trig)
	wcos_thetas = sngl(cos_thetas)
	wxx3 = sngl(xx3)
	wxx4 = sngl(xx4)
c
	if (rapidite.and.pete.and.isol.and.masscut) then
	  call remplissage(iprov,wxmgg,wpt3,wpt4,wqt,wfi34,
     #	  wy3,wy4,wyboost,wystar,wygg,wp_out3,wp_out4,wz3_trig,wz4_trig,
     #	  wcos_thetas,wxx3,wxx4,weight)
	endif
	end
