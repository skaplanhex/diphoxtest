c*****************************************************
c fonction p3di+p4di
c*****************************************************	
c*****************************************************
c partie finie pour pt5 > ptm
c*****************************************************
c*****************************************************
	subroutine ssy3ro(s,gpt4,y3,y4,x10,x20,x3,z3,ih1,ih2,ih3,sy3r)
	implicit real*8 (a-h,l-z)
	common/coup/ptm,r
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part35r(k0max),sf(k0max)
	dimension sy3r(k0max),temp(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
c avec ce choix de variable x1p et x2p ne depende plus de z3
	x1p = x10
	x2p = x20
c juste pour tester (en principe inutile)
	if (x1p.gt.un.or.x2p.gt.un) then
	  call vec_dinit(sy3r,k0max,zero)
	  write (18,*) 'attention (r) x1p > 1 ou x2p > 1',
     #	              x1p,x2p
	  return
	endif
c attention les elements de matrice 2 -> 2 sont independant de 
c l'echelle d'energie
	sc = x10*x20*s
	tc = -x20*dsqrt(s)*gpt4*dexp(y4)
	uc = -x10*dsqrt(s)*gpt4*dexp(-y4)
	call strfrao(x1p,ih1,x2p,ih2,x3,ih3,sf)
	call spart35ro(sc,tc,uc,z3,r,part35r)
	call vec_dmult_vector(part35r,sf,k0max,temp)
	tyu = 1.d0/z3
	call vec_dmult_constant(temp,k0max,tyu,sy3r)
	return
	end
c**********************************************************************
	subroutine ssy4ro(s,gpt3,y3,y4,gpt4,x10,x20,x3,z4,ih1,ih2,
     #	ih3,sy4r)
	implicit real*8 (a-h,l-z)
	common/coup/ptm,r
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part45r(k0max),sf(k0max)
	dimension sy4r(k0max),temp(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt4 = gpt4
	x1s=x10/z4
	x2s=x20/z4
c juste pour tester (en principe inutile)
	if (x1s.gt.un.or.x2s.gt.un) then
	  call vec_dinit(sy4r,k0max,zero)
	  write (18,*) 'attention (r) x1s > 1 ou x2s > 1',
     #	              x1s,x2s
	  return
	endif
c attention les elements de matrice 2 -> 2 sont independant de 
c l'echelle d'energie
	sc = x10*x20*s
	tc = -x20*dsqrt(s)*pt4*dexp(y4)
	uc = -x10*dsqrt(s)*pt4*dexp(-y4)
	call strfrao(x1s,ih1,x2s,ih2,x3,ih3,sf)
	call spart45ro(sc,tc,uc,z4,r,part45r)
	tyu = 1.d0/z4
	call vec_dmult_vector(sf,part45r,k0max,temp)
	call vec_dmult_constant(temp,k0max,tyu,sy4r)
	return
	end
c*****************************************************
c fonctions pour p3di+p4di	
c*****************************************************	
	subroutine ssv12o(s,gpt3,y3,y4,m,mu,x10,x20,x3,ih1,ih2,
     #	ih3,sv12)
c integrale a 1 dimensions
	implicit real*8 (a-h,l-z)
	common/coup/ptmm,r
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part15z(k0max),part15l(k0max),part15d(k0max)
	dimension part25z(k0max),part25l(k0max),part25d(k0max)
	dimension partvi(k0max),sf(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	dimension au1(k0max),au2(k0max),au3(k0max),sv12(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt3 = gpt3/x3
	pt4 = pt3
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
	sc = x10*x20*s
	tc = -x10*dsqrt(s)*pt3*dexp(-y3)
	uc = -x20*dsqrt(s)*pt3*dexp(y3)
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf)
	call spart15zo(sc,tc,uc,un,m,ptm,part15z)
	call spart15lo(sc,tc,uc,un,part15l)
	call spart15do(sc,tc,uc,m,ptm,part15d)
	call spart25zo(sc,tc,uc,un,m,ptm,part25z)
	call spart25lo(sc,tc,uc,un,part25l)
	call spart25do(sc,tc,uc,m,ptm,part25d)
	call spartvio(s,sc,tc,uc,mu,pt3,ptm,partvi)
c
	call vec_dmult_vector(part15z,sf,k0max,temp0)
	call vec_dmult_vector(part15l,sf,k0max,temp1)
	call vec_dmult_vector(part15d,sf,k0max,temp2)
	ty1 = dlog(1-x10)
	ty2 = dlog(1-x10)**2/2.d0
	call vec_dmult_add(temp2,temp1,k0max,ty2,temp3)
	call vec_dmult_add(temp3,temp0,k0max,ty1,au1)
c
	call vec_dmult_vector(part25z,sf,k0max,temp0)
	call vec_dmult_vector(part25l,sf,k0max,temp1)
	call vec_dmult_vector(part25d,sf,k0max,temp2)
	ty1 = dlog(1-x20)
	ty2 = dlog(1-x20)**2/2.d0
	call vec_dmult_add(temp2,temp1,k0max,ty2,temp3)
	call vec_dmult_add(temp3,temp0,k0max,ty1,au2)
c
	call vec_dmult_vector(partvi,sf,k0max,au3)
c
	call vec_dadd_vector(au1,au2,k0max,temp0)
	ty3 = 1.d0
	call vec_dadd_mult_constant(temp0,au3,k0max,ty3,sv12)
	return
	end
c*****************************************************	
	subroutine ssv2o(s,gpt3,y3,y4,x10,x20,x3,u2,ih1,ih2,ih3,sv2)
c integrale a 2 dimensions
	implicit real*8 (a-h,l-z)
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension partvi2(k0max),sf(k0max)
	dimension sv2(k0max),temp(k0max)
	pi=4.d0*datan(1.d0)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	fi = pi*u2
	pt3 = gpt3/x3
	pt4 = pt3
c
	sc = x10*x20*s
	tc = -x10*dsqrt(s)*pt3*dexp(-y3)
	uc = -x20*dsqrt(s)*pt3*dexp(y3)
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf)
	call spartvi2o(s,sc,tc,uc,fi,partvi2)
c	call spartvi2o(s,sc,tc,uc,ys,fi,partvi2)
	call vec_dmult_vector(partvi2,sf,k0max,temp)
	ty3 = pi
	call vec_dmult_constant(temp,k0max,ty3,sv2)
	return
	end
c*****************************************************	
	subroutine ssv3o(s,gpt3,gpt4,y3,y4,mf,x10,x20,x3,ih1,ih2,
     #	ih3,sv3)
c integrale a 1 dimensions
	implicit real*8 (a-h,l-z)
	common/coup/ptmm,r
	common/flagexp/icoup_exp
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part35z(k0max),part35l(k0max),part35d(k0max),sf(k0max)
	dimension sv3(k0max),au3(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt3 = gpt3/x3
	pt4 = gpt4
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
	z3min = gpt3/gpt4
c
	sc = x10*x20*s
	tc = -x20*dsqrt(s)*gpt4*dexp(y4)
	uc = -x10*dsqrt(s)*gpt4*dexp(-y4)
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf)
	call spart35zo(sc,tc,uc,un,mf,pt3,part35z)
	call spart35lo(sc,tc,uc,un,part35l)
	call spart35do(sc,tc,uc,mf,pt3,part35d)
	call vec_dmult_vector(part35z,sf,k0max,temp0)
	call vec_dmult_vector(part35l,sf,k0max,temp1)
	call vec_dmult_vector(part35d,sf,k0max,temp2)
	ty1 = dlog(1-z3min)
	ty2 = dlog(1-z3min)**2/2.d0
	call vec_dmult_add(temp2,temp1,k0max,ty2,temp3)
	call vec_dmult_add(temp3,temp0,k0max,ty1,au3)
	ty3 = 1.d0
	call vec_dmult_constant(au3,k0max,ty3,sv3)
	return
	end
c*****************************************************	
	subroutine ssv4o(s,gpt3,gpt4,etmax,y3,y4,mf,x10,x20,x3,ih1,ih2,
     #	ih3,sv4)
c integrale a 1 dimensions
	implicit real*8 (a-h,l-z)
	common/coup/ptmm,r
	common/flagexp/icoup_exp
	common/flag_iso/i_flag_iso
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part45z(k0max),part45l(k0max),part45d(k0max),sf(k0max)
	dimension sv4(k0max),au4(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
** 	gpt4 = gpt3/x3
	pt4 = gpt4
	pt3 = gpt3
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
c z4 doit etre < a 2*gpt4/dsqrts)*dcosh(y3) et 2*gpt4/dsqrts)*dcosh(y4)
	x4min = 2.d0*gpt4/dsqrt(s)*dcosh(y4)
	x3min = 2.d0*gpt4/dsqrt(s)*dcosh(y3)
	if (i_flag_iso.eq.1) then
	  cut_isol1 = gpt4/(gpt3+etmax)
	  cut_isol2 = gpt4/(gpt4+etmax)
	else if (i_flag_iso.eq.2) then
	  cut_isol1 = gpt4/(gpt3*(1.d0+etmax))
	  cut_isol2 = 1.d0/(1.d0+etmax)
	endif
	z4min = dmax1(x4min,x3min,x10,x20,cut_isol1,cut_isol2)
c
	sc = x10*x20*s
	tc = -x20*dsqrt(s)*pt4*dexp(y4)
	uc = -x10*dsqrt(s)*pt4*dexp(-y4)
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf)
	call spart45zo(sc,tc,uc,un,mf,pt4,part45z)
	call spart45lo(sc,tc,uc,un,part45l)
	call spart45do(sc,tc,uc,mf,pt4,part45d)
	call vec_dmult_vector(part45z,sf,k0max,temp0)
	call vec_dmult_vector(part45l,sf,k0max,temp1)
	call vec_dmult_vector(part45d,sf,k0max,temp2)
	ty1 = dlog(1-z4min)
	ty2 = dlog(1-z4min)**2/2.d0
	call vec_dmult_add(temp2,temp1,k0max,ty2,temp3)
	call vec_dmult_add(temp3,temp0,k0max,ty1,au4)
	ty4 = 1.d0
	call vec_dmult_constant(au4,k0max,ty4,sv4)
	return
	end
c*****************************************************	
	subroutine ssy3o(s,gpt3,y3,y4,m,x10,x20,x3,u2,ih1,ih2,ih3,sy3)
	implicit real*8 (a-h,l-z)
	common/coup/ptmm,r
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part15z(k0max),part15l(k0max)
	dimension part15z0(k0max),part15l0(k0max)
	dimension part25z(k0max),part25l(k0max)
	dimension part25z0(k0max),part25l0(k0max)
	dimension sf(k0max),sf0(k0max)
	dimension sy3(k0max),bu1(k0max),bu2(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	dimension temp4(k0max),temp5(k0max),temp6(k0max),temp7(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt3 = gpt3/x3
	pt5x = dsqrt(s)/(2.d0*dcosh(ys))*(1.d0-x10)*(1.d0-x20)/
     #  ((1.d0-x10)*dexp(-y)+(1.d0-x20)*dexp(y))
	ptm = dmin1(ptmm,pt5x)
cc	ptm = ptmm
	sc = x10*x20*s
	tc = -x10*dsqrt(s)*pt3*dexp(-y3)
	uc = -x20*dsqrt(s)*pt3*dexp(y3)
	z1max = 1.d0
	z1min = x10
	z1 = z1min + (z1max-z1min)*u2
	call strfrao(x10/z1,ih1,x20,ih2,x3,ih3,sf)
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf0)
	call spart15zo(sc,tc,uc,z1,m,ptm,part15z)
	call spart15zo(sc,tc,uc,un,m,ptm,part15z0)
	call spart15lo(sc,tc,uc,z1,part15l)
	call spart15lo(sc,tc,uc,un,part15l0)
	call vec_dmult_vector(part15z,sf,k0max,temp0)
	call vec_dmult_vector(part15z0,sf0,k0max,temp1)
	ty1 = 1.d0/z1**2
	call vec_dmult_sub(temp1,temp0,k0max,ty1,temp2)
	ty2 = 1.d0/(1.d0-z1)
	call vec_dmult_constant(temp2,k0max,ty2,temp3)
c
	call vec_dmult_vector(part15l,sf,k0max,temp4)
	call vec_dmult_vector(part15l0,sf0,k0max,temp5)
	call vec_dmult_sub(temp5,temp4,k0max,ty1,temp6)
	ty3 = dlog(1.d0-z1)/(1.d0-z1)
	call vec_dmult_constant(temp6,k0max,ty3,temp7)
c
	ty4 = (z1max-z1min)
	call vec_dadd_mult_constant(temp3,temp7,k0max,ty4,bu1)
	z2max = 1.d0
	z2min = x20
	z2 = z2min + (z2max-z2min)*u2
	call strfrao(x10,ih1,x20/z2,ih2,x3,ih3,sf)
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf0)
	call spart25zo(sc,tc,uc,z2,m,ptm,part25z)
	call spart25zo(sc,tc,uc,un,m,ptm,part25z0)
	call spart25lo(sc,tc,uc,z2,part25l)
	call spart25lo(sc,tc,uc,un,part25l0)
	call vec_dmult_vector(part25z,sf,k0max,temp0)
	call vec_dmult_vector(part25z0,sf0,k0max,temp1)
	ty1 = 1.d0/z2**2
	call vec_dmult_sub(temp1,temp0,k0max,ty1,temp2)
	ty2 = 1.d0/(1.d0-z2)
	call vec_dmult_constant(temp2,k0max,ty2,temp3)
c
	call vec_dmult_vector(part25l,sf,k0max,temp4)
	call vec_dmult_vector(part25l0,sf0,k0max,temp5)
	call vec_dmult_sub(temp5,temp4,k0max,ty1,temp6)
	ty3 = dlog(1.d0-z2)/(1.d0-z2)
	call vec_dmult_constant(temp6,k0max,ty3,temp7)
c
	ty4 = (z2max-z2min)
	call vec_dadd_mult_constant(temp3,temp7,k0max,ty4,bu2)
c
	call vec_dadd_vector(bu1,bu2,k0max,sy3)
	return
	end
c*****************************************************	
	subroutine ssy3po(s,gpt3,gpt4,y3,y4,mf,x10,x20,x3,z3,ih1,ih2,
     #	ih3,sy3p)
	implicit real*8 (a-h,l-z)
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part35z(k0max),part35l(k0max),sf(k0max)
	dimension part35z0(k0max),part35l0(k0max),sf0(k0max)
	dimension sy3p(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	dimension temp4(k0max),temp5(k0max),temp6(k0max),temp7(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt3 = z3*gpt4
	x30 = gpt3/gpt4
c avec ce choix de variable x1p et x2p ne depende plus de z3
	x1p = x10
	x2p = x20
c juste pour tester (en principe inutile)
	if (x1p.gt.un.or.x2p.gt.un) then
	  call vec_dinit(sy3p,k0max,zero)
	  write (18,*) 'attention x1p > 1 ou x2p > 1',
     #	              x1p,x2p
	  return
	endif
c attention les elements de matrice 2 -> 2 sont independant de 
c l'echelle d'energie
	sc = x10*x20*s
	tc = -x20*dsqrt(s)*gpt4*dexp(y4)
	uc = -x10*dsqrt(s)*gpt4*dexp(-y4)
	call strfrao(x1p,ih1,x2p,ih2,x3,ih3,sf)
	call strfrao(x10,ih1,x20,ih2,x30,ih3,sf0)
	call spart35zo(sc,tc,uc,z3,mf,pt3,part35z)
	call spart35zo(sc,tc,uc,un,mf,gpt4,part35z0)
	call spart35lo(sc,tc,uc,z3,part35l)
	call spart35lo(sc,tc,uc,un,part35l0)
c
	call vec_dmult_vector(part35z,sf,k0max,temp0)
	call vec_dmult_vector(part35z0,sf0,k0max,temp1)
	ty1 = 1.d0/z3
	call vec_dmult_sub(temp1,temp0,k0max,ty1,temp2)
	ty2 = 1.d0/(1.d0-z3)
	call vec_dmult_constant(temp2,k0max,ty2,temp3)
c
	call vec_dmult_vector(part35l,sf,k0max,temp4)
	call vec_dmult_vector(part35l0,sf0,k0max,temp5)
	call vec_dmult_sub(temp5,temp4,k0max,ty1,temp6)
	ty3 = dlog(1.d0-z3)/(1.d0-z3)
	call vec_dmult_constant(temp6,k0max,ty3,temp7)
c
	ty4 = un
	call vec_dadd_mult_constant(temp3,temp7,k0max,ty4,sy3p)
	return
	end
c*****************************************************	
	subroutine ssy4po(s,gpt3,y3,y4,mf,gpt4,x10,x20,x3,z4,ih1,ih2,
     #	ih3,sy4p)
	implicit real*8 (a-h,l-z)
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part45z(k0max),part45l(k0max),sf(k0max)
	dimension part45z0(k0max),part45l0(k0max),sf0(k0max)
	dimension sy4p(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	dimension temp4(k0max),temp5(k0max),temp6(k0max),temp7(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt4 = gpt4
	x1s = x10/z4
	x2s = x20/z4
c juste pour tester (en principe inutile)
	if (x1s.gt.un.or.x2s.gt.un) then
	  call vec_dinit(sy4p,k0max,zero)
	  write (18,*) 'attention x1s > 1 ou x2s > 1',
     #	              x1s,x2s
	  return
	endif
c attention les elements de matrice 2 -> 2 sont independant de 
c l'echelle d'energie
	sc = x10*x20*s
	tc = -x20*dsqrt(s)*pt4*dexp(y4)
	uc = -x10*dsqrt(s)*pt4*dexp(-y4)
	call strfrao(x1s,ih1,x2s,ih2,x3,ih3,sf)
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf0)
	call spart45zo(sc,tc,uc,z4,mf,pt4,part45z)
	call spart45zo(sc,tc,uc,un,mf,pt4,part45z0)
	call spart45lo(sc,tc,uc,z4,part45l)
	call spart45lo(sc,tc,uc,un,part45l0)
c
	call vec_dmult_vector(part45z,sf,k0max,temp0)
	call vec_dmult_vector(part45z0,sf0,k0max,temp1)
	ty1 = 1.d0/z4
	call vec_dmult_sub(temp1,temp0,k0max,ty1,temp2)
	ty2 = 1.d0/(1.d0-z4)
	call vec_dmult_constant(temp2,k0max,ty2,temp3)
c
	call vec_dmult_vector(part45l,sf,k0max,temp4)
	call vec_dmult_vector(part45l0,sf0,k0max,temp5)
	call vec_dmult_sub(temp5,temp4,k0max,ty1,temp6)
	ty3 = dlog(1.d0-z4)/(1.d0-z4)
	call vec_dmult_constant(temp6,k0max,ty3,temp7)
c
	ty4 = 1.d0
	call vec_dadd_mult_constant(temp3,temp7,k0max,ty4,sy4p)
	return
	end
c*****************************************************	
	subroutine ssy4pop(s,gpt3,y3,y4,mf,gpt4,x10,x20,x3,z4,ih1,ih2,
     #	ih3,sy4p)
	implicit real*8 (a-h,l-z)
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension part45z(k0max),part45l(k0max),sf(k0max)
	dimension sy4p(k0max)
	dimension temp0(k0max),temp1(k0max),temp2(k0max),temp3(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt4 = gpt4
	x1s = x10/z4
	x2s = x20/z4
c juste pour tester (en principe inutile)
	if (x1s.gt.un.or.x2s.gt.un) then
	  call vec_dinit(sy4p,k0max,zero)
	  write (18,*) 'attention (p) x1s > 1 ou x2s > 1',
     #	              x1s,x2s
	  return
	endif
c attention les elements de matrice 2 -> 2 sont independant de 
c l'echelle d'energie
	sc = x10*x20*s
	tc = -x20*dsqrt(s)*pt4*dexp(y4)
	uc = -x10*dsqrt(s)*pt4*dexp(-y4)
	call strfrao(x1s,ih1,x2s,ih2,x3,ih3,sf)
	call spart45zo(sc,tc,uc,z4,mf,pt4,part45z)
	call spart45lo(sc,tc,uc,z4,part45l)
	call vec_dmult_vector(part45z,sf,k0max,temp0)
	ty2 = 1.d0/(1.d0-z4)
	call vec_dmult_constant(temp0,k0max,ty2,temp1)
c
	call vec_dmult_vector(part45l,sf,k0max,temp2)
	ty3 = dlog(1.d0-z4)/(1.d0-z4)
	call vec_dmult_constant(temp2,k0max,ty3,temp3)
c
	ty4 = 1.d0/z4
	call vec_dadd_mult_constant(temp1,temp3,k0max,ty4,sy4p)
	return
	end
c*****************************************************
c terme a l'ordre de born	
c*****************************************************	
	subroutine sfborno(s,n,gpt3,y3,y4,x3,ih1,ih2,ih3,fborn)
	implicit real*8 (a-h,l-z)
	common/scale/m,mf,mu
	common/aurenche/iauren
cc	parameter (k0max=72)
	parameter (k0max=36)
	dimension sf1(k0max),rborn1(k0max)
	dimension fborn(k0max),temp0(k0max)
	zero = 0.d0
	un = 1.d0
	ys = (y3-y4)/2.d0
	y = (y3+y4)/2.d0
	pt = gpt3/x3
	x10 = 2.d0*pt/dsqrt(s)*dexp(y)*dcosh(ys)	
	x20 = 2.d0*pt/dsqrt(s)*dexp(-y)*dcosh(ys)
	if (x10.gt.un.or.x20.gt.un) then
	  call vec_dinit(fborn,k0max,zero)
	  return
	endif
	sc = 4.d0*pt**2*dcosh(ys)**2
	tc = -2.d0*pt**2*dcosh(ys)*dexp(-ys)
	uc = -2.d0*pt**2*dcosh(ys)*dexp(ys)
	if (iauren.eq.1) then
ccc
	  mf = dsqrt(sc)
ccc
	endif
	call strfrao(x10,ih1,x20,ih2,x3,ih3,sf1)
	call borno(sc,tc,uc,n,rborn1)
	call vec_dmult_vector(rborn1,sf1,k0max,temp0)
	ty3 = 1.d0
	call vec_dmult_constant(temp0,k0max,ty3,fborn)
	return
	end 
c	
