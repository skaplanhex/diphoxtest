	subroutine choiscale(p3,p4,cm,cmu,cmf,ichoi_scale,m,mu,mf)
	implicit real*8 (a-h,l-z)
	dimension p3(4),p4(4)
	pt3 = dsqrt(p3(2)**2+p3(3)**2)
	pt4 = dsqrt(p4(2)**2+p4(3)**2)
	s1 = (pt3+pt4)
	s2 = dsqrt(pt3**2+pt4**2)
	s3 = dsqrt(2.d0*sca(p3,p4))
	s4 = dmax1(pt3,pt4)
	if (ichoi_scale.eq.1) then
	  m = s1*cm
	  mu = s1*cmu
	  mf = s1*cmf
	else if (ichoi_scale.eq.2) then
	  m = s2*cm
	  mu = s2*cmu
	  mf = s2*cmf
	else if (ichoi_scale.eq.3) then
	  m = s3*cm
	  mu = s3*cmu
	  mf = s3*cmf
	else if (ichoi_scale.eq.4) then
	  m = s4*cm
	  mu = s4*cmu
	  mf = s4*cmf
	endif
	return
	end
	
