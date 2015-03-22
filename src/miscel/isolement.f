	subroutine isolement(y4,y5,fi45,x4,pt5,gpt4,r_isol,etmax,
     #	isol)
	implicit real*8 (a-h,l-z)
	common/flag_iso/i_flag_iso
	logical isol
	r45s = (y4-y5)**2 + fi45**2
	r45 = dsqrt(r45s)
	if (r45.le.r_isol) then
	  et_dep = pt5 + (1.d0-x4)/x4*gpt4
	else
	  et_dep = (1.d0-x4)/x4*gpt4
	endif
c critere standard
	if (i_flag_iso.eq.1) then
	  etcut = etmax
c critere fraction du pt du photon
	else if (i_flag_iso.eq.2) then
	  etcut = etmax*gpt4
	endif
	if (et_dep.le.etcut) then
	  isol = .true.
	else
	  isol = .false.
	endif
	end
