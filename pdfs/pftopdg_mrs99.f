c subroutine pdfset pour mrs99
	subroutine pdfset(parm,value)
	implicit real*8 (a-h,l-z)
        character*20 parm(20)
	dimension value(20)
        common/w50513/xmin,xmax,q2min,q2max
	common/w50512/qcdl4,qcdl5
	common/mrss/imode
	logical lo1,lo2
	ival1 = int(value(1)+1.d-3)
	ival2 = int(value(2)+1.d-3)
	ival3 = int(value(3)+1.d-3)
	lo1 = ival1.eq.1.and.ival2.eq.3
	lo2 = ival3.ge.199.and.ival3.le.210
	if (lo1.and.lo2) then
C  1     COR01  central gluon, a_s    300      0.1175   0.00537  C
C  2     COR02  higher gluon          300      0.1175   0.00497  C
C  3     COR03  lower gluon           300      0.1175   0.00398  C
C  4     COR04  lower a_s             229      0.1125   0.00585  C
C  5     COR05  higher a_s            383      0.1225   0.00384  C
C  6     COR06  quarks up             303.3    0.1178   0.00497  C
C  7     COR07  quarks down           290.3    0.1171   0.00593  C
C  8     COR08  strange up            300      0.1175   0.00524  C
C  9     COR09  strange down          300      0.1175   0.00524  C
C  10    C0R10  charm up              300      0.1175   0.00525  C
C  11    COR11  charm down            300      0.1175   0.00524  C
C  12    COR12  larger d/u            300      0.1175   0.00515  C
	  xmin = 1.d-5
	  xmax = 1.d0
	  q2min = 1.25d0
	  q2max = 1.d7
          imode = ival3-198
	  if (imode.eq.1.or.imode.eq.2.or.imode.eq.3) then
            qcdl4 = (.300)
	  else if (imode.eq.4) then
            qcdl4 = (.229)
	  else if (imode.eq.5) then
            qcdl4 = (.383)
	  else if (imode.eq.6) then
            qcdl4 = (.3033)
	  else if (imode.eq.7) then
            qcdl4 = (.2903)
	  else if (imode.eq.8.or.imode.eq.9.or.imode.eq.10.
     #	           or.imode.eq.11.or.imode.eq.12) then
            qcdl4 = (.300)
	  endif
c attention c'est completement faux
	  qcdl5 = qcdl4
c
	else
	  write (8,*) 'attention mauvais choix pour value',value
	endif
	return
	end
c subroutine pftopdg pour mrs99
	subroutine pftopdg(x,m,dxpdf)
	implicit real*8 (a-h,l-z)
	dimension dxpdf(-6:6)
	common/mrss/imode
	if (m*m.le.1.d+7) then
	  q = m
	else
	  q = dsqrt(1.d+7)
	endif
c
	call mrs99(x,q,imode,upv,dnv,usea,dsea,str,chm,bot,gl)
	dxpdf(-6) = 0.d0
	dxpdf(-5) = bot
	dxpdf(-4) = chm
	dxpdf(-3) = str
	dxpdf(-2) = usea
	dxpdf(-1) = dsea
	dxpdf(0) = gl
	dxpdf(1) = dnv+dsea
	dxpdf(2) = upv+usea
	dxpdf(3) = str
	dxpdf(4) = chm
	dxpdf(5) = bot
	dxpdf(6) = 0.d0	
	return
	end
