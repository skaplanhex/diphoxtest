c analog to pftopdg_mrs99 for mrst01, 31.5.02
c fixes qcdl4 
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
	lo2 = ival3.gt.199.and.ival3.lt.204
cccccccccccccccccccccccccccccccccccccccccccccccc
c  !! ival3=199+Mode !!
cccccccccccccccccccccccccccccccccccccccccccccccc	
C  Mode=1 gives the default set with Lambda(MSbar,nf=4) = 0.323 C
C  corresponding to alpha_s(M_Z) of 0.119                       C
C  This set reads a grid whose first number is 0.00927          C
C                                                               C
C  Mode=2 gives the set with Lambda(MSbar,nf=4) = 0.290         C
C  corresponding to alpha_s(M_Z) of 0.117                       C
C  This set reads a grid whose first number is 0.00953          C
C                                                               C
C  Mode=3 gives the set with Lambda(MSbar,nf=4) = 0.362         C
C  corresponding to alpha_s(M_Z) of 0.121                       C
C  This set reads a grid whose first number is 0.00889          C
C                                                               C
C  Mode=4 gives the set MRST2001J which gives better agreement  C
C  with the Tevatron inclusive jet data but has unattractive    C
C  gluon behaviour at large x (see discussion in paper)         C
C  This set has Lambda(MSbar,nf=4) = 0.353(alpha_s(M_Z) = 0.121 C 
C  This set reads a grid whose first number is 0.00826          C
	if (lo1.and.lo2) then	
	  xmin = 1.d-5
	  xmax = 1.d0
	  q2min = 1.25d0
	  q2max = 1.d7
          imode = ival3-199
	  if (imode.eq.1) then
            qcdl4 = (.323)
	  else if (imode.eq.2) then
            qcdl4 = (.290)
	  else if (imode.eq.3) then
            qcdl4 = (.362)
	  else if (imode.eq.4) then
            qcdl4 = (.353)
	  endif
c
	else
	  write (8,*) 'the set you chose for mrst01 does not exist'
	  write (8,*) 'ival1 = ', ival1
	  write (8,*) 'ival2 = ', ival2
	  write (8,*) 'ival3 = ', ival3
	endif
	return
	end
c subroutine pftopdg pour mrst01
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
	call mrst2001(x,q,imode,upv,dnv,usea,dsea,str,chm,bot,gl)
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
