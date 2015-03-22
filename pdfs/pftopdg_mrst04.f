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
	lo2 = ival3.gt.0.and.ival3.lt.3
cccccccccccccccccccccccccccccccccccccccccccccccc	
C***************************************************************C
C                                                               C
C  This is a package for the new MRST 2004 'physical gluon' NLO C
C  and NNLO parton distributions.                               C
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0410230                                  C
C                                                               C
C  There are 2 pdf sets corresponding to mode = 1, 2            C
C                                                               C
C  Mode=1 gives the NLO set with Lambda(4) = 347 MeV            C
C  This set reads a grid called mrst2004nlo.dat                 C
C  whose first number is 0.00910                                C
C                                                               C
C  Mode=2 gives the NNLO set with Lambda(4) = 251 MeV           C
C  This set reads a grid called mrst2004nnlo.dat                C
C  whose first number is 0.00673                                C
C                                                               C
C  These fits use a new, physically motivated parametrisation   C
C  for the gluon at the starting scale, Q_0^2 = 1 GeV^2         C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
	if (lo1.and.lo2) then	
	  xmin = 1.d-5
	  xmax = 1.d0
	  q2min = 1.25d0
	  q2max = 1.d7
**           imode = ival3-199
          imode = ival3
	  if (imode.eq.1) then
            qcdl4 = (.347)
	  else 
	    write (8,*) 'DIPHOX is a NLO program' 
	    write (8,*) 'you cannot choose a NNLO pdf' 
	    stop
	  endif
c
	else
	  write (8,*) 'the set you chose for mrst04 does not exist'
	  write (8,*) 'ival1 = ', ival1
	  write (8,*) 'ival2 = ', ival2
	  write (8,*) 'ival3 = ', ival3
	endif
	return
	end
c subroutine pftopdg pour mrst04
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
	call mrst2004(x,q,imode,upv,dnv,usea,dsea,str,chm,bot,gl)
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
