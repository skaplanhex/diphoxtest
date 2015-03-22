***********************************************************************
*                                                                     *
*    Program for generating charged Current Structure Functions using *
*    consistent treatment of charm and bottom structure functions     *
*    Charm mass = 1.35 GeV   Bottom mass = 4.3 GeV                    *
*                                                                     *
*    The program should be run only with iord set to 1                *
*    The program is self contained, only requiring the subroutine     *
*    mrst2002.f and the grid file mrst2002nlo.dat to be accessible    *
*                                                                     *
***********************************************************************

      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      COMMON/INPUT/alambda,flavor,qsct,qsdt,iord
      dimension xccfr2(11),add(11)
      data xccfr2/0.015,0.045,0.08,0.125,0.175,0.225,0.275,
     .0.35,0.45,0.55,0.65/
      data add/0.4,0.3,0.2,8*0.1/

      call wate96
      iord=1
      flavor=3
      qsdt=8.18
      qsct=74.
      emw=80.41
      mode=1
      alambda=0.3342
   99 read(5,*) x,q2

      call sfun(x,q2,mode,f2wpp,f3wpp,f2wmp,f3wmp,f2wpn,f3wpn,
     .f2wmn,f3wmn,flwpp,flwmp,flwpn,flwmn)      

      print 100,x,q2,f2wpp,f3wpp,f2wmp,f3wmp,flwpp,flwmp,
     .f2wpn,f3wpn,f2wmn,f3wmn,flwpn,flwmn

      f2ccfr=0.25*(f2wpp+f2wpn+f2wmp+f2wmn)
      f3ccfr=0.25*(f3wpp+f3wpn+f3wmp+f3wmn)
      delxf3=0.25*(f3wpp+f3wpn-f3wmp-f3wmn)
c      print 98, x,q2,f2ccfr,f3ccfr

   98 format(8(1x,e12.6))
      go to 99
  100 format(1x,'x= ',f10.6,' Q2= ',f14.3//
     11x,' Proton Structure Functions '/
     25x,' F2(W+p)= ',f9.5,' xF3(W+p)= ',f9.5,
     3' F2(W-p)= ',f9.5,' xF3(W-p)= ',f9.5/
     25x,' FL(W+p)= ',f8.4,' FL(W-p)= ',f8.4//
     41x,' Neutron Structure Functions '/
     55x,' F2(W+n)= ',f8.4,' xF3(W+n)= ',f8.4,
     6' F2(W-n)= ',f8.4,' xF3(W-n)= ',f8.4/
     25x,' FL(W+n)= ',f8.4,' FL(W-n)= ',f8.4)
  234 format(10x,2(1x,e10.4))
      end


      subroutine sfun(x,q2,mode,f2wpp,f3wpp,f2wmp,f3wmp,
     .f2wpn,f3wpn,f2wmn,f3wmn,flwpp,flwmp,flwpn,flwmn)

CC  NOTATION

cc  f2wpp is the  F2 str. fn. for W+ proton scattering  
cc  f3wpp is the xF3 str. fn. for W+ proton scattering  
cc  f2wmp is the  F2 str. fn. for W- proton scattering  
cc  f3wpp is the xF3 str. fn. for W- proton scattering  
cc  f2wpn is the  F2 str. fn. for W+ neutron scattering  
cc  f3wpn is the xF3 str. fn. for W+ neutron scattering  
cc  f2wmn is the  F2 str. fn. for W- neutron scattering  
cc  f3wpn is the xF3 str. fn. for W- neutron scattering
   
cc  To compute the quantities for neutrino scattering off heavy 
cc  isoscalar target, for example in the CCFR experiment, take the
cc  average of all four F2 values to get F2 and likewise for xF3,
cc  then apply a heavy target correction (EMC effect). A possible
cc  correction factor is given in MRST.

      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      COMMON/INPUT/alambda,flavor,qsct,qsdt,iord
      data pi,pi2/3.14159,9.8696/
      XLAM=alambda
      xlam2=xlam*xlam
      t=dlog(q2/xlam2)
      al=alpha(t)/(4.*pi)
      scale=dsqrt(q2)

      sc2=0.22*0.22
      cc2=1.-sc2
      cf=4./3.
      ca=3.
      enf=flavor

      ww2=(1.-x)*q2/x
      epsc4=qsdt/q2
      epsc=epsc4/4.
      xcc=(1.+epsc)*x
      if(xcc.gt.0.95) xcc=0.95
      f2wpp=0.
      f3wpp=0.
      f2wmp=0.
      f3wmp=0.
      f2wpn=0.
      f3wpn=0.
      f2wmn=0.
      f3wmn=0.
      flwpp=0.
      flwmp=0.
      flwpn=0.
      flwmn=0.

      call mrst2002(x,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      if(epsc.le.1.) chm=(1.-epsc)*chm
      f2wppx=2.*((dnv+dsea)+str+usea)+2.*chm
      f3wppx=2.*((dnv+dsea)+str-usea)-2.*chm
      f2wmpx=2.*( dsea+str+(upv+usea))+2.*chm
      f3wmpx=2.*(-dsea-str+(upv+usea))+2.*chm
      f2wppx0=2.*(cc2*(dnv+dsea)+sc2*str+usea)
      f3wppx0=2.*(cc2*(dnv+dsea)+sc2*str-usea)
      f2wmpx0=2.*( cc2*dsea+sc2*str+(upv+usea))
      f3wmpx0=2.*(-cc2*dsea-sc2*str+(upv+usea))
      f2wpnx=2.*((upv+usea)+str+dsea)+2.*chm
      f3wpnx=2.*((upv+usea)+str-dsea)-2.*chm
      f2wmnx=2.*( usea+str+(dnv+dsea))+2.*chm
      f3wmnx=2.*(-usea-str+(dnv+dsea))+2.*chm
      f2wpnx0=2.*(cc2*(upv+usea)+sc2*str+dsea)
      f3wpnx0=2.*(cc2*(upv+usea)+sc2*str-dsea)
      f2wmnx0=2.*( cc2*usea+sc2*str+(dnv+dsea))
      f3wmnx0=2.*(-cc2*usea-sc2*str+(dnv+dsea))

      call mrst2002(xcc,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      f2wppx0=f2wppx0+2.*(sc2*(dnv+dsea)+cc2*str)
      f3wppx0=f3wppx0+2.*(sc2*(dnv+dsea)+cc2*str)*x/xcc
      f2wmpx0=f2wmpx0+2.*( sc2*dsea+cc2*str)
      f3wmpx0=f3wmpx0+2.*(-sc2*dsea-cc2*str)*x/xcc
      f2wpnx0=f2wpnx0+2.*(sc2*(upv+usea)+cc2*str)
      f3wpnx0=f3wpnx0+2.*(sc2*(upv+usea)+cc2*str)*x/xcc
      f2wmnx0=f2wmnx0+2.*( sc2*usea+cc2*str)
      f3wmnx0=f3wmnx0+2.*(-sc2*usea-cc2*str)*x/xcc

      AL1=dLOG(1.-X)

      f2wpp=f2wppx0+f2wppx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      f3wpp=f3wppx0+f3wppx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      f2wmp=f2wmpx0+f2wmpx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      f3wmp=f3wmpx0+f3wmpx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      f2wpn=f2wpnx0+f2wpnx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      f3wpn=f3wpnx0+f3wpnx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      f2wmn=f2wmnx0+f2wmnx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      f3wmn=f3wmnx0+f3wmnx*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))


      DO 23 I=1,NTERMS
      Y=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      XY=X/Y
      AL1=dLOG(1.-Y)
      call mrst2002(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      if(epsc.le.1.) chm=(1.-epsc)*chm
      f2wppxy=2.*((dnv+dsea)+str+usea)+2.*chm
      f3wppxy=2.*((dnv+dsea)+str-usea)-2.*chm
      f2wmpxy=2.*( dsea+str+(upv+usea))+2.*chm
      f3wmpxy=2.*(-dsea-str+(upv+usea))+2.*chm
      f2wpnxy=2.*((upv+usea)+str+dsea)+2.*chm
      f3wpnxy=2.*((upv+usea)+str-dsea)-2.*chm
      f2wmnxy=2.*( usea+str+(dnv+dsea))+2.*chm
      f3wmnxy=2.*(-usea-str+(dnv+dsea))+2.*chm
      gluxy=glu
      C222=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*dLOG(Y)-2.*(1.+Y)*AL1)
      C223=C222-2.*(1.+Y)
      C23=CF*(-3.+4.*AL1)/(1.-Y)
      CG2=4.*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*dLOG(1./Y-1.))

      f1lq=4.*cf*y
      f1lg=8.*enf*y*(1.-y)

      f2wpp=f2wpp
     ++0.5*(1.-X)*WI(I)*AL*(C222*f2wppxy+C23*(f2wppxy-f2wppx))
      f3wpp=f3wpp
     ++0.5*(1.-X)*WI(I)*AL*(C223*f3wppxy+C23*(f3wppxy-f3wppx))
      f2wmp=f2wmp
     ++0.5*(1.-X)*WI(I)*AL*(C222*f2wmpxy+C23*(f2wmpxy-f2wmpx))
      f3wmp=f3wmp
     ++0.5*(1.-X)*WI(I)*AL*(C223*f3wmpxy+C23*(f3wmpxy-f3wmpx))
      f2wpn=f2wpn
     ++0.5*(1.-X)*WI(I)*AL*(C222*f2wpnxy+C23*(f2wpnxy-f2wpnx))
      f3wpn=f3wpn
     ++0.5*(1.-X)*WI(I)*AL*(C223*f3wpnxy+C23*(f3wpnxy-f3wpnx))
      f2wmn=f2wmn
     ++0.5*(1.-X)*WI(I)*AL*(C222*f2wmnxy+C23*(f2wmnxy-f2wmnx))
      f3wmn=f3wmn
     ++0.5*(1.-X)*WI(I)*AL*(C223*f3wmnxy+C23*(f3wmnxy-f3wmnx))

      flwpp=flwpp+0.5*(1.-x)*wi(i)*al*f1lq*f2wpp
      flwmp=flwmp+0.5*(1.-x)*wi(i)*al*f1lq*f2wmp
      flwpn=flwpn+0.5*(1.-x)*wi(i)*al*f1lq*f2wpn
      flwmn=flwmn+0.5*(1.-x)*wi(i)*al*f1lq*f2wmn

   24 CONTINUE

      f2wpp=f2wpp+0.5*(1.-X)*WI(I)*AL*CG2*gluxy
      f2wmp=f2wmp+0.5*(1.-X)*WI(I)*AL*CG2*gluxy
      f2wpn=f2wpn+0.5*(1.-X)*WI(I)*AL*CG2*gluxy
      f2wmn=f2wmn+0.5*(1.-X)*WI(I)*AL*CG2*gluxy

      flwpp=flwpp+0.5*(1.-x)*wi(i)*al*f1lg*gluxy
      flwmp=flwmp+0.5*(1.-x)*wi(i)*al*f1lg*gluxy
      flwpn=flwpn+0.5*(1.-x)*wi(i)*al*f1lg*gluxy
      flwmn=flwmn+0.5*(1.-x)*wi(i)*al*f1lg*gluxy

   23 CONTINUE
   21 CONTINUE

      r7=dsqrt(7.d0)
      xcmax=1./(1.+epsc)
      if(xcmax.le.x) go to 321
      call mrst2002(xcc,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1) chm=0.
      gluxcc=glu
      if(epsc.gt.1) gluxcc=0.
      f2wppxcc= 2.*chm
      f3wppxcc=-2.*chm
      f2wmpxcc= 2.*chm
      f3wmpxcc= 2.*chm
      del=0.01
      xccu=xcc*(1.+0.5*del)
      if(xccu.gt.1.d0) then
      df2wppx=0.
      df3wppx=0.
      df2wmpx=0.
      df3wmpx=0.
      go to 1325
      endif
      call mrst2002(xccu,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      f2wppxu= 2.*chm
      f3wppxu=-2.*chm
      f2wmpxu= 2.*chm
      f3wmpxu= 2.*chm
      xccl=xcc*(1.-0.5*del)
      call mrst2002(xccl,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      f2wppxl= 2.*chm
      f3wppxl=-2.*chm
      f2wmpxl= 2.*chm
      f3wmpxl= 2.*chm
      df2wppxcc=(f2wppxu-f2wppxl)/del
      df3wppxcc=(f3wppxu-f3wppxl)/del
      df2wmpxcc=(f2wmpxu-f2wmpxl)/del
      df3wmpxcc=(f3wmpxu-f3wmpxl)/del
      c0cfac=-epsc*xcmax*dlog(epsc*xcmax*xcmax)
 1325 c0c2l=(2.*epsc*xcmax*dlog(1.-xcc)+c0cfac)*(1.
     .-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc))
      c0c3l=c0c2l*xcmax
      f2wpp=f2wpp+c0c2l*(-df2wppxcc+3.*f2wppxcc)
      f2wpp=f2wpp-al*c0c2l*gluxcc*2.*dlog(1./epsc)     
      f2wmp=f2wmp+c0c2l*(-df2wmpxcc+3.*f2wmpxcc)
      f2wmp=f2wmp-al*c0c2l*gluxcc*2.*dlog(1./epsc)     
      f2wpn=f2wpn+c0c2l*(-df2wppxcc+3.*f2wppxcc)
      f2wpn=f2wpn-al*c0c2l*gluxcc*2.*dlog(1./epsc)     
      f2wmn=f2wmn+c0c2l*(-df2wmpxcc+3.*f2wmpxcc)
      f2wmn=f2wmn-al*c0c2l*gluxcc*2.*dlog(1./epsc)     
      f3wpp=f3wpp-c0c3l*(-df3wppxcc+3.*f3wppxcc)
      f3wpp=f3wpp-al*c0c3l*gluxcc*2.*dlog(1./epsc)     
      f3wmp=f3wmp-c0c3l*(-df3wmpxcc+3.*f3wmpxcc)
      f3wmp=f3wmp+al*c0c3l*gluxcc*2.*dlog(1./epsc)     
      f3wpn=f3wpn-c0c3l*(-df3wppxcc+3.*f3wppxcc)
      f3wpn=f3wpn-al*c0c3l*gluxcc*2.*dlog(1./epsc)     
      f3wmn=f3wmn-c0c3l*(-df3wmpxcc+3.*f3wmpxcc)
      f3wmn=f3wmn+al*c0c3l*gluxcc*2.*dlog(1./epsc)     

      DO 323 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      ycc=(1.+epsc)*y
      call mrst2002(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      gluxy=glu
      f2wppxy= 2.*chm
      f3wppxy=-2.*chm
      f2wmpxy= 2.*chm
      f3wmpxy= 2.*chm
      del=0.01
      xyu=xy*(1.+0.5*del)
      if(xyu.gt.1.d0) then
      df2wppxy=0.
      df3wppxy=0.
      df2wmpxy=0.
      df3wmpxy=0.
      go to 1324
      endif
      call mrst2002(xyu,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      f2wppxyu= 2.*chm
      f3wppxyu=-2.*chm
      f2wmpxyu= 2.*chm
      f3wmpxyu= 2.*chm
      xyl=xy*(1.-0.5*del)
      call mrst2002(xyl,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      f2wppxyl= 2.*chm
      f3wppxyl=-2.*chm
      f2wmpxyl= 2.*chm
      f3wmpxyl= 2.*chm
      df2wppxy=(f2wppxyu-f2wppxyl)/del
      df3wppxy=(f3wppxyu-f3wppxyl)/del
      df2wmpxy=(f2wmpxyu-f2wmpxyl)/del
      df3wmpxy=(f3wmpxyu-f3wmpxyl)/del
 1324 c0c2=cheavy(4,y,epsc)
      if(epsc.gt.1.d0) c0c2=0.d0
      c0c3=cheavy(6,y,epsc)
      if(epsc.gt.1.d0) c0c3=0.d0
      cg21c=2.*cheavy(3,y,epsc)
      cg22c=2.*c0c2*dlog(1./epsc) 
      cg31c=2.*cheavy(5,y,epsc)
      cg32c=2.*c0c3*dlog(1./epsc)
c      print 333,y,al,cg21c,cg22c,c0c2,c0c3
c  333 format(1x,f10.6,6(1x,e9.3))
      clg2c=2.*clheavy(3,y,epsc)
      f1lq=clheavy(4,y,epsc)

      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*1.*c0c2*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df2wppxy+3.*f2wppxy)
      f3wpp=f3wpp-0.5*(xcmax-x)*wi(i)*1.*c0c3*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df3wppxy+3.*f3wppxy)

      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*1.*c0c2*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df2wmpxy+3.*f2wmpxy)
      f3wmp=f3wmp-0.5*(xcmax-x)*wi(i)*1.*c0c3*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df3wmpxy+3.*f3wmpxy)

      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*1.*c0c2*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df2wppxy+3.*f2wppxy)
      f3wpn=f3wpn-0.5*(xcmax-x)*wi(i)*1.*c0c3*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df3wppxy+3.*f3wppxy)

      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*1.*c0c2*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df2wmpxy+3.*f2wmpxy)
      f3wmn=f3wmn-0.5*(xcmax-x)*wi(i)*1.*c0c3*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)
     .-1.6-2*xy))*(-df3wmpxy+3.*f3wmpxy)

      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      f3wpp=f3wpp+0.5*(xcmax-x)*wi(i)*al*(cg31c-cg32c)*gluxy
      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      f3wmp=f3wmp-0.5*(xcmax-x)*wi(i)*al*(cg31c-cg32c)*gluxy
      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      f3wpn=f3wpn+0.5*(xcmax-x)*wi(i)*al*(cg31c-cg32c)*gluxy
      f3wmn=f3wmn-0.5*(xcmax-x)*wi(i)*al*(cg31c-cg32c)*gluxy

      ff2y=2.*(1.-2.*ycc+2.*ycc*ycc)/xcmax
      ff3y=ff2y*xcmax
      aldl=2*al*dlog(1./epsc)
      if(epsc.gt.1.d0) aldl=0.d0

      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff2y*(-df2wppxy
     .+3.*f2wppxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*(-df2wppxcc+3.*f2wppxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f2wpp=f2wpp-0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff2y*gluxy
     .-2.*gluxcc)/(1.-ycc)

      f3wpp=f3wpp-0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff3y*(-df3wppxy
     .+3.*f3wppxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*xcmax*(-df3wppxcc+3.*f3wppxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f3wpp=f3wpp-0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff3y*gluxy
     .-2.*xcmax*gluxcc)/(1.-ycc)

      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff2y*(-df2wmpxy
     .+3.*f2wmpxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*(-df2wmpxcc+3.*f2wmpxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f2wmp=f2wmp-0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff2y*gluxy
     .-2.*gluxcc)/(1.-ycc)

      f3wmp=f3wmp-0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff3y*(-df3wmpxy
     .+3.*f3wmpxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*xcmax*(-df3wmpxcc+3.*f3wmpxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f3wmp=f3wmp+0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff3y*gluxy
     .-2.*xcmax*gluxcc)/(1.-ycc)


      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff2y*(-df2wppxy
     .+3.*f2wppxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*(-df2wppxcc+3.*f2wppxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f2wpn=f2wpn-0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff2y*gluxy
     .-2.*gluxcc)/(1.-ycc)

      f3wpn=f3wpn-0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff3y*(-df3wppxy
     .+3.*f3wppxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*xcmax*(-df3wppxcc+3.*f3wppxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f3wpn=f3wpn-0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff3y*gluxy
     .-2.*xcmax*gluxcc)/(1.-ycc)

      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff2y*(-df2wmpxy
     .+3.*f2wmpxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*(-df2wmpxcc+3.*f2wmpxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f2wmn=f2wmn-0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff2y*gluxy
     .-2.*gluxcc)/(1.-ycc)

      f3wmn=f3wmn-0.5*(xcmax-x)*wi(i)*1.*epsc*(y*ff3y*(-df3wmpxy
     .+3.*f3wmpxy)*(1.-38.*al*epsc
     .*(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
     .-2.*xcmax*(-df3wmpxcc+3.*f3wmpxcc)*(1.-38.*al*epsc*(1.-
     .0.0*dlog(epsc))*(dlog(4+1/xcc**0.25)-1.6-2*xcc)))/(1.-ycc)
      f3wmn=f3wmn+0.5*(xcmax-x)*wi(i)*aldl*epsc*(y*ff3y*gluxy
     .-2.*xcmax*gluxcc)/(1.-ycc)

      flwpp=flwpp+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wppxy
      flwmp=flwmp+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wmpxy
      flwpn=flwpn+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wpnxy
      flwmn=flwmn+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wmnxy
      flwpp=flwpp+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy
      flwmp=flwmp+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy
      flwpn=flwpn+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy
      flwmn=flwmn+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy

      rrycc=rrz(ycc)
      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*1./xcmax*f2wppxy*(-2.*c0c2l*rrycc)
      f3wpp=f3wpp-0.5*(xcmax-x)*wi(i)*1./xcmax*f3wppxy*(-2.*c0c3l*rrycc)
      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*1./xcmax*f2wmpxy*(-2.*c0c2l*rrycc)
      f3wmp=f3wmp-0.5*(xcmax-x)*wi(i)*1./xcmax*f3wmpxy*(-2.*c0c3l*rrycc)
      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*1./xcmax*f2wppxy*(-2.*c0c2l*rrycc)
      f3wpn=f3wpn-0.5*(xcmax-x)*wi(i)*1./xcmax*f3wppxy*(-2.*c0c3l*rrycc)
      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*1./xcmax*f2wmpxy*(-2.*c0c2l*rrycc)
      f3wmn=f3wmn-0.5*(xcmax-x)*wi(i)*1./xcmax*f3wmpxy*(-2.*c0c3l*rrycc)

      cf21=coeff21(y,epsc)
      cf31=coeff31(y,epsc)
      cf22=coeff22(y,epsc)
      cf32=coeff32(y,epsc)

      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*1.*f2wppxy*(cf21+cf22)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wpp=f3wpp-0.5*(xcmax-x)*wi(i)*1.*f3wppxy*(cf31+cf32)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))
     .*(dlog(4+1/xy**0.25)-1.6-2*xy))


      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*1.*f2wmpxy*(cf21+cf22)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wmp=f3wmp-0.5*(xcmax-x)*wi(i)*1.*f3wmpxy*(cf31+cf32)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))
     .*(dlog(4+1/xy**0.25)-1.6-2*xy))


      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*1.*f2wppxy*(cf21+cf22)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wpn=f3wpn-0.5*(xcmax-x)*wi(i)*1.*f3wppxy*(cf31+cf32)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))
     .*(dlog(4+1/xy**0.25)-1.6-2*xy))


      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*1.*f2wmpxy*(cf21+cf22)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wmn=f3wmn-0.5*(xcmax-x)*wi(i)*1.*f3wmpxy*(cf31+cf32)*
     .(1.-38.*al*epsc*(1.-0.0*dlog(epsc))
     .*(dlog(4+1/xy**0.25)-1.6-2*xy))


      argcc=(1.-ycc)/(1.-xcc)
      fextra=rrycc*dlog(argcc)
      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*1.*f2wppxy*fextra*(-2.)*
     .2.*epsc*(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wpp=f3wpp-0.5*(xcmax-x)*wi(i)*1.*f3wppxy*fextra*(-2.)*
     .2.*epsc*xcmax*(1.-38.*al*epsc*
     .(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*1.*f2wmpxy*fextra*(-2.)*
     .2.*epsc*(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wmp=f3wmp-0.5*(xcmax-x)*wi(i)*1.*f3wmpxy*fextra*(-2.)*
     .2.*epsc*xcmax*(1.-38.*al*epsc*
     .(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*1.*f2wppxy*fextra*(-2.)*
     .2.*epsc*(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wpn=f3wpn-0.5*(xcmax-x)*wi(i)*1.*f3wppxy*fextra*(-2.)*
     .2.*epsc*xcmax*(1.-38.*al*epsc*
     .(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))
      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*1.*f2wmpxy*fextra*(-2.)*
     .2.*epsc*(1.-38.*al*epsc*(1.-0.0*dlog(epsc))*
     .(dlog(4+1/xy**0.25)-1.6-2*xy))
      f3wmn=f3wmn-0.5*(xcmax-x)*wi(i)*1.*f3wmpxy*fextra*(-2.)*
     .2.*epsc*xcmax*(1.-38.*al*epsc*
     .(1.-0.0*dlog(epsc))*(dlog(4+1/xy**0.25)-1.6-2*xy))


  323 CONTINUE
  321 CONTINUE


c  computes contribution proportional to Vcb**2

      vcb=0.047
      vcb2=vcb*vcb
      epsb4=qsct/q2
      epsb=epsb4/4.
      xbb=(1.+epsb)*x
      xbmax=1./(1.+epsb)
      if(xbmax.le.x) go to 421
      if(xbb.gt.0.95) xbb=0.95
      call mrst2002(xbb,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      bf2wpp= 2.*chm
      bf3wpp=-2.*chm
      bf2wmp= 2.*chm
      bf3wmp= 2.*chm
      gluxbb=glu
      if(epsb.gt.1) gluxbb=0.
      f2wppxbb= 2.*bot
      f3wppxbb= 2.*bot
      f2wmpxbb= 2.*bot
      f3wmpxbb=-2.*bot
      del=0.01
      xbbu=xbb*(1.+0.5*del)
      if(xbbu.gt.1.d0) then
      df2wppx=0.
      df3wppx=0.
      df2wmpx=0.
      df3wmpx=0.
      go to 4325
      endif
      call mrst2002(xbbu,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      f2wppxu= 2.*bot
      f3wppxu= 2.*bot
      f2wmpxu= 2.*bot
      f3wmpxu=-2.*bot
      xbbl=xbb*(1.-0.5*del)
      call mrst2002(xbbl,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      f2wppxl= 2.*bot
      f3wppxl= 2.*bot
      f2wmpxl= 2.*bot
      f3wmpxl=-2.*bot
      df2wppxbb=(f2wppxu-f2wppxl)/del
      df3wppxbb=(f3wppxu-f3wppxl)/del
      df2wmpxbb=(f2wmpxu-f2wmpxl)/del
      df3wmpxbb=(f3wmpxu-f3wmpxl)/del
 4325 c0b2l=2.*epsb*xbmax*dlog(1.-xbb)
      c0b3l=c0b2l*xbmax
      bf2wpp=bf2wpp+c0b2l*(-df2wppxbb+3.*f2wppxbb)
      bf2wpp=bf2wpp-al*c0b2l*gluxbb*2.*dlog(1./epsb)     
      bf2wmp=bf2wmp+c0b2l*(-df2wmpxbb+3.*f2wmpxbb)
      bf2wmp=bf2wmp-al*c0b2l*gluxbb*2.*dlog(1./epsb)     
      bf3wpp=bf3wpp+c0b3l*(-df3wppxbb+3.*f3wppxbb)
      bf3wpp=bf3wpp+al*c0b3l*gluxbb*2.*dlog(1./epsb)     
      bf3wmp=bf3wmp+c0b3l*(-df3wmpxbb+3.*f3wmpxbb)
      bf3wmp=bf3wmp-al*c0b3l*gluxbb*2.*dlog(1./epsb)     
      DO 423 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      ybb=(1.+epsb)*y
      call mrst2002(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      gluxy=glu
      f2wppxy= 2.*bot
      f3wppxy= 2.*bot
      f2wmpxy= 2.*bot
      f3wmpxy=-2.*bot
      del=0.01
      xyu=xy*(1.+0.5*del)
      if(xyu.gt.1.d0) then
      df2wppxy=0.
      df3wppxy=0.
      df2wmpxy=0.
      df3wmpxy=0.
      go to 4324
      endif
      call mrst2002(xyu,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      f2wppxyu= 2.*bot
      f3wppxyu= 2.*bot
      f2wmpxyu= 2.*bot
      f3wmpxyu=-2.*bot
      xyl=xy*(1.-0.5*del)
      call mrst2002(xyl,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      f2wppxyl= 2.*bot
      f3wppxyl= 2.*bot
      f2wmpxyl= 2.*bot
      f3wmpxyl=-2.*bot
      df2wppxy=(f2wppxyu-f2wppxyl)/del
      df3wppxy=(f3wppxyu-f3wppxyl)/del
      df2wmpxy=(f2wmpxyu-f2wmpxyl)/del
      df3wmpxy=(f3wmpxyu-f3wmpxyl)/del
 4324 c0b2=cheavy(4,y,epsb)
      if(epsb.gt.1.d0) c0b2=0.d0
      c0b3=cheavy(6,y,epsb)
      if(epsb.gt.1.d0) c0b3=0.d0
      cg21b=2.*cheavy(3,y,epsb)
      cg22b=2.*c0b2*dlog(1./epsb) 
      cg31b=2.*cheavy(5,y,epsb)
      cg32b=2.*c0b3*dlog(1./epsb) 
      bf2wpp=bf2wpp+0.5*(xbmax-x)*wi(i)*c0b2*(-df2wppxy+3.*f2wppxy)
      bf3wpp=bf3wpp+0.5*(xbmax-x)*wi(i)*c0b3*(-df3wppxy+3.*f3wppxy)
      bf2wmp=bf2wmp+0.5*(xbmax-x)*wi(i)*c0b2*(-df2wmpxy+3.*f2wmpxy)
      bf3wmp=bf3wmp+0.5*(xbmax-x)*wi(i)*c0b3*(-df3wmpxy+3.*f3wmpxy)
      bf2wpp=bf2wpp+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b)*gluxy
      bf2wmp=bf2wmp+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b)*gluxy
      bf3wpp=bf3wpp-0.5*(xbmax-x)*wi(i)*al*(cg31b-cg32b)*gluxy
      bf3wmp=bf3wmp+0.5*(xbmax-x)*wi(i)*al*(cg31b-cg32b)*gluxy
      ff2y=2.*(1.-2.*ybb+2.*ybb*ybb)/xbmax
      ff3y=ff2y*xbmax
      aldl=2.*al*dlog(1./epsb)
      if(epsb.gt.1.d0) aldl=0.d0

      bf2wpp=bf2wpp+0.5*(xbmax-x)*wi(i)*epsb*(y*ff2y*(-df2wppxy
     .+3.*f2wppxy)-2.*(-df2wppxbb+3.*f2wppxbb))/(1.-ybb)
      bf2wpp=bf2wpp-0.5*(xbmax-x)*wi(i)*aldl*epsb*(y*ff2y*gluxy
     .-2.*gluxbb)/(1.-ybb)
      bf3wpp=bf3wpp+0.5*(xbmax-x)*wi(i)*epsb*(y*ff3y*(-df3wppxy
     .+3.*f3wppxy)-2.*xbmax*(-df3wppxbb+3.*f3wppxbb))/(1.-ybb)
      bf3wpp=bf3wpp+0.5*(xbmax-x)*wi(i)*aldl*epsb*(y*ff3y*gluxy
     .-2.*xbmax*gluxbb)/(1.-ybb)

      bf2wmp=bf2wmp+0.5*(xbmax-x)*wi(i)*epsb*(y*ff2y*(-df2wmpxy
     .+3.*f2wmpxy)-2.*(-df2wmpxbb+3.*f2wmpxbb))/(1.-ybb)
      bf2wmp=bf2wmp-0.5*(xbmax-x)*wi(i)*aldl*epsb*(y*ff2y*gluxy
     .-2.*gluxbb)/(1.-ybb)
      bf3wmp=bf3wmp+0.5*(xbmax-x)*wi(i)*epsb*(y*ff3y*(-df3wmpxy
     .+3.*f3wmpxy)-2.*xbmax*(-df3wmpxbb+3.*f3wmpxbb))/(1.-ybb)
      bf3wmp=bf3wmp-0.5*(xbmax-x)*wi(i)*aldl*epsb*(y*ff3y*gluxy
     .-2.*xbmax*gluxbb)/(1.-ybb)
 
      rrybb=rrz(ybb)
      bf2wpp=bf2wpp+0.5*(xbmax-x)*wi(i)*f2wppxy*(-2.*c0b2l*rrybb)/xbmax
      bf3wpp=bf3wpp+0.5*(xbmax-x)*wi(i)*f3wppxy*(-2.*c0b3l*rrybb)/xbmax
      bf2wmp=bf2wmp+0.5*(xbmax-x)*wi(i)*f2wmpxy*(-2.*c0b2l*rrybb)/xbmax
      bf3wmp=bf3wmp+0.5*(xbmax-x)*wi(i)*f3wmpxy*(-2.*c0b3l*rrybb)/xbmax

      cf21=coeff21(y,epsb)
      cf31=coeff31(y,epsb)
      cf22=coeff22(y,epsb)
      cf32=coeff32(y,epsb)

      bf2wpp=bf2wpp+0.5*(xbmax-x)*wi(i)*f2wppxy*(cf21+cf22)
      bf3wpp=bf3wpp+0.5*(xbmax-x)*wi(i)*f3wppxy*(cf31+cf32)
      bf2wmp=bf2wmp+0.5*(xbmax-x)*wi(i)*f2wmpxy*(cf21+cf22)
      bf3wmp=bf3wmp+0.5*(xbmax-x)*wi(i)*f3wmpxy*(cf31+cf32)

      argbb=(1.-ybb)/(1.-xbb)
      fextra=rrybb*dlog(argbb)
      bf2wpp=bf2wpp+0.5*(xbmax-x)*wi(i)*f2wppxy*fextra*(-2.)*
     .2.*epsb
      bf3wpp=bf3wpp+0.5*(xbmax-x)*wi(i)*f3wppxy*fextra*(-2.)*
     .2.*epsb*xbmax
      bf2wmp=bf2wmp+0.5*(xbmax-x)*wi(i)*f2wmpxy*fextra*(-2.)*
     .2.*epsb
      bf3wmp=bf3wmp+0.5*(xbmax-x)*wi(i)*f3wmpxy*fextra*(-2.)*
     .2.*epsb*xbmax

  423 CONTINUE
  421 CONTINUE

      f2wpp=f2wpp+vcb2*bf2wpp
      f3wpp=f3wpp+vcb2*bf3wpp
      f2wmp=f2wmp+vcb2*bf2wmp
      f3wmp=f3wmp+vcb2*bf3wmp
      f2wpn=f2wpn+vcb2*bf2wpp
      f3wpn=f3wpn+vcb2*bf3wpp
      f2wmn=f2wmn+vcb2*bf2wmp
      f3wmn=f3wmn+vcb2*bf3wmp

   27 RETURN
      END

      real*8 function rrz(z)
      implicit real*8(a-h,o-z)
      r7=dsqrt(7.d0)
      argz=-0.5*r7*dlog(z)
      rrz=dsqrt(z)*(dcos(argz)+(3./r7)*dsin(argz))
      return
      end

      real*8 function coeff21(y,eps)
      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      r7=dsqrt(7.d0)
      xcmax=1./(1.+eps)
      ycc=y/xcmax
      fcextra=0.
      do 1323 ii=1,nterms
      z=0.5*(xcmax-y)*xi(ii)+0.5*(xcmax+y)
      zcc=z/xcmax
      yz=y/z
      argz=-0.5*r7*dlog(yz)
      arg0=-0.5*r7*dlog(ycc)
      rz=dsqrt(yz)*(dcos(argz)+(3./r7)*dsin(argz))
      r0=dsqrt(ycc)*(dcos(arg0)+(3./r7)*dsin(arg0))
      ff2z=2.*(1.-2.*zcc*(1.-zcc))
      term=eps*(ff2z*rz-2.*r0)*(1.+eps)
      fcextra=fcextra+0.5*(xcmax-y)*wi(ii)*(-2.*term)/(1.-zcc)
 1323 continue
      coeff21=fcextra
      return
      end

      real*8 function coeff31(y,eps)
      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      r7=dsqrt(7.d0)
      xcmax=1./(1.+eps)
      ycc=y/xcmax
      fcextra=0.
      do 1323 ii=1,nterms
      z=0.5*(xcmax-y)*xi(ii)+0.5*(xcmax+y)
      zcc=z/xcmax
      yz=y/z
      argz=-0.5*r7*dlog(yz)
      arg0=-0.5*r7*dlog(ycc)
      rz=dsqrt(yz)*(dcos(argz)+(3./r7)*dsin(argz))
      r0=dsqrt(ycc)*(dcos(arg0)+(3./r7)*dsin(arg0))
      ff3z=2.*(1.-2.*zcc*(1.-zcc))
      term=eps*(ff3z*rz-2.*r0)
      fcextra=fcextra+0.5*(xcmax-y)*wi(ii)*(-2.*term)/(1.-zcc)
 1323 continue
      coeff31=fcextra
      return
      end

      real*8 function coeff22(y,eps)
      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      r7=dsqrt(7.d0)
      xcmax=1./(1.+eps)
      ycc=y/xcmax
      fcextra=0.
      do 1323 ii=1,nterms
      z=0.5*(xcmax-y)*xi(ii)+0.5*(xcmax+y)
      yz=y/z
      c0c2=cheavy(4,z,eps)
      if(eps.gt.1.d0) c0c2=0.
      argz=-0.5*r7*dlog(yz)
      rz=dsqrt(yz)*(dcos(argz)+(3./r7)*dsin(argz))
      fcextra=fcextra+0.5*(xcmax-y)*wi(ii)*c0c2*(-2.*rz)/z
 1323 continue
      coeff22=fcextra
      return
      end

      real*8 function coeff32(y,eps)
      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      r7=dsqrt(7.d0)
      xcmax=1./(1.+eps)
      ycc=y/xcmax
      fcextra=0.
      do 1323 ii=1,nterms
      z=0.5*(xcmax-y)*xi(ii)+0.5*(xcmax+y)
      yz=y/z
      c0c3=cheavy(6,z,eps)
      if(eps.gt.1.d0) c0c3=0.
      argz=-0.5*r7*dlog(yz)
      rz=dsqrt(yz)*(dcos(argz)+(3./r7)*dsin(argz))
      fcextra=fcextra+0.5*(xcmax-y)*wi(ii)*c0c3*(-2.*rz)/z
 1323 continue
      coeff32=fcextra
      return
      end

      DOUBLE PRECISION FUNCTION	ALPHA(T)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/INPUT/alambda,flavor,qsct,qsdt,iord
      DATA PI/3.14159/
      DATA TOL/.0005/
      ITH=0
      TT=T
      qsdtt=qsdt/4.
      qsctt=qsct/4.
      AL=ALAMBDA
      AL2=AL*AL
      FLAV=4.
      QS=AL2*dEXP(T)

      if(qs.lt.0.5d0) then   !!  running stops below 0.5
          qs=0.5d0
          t=dlog(qs/al2)
          tt=t
      endif

      IF(QS.gt.QSCTT) GO	TO 12  
      IF(QS.lt.QSDTT) GO	TO 312  
   11 CONTINUE
      B0=11-2.*FLAV/3. 
      IF(IORD)1,1,2
c     IF(IORD)2,2,2	!TAKE CARE !!
    1 CONTINUE
      ALPHA=4.*PI/B0/T
      RETURN
    2 CONTINUE
      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*dLOG(T)/T)
    5 CONTINUE
      F=-T+X1/AS-X2*dLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF(DEL-TOL)3,3,4
    3 CONTINUE
      ALPHA=AS2
      IF(ITH.EQ.0) RETURN
      GO TO (13,14,15) ITH
    4 CONTINUE
      AS=AS2
      GO TO 5
   12 ITH=1
      T=dLOG(QSCTT/AL2)
      GO TO 11
   13 ALFQC4=ALPHA
      FLAV=5.   
      ITH=2
      GO TO 11
   14 ALFQC5=ALPHA
      ITH=3
      T=TT
      GO TO 11
   15 ALFQS5=ALPHA
      ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
      ALPHA=1./ALFINV
      RETURN

  311 CONTINUE
      B0=11-2.*FLAV/3. 
      IF(IORD)31,31,32
c     IF(IORD)32,32,32	!TAKE CARE !!
   31 CONTINUE
      ALPHA=4.*PI/B0/T
      RETURN
   32 CONTINUE
      X1=4.*PI/B0
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS=X1/T*(1.-X2*dLOG(T)/T)
   35 CONTINUE
      F=-T+X1/AS-X2*dLOG(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF(DEL-TOL)33,33,34
   33 CONTINUE
      ALPHA=AS2
      IF(ITH.EQ.0) RETURN
      GO TO (313,314,315) ITH
   34 CONTINUE
      AS=AS2
      GO TO 35
  312 ITH=1
      T=dLOG(QSDTT/AL2)
      GO TO 311
  313 ALFQC4=ALPHA
      FLAV=3.   
      ITH=2
      GO TO 311
  314 ALFQC3=ALPHA
      ITH=3
      T=TT
      GO TO 311
  315 ALFQS3=ALPHA
      ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
      ALPHA=1./ALFINV
      RETURN
      END


      SUBROUTINE WATE96
C*******************************************************************
C*****							       *****
C***** THE X(I)	AND W(I) ARE THE DIRECT	OUTPUT FROM A PROGRAM  *****
C***** USING NAG ROUTINE D01BCF	TO CALCULATE THE	       *****
C***** GAUSS-LEGENDRE WEIGHTS FOR 96 POINT INTEGRATION.	       *****
C***** THEY AGREE TO TYPICALLY 14 DECIMAL PLACES WITH THE      *****
C***** TABLE IN	ABRAMOWITZ & STEGUN, PAGE 919.		       *****
C*****							       *****
C***** ---->   PETER HARRIMAN, APRIL 3RD 1990.		       *****
C*****							       *****
C*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION	X(48),W(48)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      NTERMS=96

      X( 1)=   0.01627674484960183561
      X( 2)=   0.04881298513604856015
      X( 3)=   0.08129749546442434360
      X( 4)=   0.11369585011066471632
      X( 5)=   0.14597371465489567682
      X( 6)=   0.17809688236761733390
      X( 7)=   0.21003131046056591064
      X( 8)=   0.24174315616383866556
      X( 9)=   0.27319881259104774468
      X(10)=   0.30436494435449495954
      X(11)=   0.33520852289262397655
      X(12)=   0.36569686147231213885
      X(13)=   0.39579764982890709712
      X(14)=   0.42547898840729897474
      X(15)=   0.45470942216774136446
      X(16)=   0.48345797392059470382
      X(17)=   0.51169417715466604391
      X(18)=   0.53938810832435567233
      X(19)=   0.56651041856139533470
      X(20)=   0.59303236477757022282
      X(21)=   0.61892584012546672523
      X(22)=   0.64416340378496526886
      X(23)=   0.66871831004391424358
      X(24)=   0.69256453664216964528
      X(25)=   0.71567681234896561582
      X(26)=   0.73803064374439816819
      X(27)=   0.75960234117664555964
      X(28)=   0.78036904386743123629
      X(29)=   0.80030874413913884180
      X(30)=   0.81940031073792957139
      X(31)=   0.83762351122818502758
      X(32)=   0.85495903343459936363
      X(33)=   0.87138850590929436968
      X(34)=   0.88689451740241818933
      X(35)=   0.90146063531585023110
      X(36)=   0.91507142312089592706
      X(37)=   0.92771245672230655266
      X(38)=   0.93937033975275308073
      X(39)=   0.95003271778443564022
      X(40)=   0.95968829144874048809
      X(41)=   0.96832682846326217918
      X(42)=   0.97593917458513455843
      X(43)=   0.98251726356301274934
      X(44)=   0.98805412632962202890
      X(45)=   0.99254390032376081654
      X(46)=   0.99598184298720747465
      X(47)=   0.99836437586317963722
      X(48)=   0.99968950388322870559
      W( 1)=   0.03255061449236316962
      W( 2)=   0.03251611871386883307
      W( 3)=   0.03244716371406427668
      W( 4)=   0.03234382256857594104
      W( 5)=   0.03220620479403026124
      W( 6)=   0.03203445623199267876
      W( 7)=   0.03182875889441101874
      W( 8)=   0.03158933077072719007
      W( 9)=   0.03131642559686137819
      W(10)=   0.03101033258631386231
      W(11)=   0.03067137612366917839
      W(12)=   0.03029991542082762553
      W(13)=   0.02989634413632842385
      W(14)=   0.02946108995816795100
      W(15)=   0.02899461415055528410
      W(16)=   0.02849741106508543861
      W(17)=   0.02797000761684838950
      W(18)=   0.02741296272602931385
      W(19)=   0.02682686672559184485
      W(20)=   0.02621234073567250055
      W(21)=   0.02557003600534944960
      W(22)=   0.02490063322248370695
      W(23)=   0.02420484179236479915
      W(24)=   0.02348339908592633665
      W(25)=   0.02273706965832950717
      W(26)=   0.02196664443874448477
      W(27)=   0.02117293989219144572
      W(28)=   0.02035679715433347898
      W(29)=   0.01951908114014518992
      W(30)=   0.01866067962741165898
      W(31)=   0.01778250231604547316
      W(32)=   0.01688547986424539715
      W(33)=   0.01597056290256253144
      W(34)=   0.01503872102699521608
      W(35)=   0.01409094177231515264
      W(36)=   0.01312822956696188190
      W(37)=   0.01215160467108866759
      W(38)=   0.01116210209983888144
      W(39)=   0.01016077053500880978
      W(40)=   0.00914867123078384552
      W(41)=   0.00812687692569928101
      W(42)=   0.00709647079115442616
      W(43)=   0.00605854550423662775
      W(44)=   0.00501420274292825661
      W(45)=   0.00396455433844564804
      W(46)=   0.00291073181793626202
      W(47)=   0.00185396078894924657
      W(48)=   0.00079679206555731759
      DO 1 I=1,48
      XI(I)=-X(49-I)
      WI(I)=W(49-I)
      XI(I+48)=X(I)
      WI(I+48)=W(I)
    1 CONTINUE
      DO 2 I=1,96
    2 XX(I)=0.5*(XI(I)+1.)
      XX(97)=1.0
      EXPON=1.0
      DO 3 I=1,96
      YI=2.*(0.5*(1.+XI(I)))**EXPON-1.
      WI(I)=WI(I)/(1.+XI(I))*(1.+YI)*EXPON
      XI(I)=YI
      XX(I)=0.5*(1.+YI)
    3 CONTINUE
      RETURN
      END

      real*8 function cheavy(i,z,eps)
      implicit real*8(a-h,o-z)

c     this function returns the values of C_g(z,Q^2) 
c     and the deriv. wrt log Q^2. Here eps=m^2/Q^2.
c     If i=1,3  C_g for F2.  If i=2,4 deriv of C_g for F2 
c     i=1,2 refer to photon current
c     i=3,4 refer to W current for F2
c     i=5,6 refer to W current for xF3

      if(i.gt.6) stop
      z1=1.-z
      z2=z*z
      z3=z2*z
      zr=z1/z
      eps1=1.+eps
      eps2=eps*eps
      beta2=1.-4.*eps*z/z1
      if(i.gt.2) beta2=1.-eps*z/z1
      if(beta2.lt.0.) go to 10
      beta=dsqrt(beta2)
      a=z2+z1*z1
      b=4.*z*(1.-3.*z)
      c=-8.*z2
      aa=8.*z*z1-1.
      bb=-4.*z*z1
      arg=(1.+beta)/(1.-beta)
      fac=dlog(arg)
      alam=1./(1.+eps)
      zp=z/alam
      alz=alam*z
      alz1=1.-alz
      al1=1.-alam
      al2=alam*alam
      al3=al2*alam
      al4=al2*al2
      z1=1.-z
      z2=z*z
      zp1=1.-zp
      go to (1,2,3,4,5,6) i
    1 cheavy=(a+b*eps+c*eps2)*fac+(aa+bb*eps)*beta
      return
    2 cheavy=(-b*eps-2.*c*eps2)*fac+(a+b*eps+c*eps2)/beta
     .      +(-bb*eps)*beta +(aa+bb*eps)*2.*z*eps/z1/beta
      return
    3 pqg=0.5*(zp*zp+zp1*zp1)
      arg0=alam*zp1*zp1/z1/z
      c0=pqg*dlog(arg0)
      arg=alam*z1/al1/z
      bigl=dlog(arg)
      f1=c0+(8.-18.*al1+12.*al1*al1)*zp*zp1+(al1/z1-1.)
      f3=al1*z*bigl*6.*(1.-2.*z)
      f2=pqg*( bigl-dlog(alam))
      cheavy=2.*(f1+f3+f2)/alam
      return
    4 arg=alam/al1/z2
      argbig=alam*z1/al1/z
      bigl=dlog(argbig)
      f1=(1.-2.*zp+2.*zp*zp)/al2
      f2=(4.*z*alam-6.*z2-al2)/al4*dlog(alam/z1/z)
      f3=4.*z/al4*(3.*z-2.*(1.+3.*z)*alam+3.*(1.+2.*z)*al2)-2.*z/al2/z1
      f4=12.*z*(1.-2.*z)/al2*(1.-bigl)
      f5=(1.-2.*zp+2.*zp*zp)/alam/al1
     .+(4.*z*alam-6.*z2-al2)/al4*dlog(z1/al1/z)
      f6=2./al4*(4.*z*alam-6.*z2-al2)*dlog(zp1)
      cheavy=(f1+f2+f3+f4+f5+f6)*eps*al2
      return
    5 alam=1./(1.+eps)
      zp=z/alam
      alz=alam*z
      alz1=1.-alz
      al1=1.-alam
      z1=1.-z
      zp1=1.-zp
      pqg=0.5*(zp*zp+zp1*zp1)
      arg0=alam*zp1*zp1/z1/z
      c0=pqg*dlog(arg0)
      arg=alam*z1/al1/z
      bigl=dlog(arg)
      f1=c0+2.*al1*zp*zp1+al1*zp*bigl*(-2.*zp1+2.*z)
      f2=pqg*(-bigl-dlog(alam))
      cheavy=2.*(f1+f2)
      return
    6 alam=1./(1.+eps)
      zp=z/alam
      alz=alam*z
      alz1=1.-alz
      al1=1.-alam
      al2=alam*alam
      al3=al2*alam
      z1=1.-z
      zp1=1.-zp
      arg=alam/al1/z2
      f1=(1.-2.*zp+2.*zp*zp)*(-1./alam/al1)
      f2=2.*z/al3*(alam-2.*z)*dlog(arg)
      f3=4.*z/al3*(3.*z-2.*alam)
      f4=4.*z/al3*(alam-2.*z)*dlog(zp1)
      cheavy=(f1+f2+f3+f4)*eps*al2
      return
   10 print 99
   99 format(1x,'x > x0')
      print 98,i,z,eps,beta2
   98 format(1x,i3,' z=',f10.6,' eps=',f10.6,' beta2=',f10.6)
      stop
      end

      real*8 function clheavy(i,z,eps)
      implicit real*8(a-h,o-z)

c     this function returns the values of C_g(z,Q^2) 
c     and the deriv. wrt log Q^2. Here eps=m^2/Q^2.
c     If i=1  C_g for F2.  If i=2 deriv of C_g for F2 
c     If i=3  C_g for FL.  If i=4 (1-m2/Q2)*beta*Clq(massless)
      if(i.gt.4) stop
      z1=1.-z
      z2=z*z
      eps2=eps*eps
      beta2=1.-eps*z/z1
      if(beta2.lt.0.) go to 10
      beta=dsqrt(beta2)
      a=z2+z1*z1
      b=4.*z*(1.-3.*z)
      c=-8.*z2
      aa=8.*z*z1-1.
      bb=-4.*z*z1
      arg=(1.+beta)/(1.-beta)
      fac=dlog(arg)
      cf=4./3.
      ca=3.
      enf=4.
      ZETA2=1.64493406684823
      ZETA3=1.20205690315959
      go to (1,2,3,4) i
    1 clheavy=(a+b*eps+c*eps2)*fac+(aa+bb*eps)*beta
      return
    2 clheavy=(-b*eps-2.*c*eps2)*fac+(a+b*eps+c*eps2)/beta
     .      +(-bb*eps)*beta +(aa+bb*eps)*2.*z*eps/z1/beta
      return
    3 clheavy=-bb*beta+c*eps*fac
      return
    4 clheavy=(1.-eps)*beta*4.*cf*z
      return
   10 print 99,i,y,eps
   99 format(1x,'x > x0',i3,2(1x,e9.3))
      stop
      end


      real*8 FUNCTION RATFE(X)
      implicit real*8(a-h,o-z)
C-----THIS IS MY LOGARITHMICALLY LIN-CON-LIN VERSION FROM EMC SEPT89
 
      IF(X.LT.0.0903) THEN
        RATFE=1.238+0.203*dLOG10(X)
      ELSEIF(X.GT.0.234) THEN
        RATFE=0.783-0.385*dLOG10(X)
      ELSE
        RATFE=1.026
      ENDIF
 
      RETURN
      END


