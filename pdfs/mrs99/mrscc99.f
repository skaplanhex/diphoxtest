***********************************************************************
*                                                                     *
*    Program for generating charged Current Structure Functions using *
*    consistent treatment of charm and bottom structure functions     *
*    Charm mass = 1.43 GeV   Bottom mass = 4.3 GeV                    *
*    A.D. Martin, R.G. Roberts, W.J. Stirling and R.S. Thorne         *
*    University of Durham preprint DTP/99/64                          *
*    RAL preprint    RAL-TR-1999-047                                  *
*    Setting MODE=1 gives MRST99      L(4)=300 MeV a_s(M_Z)=0.1175    *
*    Setting MODE=2 gives MRST99(g^^) L(4)=300 MeV a_s(M_Z)=0.1175    *
*    Setting MODE=3 gives MRST99(gvv) L(4)=300 MeV a_s(M_Z)=0.1175    *
*    Setting MODE=4 gives MRST99(avv) L(4)=229 MeV a_s(M_Z)=0.1125    *
*    Setting MODE=5 gives MRST99(a^^) L(4)=383 MeV a_s(M_Z)=0.1225    *
*                                                                     *
*    The program should be run only with iord set to 1                *
*                                                                     *
*    The program evaluates the structure functions F2, xF3 and        *
*    FL (only to Order alpha_S)                                       *
*    for W+ p , W- p , W+ n and W- n scattering                       *
*                                                                     *
***********************************************************************

      implicit real*8(a-h,o-z)
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS
      COMMON/INPUT/alambda,flavor,qsct,qsdt,iord

      call wate96
      iord=1
      flavor=3
      qsdt=8.18
      qsct=74.
   99 read(5,*) mode,x,q2

      if(mode.eq.1) alambda=0.300 
      if(mode.eq.2) alambda=0.300 
      if(mode.eq.3) alambda=0.300 
      if(mode.eq.4) alambda=0.229 
      if(mode.eq.5) alambda=0.383

      call sfun(x,q2,mode,f2wpp,f3wpp,f2wmp,f3wmp,f2wpn,f3wpn,
     .f2wmn,f3wmn,flwpp,flwmp,flwpn,flwmn)      

      print 100,x,q2,f2wpp,f3wpp,f2wmp,f3wmp,flwpp,flwmp,
     .f2wpn,f3wpn,f2wmn,f3wmn,flwpn,flwmn

      go to 99

  100 format(1x,'x= ',f10.6,' Q2= ',f14.3//
     11x,' Proton Structure Functions '/
     25x,' F2(W+p)= ',f8.4,' xF3(W+p)= ',f8.4,
     3' F2(W-p)= ',f8.4,' xF3(W-p)= ',f8.4/
     25x,' FL(W+p)= ',f8.4,' FL(W-p)= ',f8.4//
     41x,' Neutron Structure Functions '/
     55x,' F2(W+n)= ',f8.4,' xF3(W+n)= ',f8.4,
     6' F2(W-n)= ',f8.4,' xF3(W-n)= ',f8.4/
     25x,' FL(W+n)= ',f8.4,' FL(W-n)= ',f8.4)

      end


      subroutine sfun(x,q2,mode,f2wpp,f3wpp,f2wmp,f3wmp,
     .f2wpn,f3wpn,f2wmn,f3wmn,flwpp,flwmp,flwpn,flwmn)

CC  NOTATION

cc  f2wpp is the  F2 str. fn. for W+ proton scattering  
cc  f3wpp is the xF3 str. fn. for W+ proton scattering  
cc  f2wmp is the  F2 str. fn. for W- proton scattering  
cc  f3wpp is the xF3 str. fn. for W- proton scattering  
cc  flwpp is the  FL str. fn. for W+ proton scattering  
cc  flwmp is the  FL str. fn. for W- proton scattering  
cc  f2wpn is the  F2 str. fn. for W+ neutron scattering  
cc  f3wpn is the xF3 str. fn. for W+ neutron scattering  
cc  f2wmn is the  F2 str. fn. for W- neutron scattering  
cc  f3wpn is the xF3 str. fn. for W- neutron scattering
cc  flwpn is the  FL str. fn. for W+ neutron scattering  
cc  flwmn is the  FL str. fn. for W- neutron scattering  

cc  In computing the longitudinal str. fn., only the order alpha_s
cc  terms are used.
   
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

      call mrs99(x,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      f2wppx=2.*((dnv+dsea)+str+usea)
      f3wppx=2.*((dnv+dsea)+str-usea)
      f2wmpx=2.*( dsea+str+(upv+usea))
      f3wmpx=2.*(-dsea-str+(upv+usea))
      f2wppx0=2.*(cc2*(dnv+dsea)+sc2*str+usea)
      f3wppx0=2.*(cc2*(dnv+dsea)+sc2*str-usea)
      f2wmpx0=2.*( cc2*dsea+sc2*str+(upv+usea))
      f3wmpx0=2.*(-cc2*dsea-sc2*str+(upv+usea))
      f2wpnx=2.*((upv+usea)+str+dsea)
      f3wpnx=2.*((upv+usea)+str-dsea)
      f2wmnx=2.*( usea+str+(dnv+dsea))
      f3wmnx=2.*(-usea-str+(dnv+dsea))
      f2wpnx0=2.*(cc2*(upv+usea)+sc2*str+dsea)
      f3wpnx0=2.*(cc2*(upv+usea)+sc2*str-dsea)
      f2wmnx0=2.*( cc2*usea+sc2*str+(dnv+dsea))
      f3wmnx0=2.*(-cc2*usea-sc2*str+(dnv+dsea))

      call mrs99(xcc,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      f2wppx0=f2wppx0+2.*(sc2*(dnv+dsea)+cc2*str)
      f3wppx0=f3wppx0+2.*(sc2*(dnv+dsea)+cc2*str)
      f2wmpx0=f2wmpx0+2.*( sc2*dsea+cc2*str)
      f3wmpx0=f3wmpx0+2.*(-sc2*dsea-cc2*str)
      f2wpnx0=f2wpnx0+2.*(sc2*(upv+usea)+cc2*str)
      f3wpnx0=f3wpnx0+2.*(sc2*(upv+usea)+cc2*str)
      f2wmnx0=f2wmnx0+2.*( sc2*usea+cc2*str)
      f3wmnx0=f3wmnx0+2.*(-sc2*usea-cc2*str)

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
      call mrs99(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      f2wppxy=2.*((dnv+dsea)+str+usea)
      f3wppxy=2.*((dnv+dsea)+str-usea)
      f2wmpxy=2.*( dsea+str+(upv+usea))
      f3wmpxy=2.*(-dsea-str+(upv+usea))
      f2wpnxy=2.*((upv+usea)+str+dsea)
      f3wpnxy=2.*((upv+usea)+str-dsea)
      f2wmnxy=2.*( usea+str+(dnv+dsea))
      f3wmnxy=2.*(-usea-str+(dnv+dsea))
      gluxy=glu
      C222=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*dLOG(Y)-2.*(1.+Y)*AL1)
      C223=C222-2.*(1.+Y)
      C23=CF*(-3.+4.*AL1)/(1.-Y)
      CG2=6.*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*dLOG(1./Y-1.))

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
      DO 323 I=1,NTERMS
      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
      XY=X/Y
      call mrs99(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
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
      call mrs99(xyu,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      f2wppxyu= 2.*chm
      f3wppxyu=-2.*chm
      f2wmpxyu= 2.*chm
      f3wmpxyu= 2.*chm
      xyl=xy*(1.-0.5*del)
      call mrs99(xyl,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      f2wppxyl= 2.*chm
      f3wppxyl=-2.*chm
      f2wmpxyl= 2.*chm
      f3wmpxyl= 2.*chm
      df2wppxy=(f2wppxyu-f2wppxyl)/del
      df3wppxy=(f3wppxyu-f3wppxyl)/del
      df2wmpxy=(f2wmpxyu-f2wmpxyl)/del
      df3wmpxy=(f3wmpxyu-f3wmpxyl)/del
 1324 c0c=cheavy(4,y,epsc)
      if(epsc.gt.1.d0) c0c=0.d0
      cg21c=2.*cheavy(3,y,epsc)
      cg22c=2.*c0c*dlog(1./epsc) 
      clg2c=2.*clheavy(3,y,epsc)
      f1lq=clheavy(4,y,epsc)
      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df2wppxy+3.*f2wppxy)
      f3wpp=f3wpp+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df3wppxy+3.*f3wppxy)
      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df2wmpxy+3.*f2wmpxy)
      f3wmp=f3wmp+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df3wmpxy+3.*f3wmpxy)
      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df2wppxy+3.*f2wppxy)
      f3wpn=f3wpn+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df3wppxy+3.*f3wppxy)
      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df2wmpxy+3.*f2wmpxy)
      f3wmn=f3wmn+0.5*(xcmax-x)*wi(i)*2.*c0c*(-df3wmpxy+3.*f3wmpxy)
      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c)*gluxy
      flwpp=flwpp+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wppxy
      flwmp=flwmp+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wmpxy
      flwpn=flwpn+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wpnxy
      flwmn=flwmn+0.5*(xcmax-x)*wi(i)*al*f1lq*f2wmnxy
      flwpp=flwpp+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy
      flwmp=flwmp+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy
      flwpn=flwpn+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy
      flwmn=flwmn+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy
      fextra=coeff3(y,epsc)
      f2wpp=f2wpp+0.5*(xcmax-x)*wi(i)*2.*f2wppxy*fextra
      f3wpp=f3wpp+0.5*(xcmax-x)*wi(i)*2.*f3wppxy*fextra
      f2wmp=f2wmp+0.5*(xcmax-x)*wi(i)*2.*f2wmpxy*fextra
      f3wmp=f3wmp+0.5*(xcmax-x)*wi(i)*2.*f3wmpxy*fextra
      f2wpn=f2wpn+0.5*(xcmax-x)*wi(i)*2.*f2wppxy*fextra
      f3wpn=f3wpn+0.5*(xcmax-x)*wi(i)*2.*f3wppxy*fextra
      f2wmn=f2wmn+0.5*(xcmax-x)*wi(i)*2.*f2wmpxy*fextra
      f3wmn=f3wmn+0.5*(xcmax-x)*wi(i)*2.*f3wmpxy*fextra

  323 CONTINUE
  321 CONTINUE


c  computes contribution proportional to Vcb**2

      vcb=0.047
      vcb2=vcb*vcb
      epsb4=qsct/q2
      epsb=epsb4/4.
      xbb=(1.+epsb)*x
      if(xbb.gt.0.95) xbb=0.95
      call mrs99(xbb,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsc.gt.1.) chm=0.
      b2wpp= 2.*chm
      b3wpp=-2.*chm
      b2wmp= 2.*chm
      b3wmp= 2.*chm
      xbmax=1./(1.+epsb)
      if(xbmax.le.x) go to 421
      DO 423 I=1,NTERMS
      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
      XY=X/Y
      call mrs99(xy,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      gluxy=glu
      b2wppxy= 2.*bot
      b3wppxy= 2.*bot
      b2wmpxy= 2.*bot
      b3wmpxy=-2.*bot
      del=0.01
      xyu=xy*(1.+0.5*del)
      if(xyu.gt.1.d0) then
      db2wppxy=0.
      db3wppxy=0.
      db2wmpxy=0.
      db3wmpxy=0.
      go to 4324
      endif
      call mrs99(xyu,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      b2wppxyu= 2.*bot
      b3wppxyu= 2.*bot
      b2wmpxyu= 2.*bot
      b3wmpxyu=-2.*bot
      xyl=xy*(1.-0.5*del)
      call mrs99(xyl,scale,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      if(epsb.gt.1.) bot=0.
      b2wppxyl= 2.*bot
      b3wppxyl= 2.*bot
      b2wmpxyl= 2.*bot
      b3wmpxyl=-2.*bot
      db2wppxy=(b2wppxyu-b2wppxyl)/del
      db3wppxy=(b3wppxyu-b3wppxyl)/del
      db2wmpxy=(b2wmpxyu-b2wmpxyl)/del
      db3wmpxy=(b3wmpxyu-b3wmpxyl)/del
 4324 c0b=cheavy(4,y,epsb)
      if(epsb.gt.1.d0) c0b=0.d0
      cg21b=2.*cheavy(3,y,epsb)
      cg22b=2.*c0b*dlog(1./epsb) 
      b2wpp=b2wpp+0.5*(xbmax-x)*wi(i)*2.*c0b*(-db2wppxy+3.*b2wppxy)
      b3wpp=b3wpp+0.5*(xbmax-x)*wi(i)*2.*c0b*(-db3wppxy+3.*b3wppxy)
      b2wmp=b2wmp+0.5*(xbmax-x)*wi(i)*2.*c0b*(-db2wmpxy+3.*b2wmpxy)
      b3wmp=b3wmp+0.5*(xbmax-x)*wi(i)*2.*c0b*(-db3wmpxy+3.*b3wmpxy)
      b2wpp=b2wpp+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b)*gluxy
      b2wmp=b2wmp+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b)*gluxy
      b2wpp=b2wpp+0.5*(xbmax-x)*wi(i)*2.*b2wppxy*fextra
      b3wpp=b3wpp+0.5*(xbmax-x)*wi(i)*2.*b3wppxy*fextra
      b2wmp=b2wmp+0.5*(xbmax-x)*wi(i)*2.*b2wmpxy*fextra
      b3wmp=b3wmp+0.5*(xbmax-x)*wi(i)*2.*b3wmpxy*fextra

  423 CONTINUE
  421 CONTINUE

      f2wpp=f2wpp+vcb2*b2wpp
      f3wpp=f3wpp+vcb2*b3wpp
      f2wmp=f2wmp+vcb2*b2wmp
      f3wmp=f3wmp+vcb2*b3wmp
      f2wpn=f2wpn+vcb2*b2wpp
      f3wpn=f3wpn+vcb2*b3wpp
      f2wmn=f2wmn+vcb2*b2wmp
      f3wmn=f3wmn+vcb2*b3wmp

   27 RETURN
      END



      real*8 function coeff3(y,eps)
      implicit real*8(a-h,o-z)
      dimension z1(17),z2(17)
      data z1/
     .-0.216E+00,0.163E+02,-0.199E+02,0.117E+02,-0.117E+02,-0.804E+02,
     .0.607E+02,-0.368E+01,0.365E+02, 0.945E+02, 0.112E+03,-0.177E+03,
     .-0.431E+02,-0.814E+02,-0.132E+03, 0.201E+03, 0.225E+01/
      data z2/
     .-0.935E+00,0.169E+02,-0.574E+01,0.585E+01,-0.919E+02,-0.600E+03,
     .0.357E+03,-0.900E+03,0.276E+04,0.879E+05,-0.349E+05, 0.885E+05,
     .-0.151E+06,-0.117E+05,-0.131E+07,-0.129E+06, 0.582E+02/
      epsl=eps
      eps4=4.d0*eps
      x0=1.d0/(1.d0+eps)
      yx0=y/x0
      yx01=1.d0-yx0
      epsl2=epsl*epsl
      epsl3=epsl2*epsl
      if(y.lt.0.05) go to 100    
      a0=z1(1)+z1(2)*epsl+z1(3)*epsl2+z1(4)*epsl3
      a1=z1(5)+z1(6)*epsl+z1(7)*epsl2+z1(8)*epsl3
      a2=z1(9)+z1(10)*epsl+z1(11)*epsl2+z1(12)*epsl3
      a3=z1(13)+z1(14)*epsl+z1(15)*epsl2+z1(16)*epsl3
      fac=a0+a1*yx0+a2*yx0*yx0+a3*yx0*yx0*yx0
      coeff3=fac*yx01**z1(17)
      return
  100 continue
      a0=z2(1)+z2(2)*epsl+z2(3)*epsl2+z2(4)*epsl3
      a1=z2(5)+z2(6)*epsl+z2(7)*epsl2+z2(8)*epsl3
      a2=z2(9)+z2(10)*epsl+z2(11)*epsl2+z2(12)*epsl3
      a3=z2(13)+z2(14)*epsl+z2(15)*epsl2+z2(16)*epsl3
      fac=a0+a1*yx0+a2*yx0*yx0+a3*yx0*yx0*yx0
      coeff3=fac*yx01**z2(17)
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
c     i=3,4 refer to W current

      if(i.gt.4) stop
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
      go to (1,2,3,4) i
    1 cheavy=(a+b*eps+c*eps2)*fac+(aa+bb*eps)*beta
      return
    2 cheavy=(-b*eps-2.*c*eps2)*fac+(a+b*eps+c*eps2)/beta
     .      +(-bb*eps)*beta +(aa+bb*eps)*2.*z*eps/z1/beta
      return
    3 f1=2.*z*z1*(1.-eps/zr)*eps1*(1.-0.5*eps) -1.
      f2=6.*z*z1*(1.-0.5*eps*(1.-eps)) 
     .  +6.*z*(1.-2.*z)*eps*dlog(zr/eps)
      arg=zr*zr/eps 
      f3=0.5*eps1*(1.-2.*z*eps1*(1.-eps1*z))*dlog(arg)
      cheavy=f1+f2+f3
      return
    4 eps3=eps2*eps
      arg=zr*zr/eps 
      f1=2.*z*(eps-eps2*(2.-3.*z)-1.5*eps3*z)
      f2=0.5*(a+eps*(1.-4.*z+6.*z2)-2.*eps2*z*(1.-3.*z)
     .  +2.*eps3*z2) 
     .  +6.*eps*z*(1.-2.*z)*(1.-dlog(zr/eps))
      f3=0.5*(-eps*(1.-4.*z+6.*z2)+4.*eps2*z*(1.-3.*z)
     .   -6.*eps3*z2)*dlog(arg)
      cheavy=f1+f2+f3
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


