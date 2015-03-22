c choix des differents parametres
	subroutine param
	implicit real*8 (a-h,l-z)
	logical lo,nlo,gener,integ,intega
	logical dir,onef,twof
	logical distrib,resom
	logical pt_pair,pt_photon
	logical facteur_symetrie
        common/alfa/masch2,masbo2,masto2,lambda4square
	common/parami/ysup,yinf,gptsup,gptinf,s
	common/gdf/jnf
	common/coul/n,cf,gtr
	common/hadron/ih1,ih2,ih3,ih4
	common/coup/ptm,r
	common/flag/icoup
	common/coupexp/ptm_exp,r_exp
	common/flagexp/icoup_exp
	common/w50512/qcdl4,qcdl5
	common/alem/iloopem
	common/faconv/hc2
	common/nboucle/iloop
	common/fscheme/ischeme
	common/aurenche/iauren
	common/box/ibox
	common/processd/j_processd_min,j_processd_max
	common/processo/j_processo_min,j_processo_max
	common/processt/j_processt_min,j_processt_max
	common/approx/lo,nlo
	common/calcul/integ,intega,gener
	common/nbevent/inbevent
	common/baseparam/jtmx1,jtmx2,jcall,jcall1,jcall2,ixtry
	common/accuracy/accu
	common/typesc/ichoi_scale
	common/normasc/cm,cmu,cmf
	common/pclass/dir,onef,twof
	common/phocut/iphocut
	common/phocut1/mmin,mmax,dmin
	common/ficut/fimin
	common/coup_histo/r_isol,etmax
        character*20 parm(20)
	common/pdf/value(20),valuep(20)
	character*128 path_bsfile,path_rzfile
	common/cheminbs/path_bsfile
	common/long/ilen
	common/cheminrz/path_rzfile
	common/longrz/ilenrz
	common/calchoi/distrib,resom
	common/resompar/q0,g1,g2,c1,c2,c,bmax,g3
	common/resompargg/q0gg,g1gg,g2gg,c1gg,c2gg,cgg,bmaxgg
	common/diffpt/pt_pair,pt_photon
	common/symetrie/facteur_symetrie
	character*128 name_experiment
	common/name/name_experiment
	common/lname/ilname
	common/flag_iso/i_flag_iso
c nom pour le chemin ou se trouverons les fichiers .bs
	read *,path_bsfile
	ilen = icharnum(path_bsfile)	
	open(unit=8,file=path_bsfile(1:ilen)//'output.param',
     #	status='unknown')
	write (8,*) 'ilen=',ilen
	write (8,*) 'path_bsfile=',path_bsfile(1:ilen)
c nom pour le chemin ou se trouverons les fichiers .rz et .res
	read *,path_rzfile
	ilenrz = icharnum(path_rzfile)	
	write (8,*) 'ilenrz=',ilenrz
	write (8,*) 'path_rzfile=',path_rzfile(1:ilenrz)
c determination dans pdflib des fonctions de distribution dans 
c le hadron h1
	parm(1) = 'nptype'
	read *,value(1)
	write(8,*) 'value(1)=',value(1)
	parm(2) = 'ngroup'
	read *,value(2)
	write(8,*) 'value(2)=',value(2)
	parm(3) = 'nset'
	read *,value(3)
	write(8,*) 'value(3)=',value(3)
	call pdfset(parm,value)
c determination dans pdflib des fonctions de distribution dans 
c le hadron h2
	read *,valuep(1)
	write(8,*) 'valuep(1)=',valuep(1)
	read *,valuep(2)
	write(8,*) 'valuep(2)=',valuep(2)
	read *,valuep(3)
	write(8,*) 'valuep(3)=',valuep(3)
c alpha electromagnetique, iloopem = 0 alpha_em = 1/137, autrement
c alpha_em est running comme dans jetset
	read *,iloopem
	write(8,*) 'iloopem=',iloopem
c ibox = 2 le terme de born de twophod.f est en fait: born + box
c ibox = 0 le terme de born de twophod.f est: born
c ibox = 1 le terme de born de twophod.f est: box
	read *,ibox
	write(8,*) 'ibox=',ibox
c saveurs actives
	read *,jnf
	write(8,*) 'jnf=',jnf
c masses des differents quarks. ces masses rentrent dans alpha_s pour
c savoir le nombre de saveurs actives.
	masch2 = (1.5d0)**2
	if (jnf.eq.5) then
	  masbo2 = (4.5d0)**2
	else if (jnf.eq.4) then
	  masbo2 = (1.d+10)**2
	endif
	masto2 = (1.d+10)**2
	lambda4square = qcdl4**2
	if ((value(1).eq.1.d0).and.(value(2).eq.4.d0).and.
     #	(value(3).eq.34.d0)) then
	  lambda4square = (.296)**2
	endif
c facteur (hbarr*c)**2
	read *,hc2
	write(8,*) 'hc2=',hc2
c nombre de boucles dans l'evaluation de alpha_s
	read *,iloop
	write(8,*) 'iloop=',iloop
c ischeme: schema de factorisation pour les pattes initiales
c 0: msbarr; 1:dis
	read *,ischeme
	write(8,*) 'ischeme=',ischeme
c pour retrouver les calculs de aurenche et al., il faut prendre 2 etats
c de polarisation pour le gluon: iauren=1, sinon convention standard
c 2*(1-epsilon) degree de polarisation: iauren=0
	read *,iauren
	write(8,*) 'iauren=',iauren
c ici : pour choisir entre collider et cible fixe
c ici = 0 collider, ici = 1 ciblefixe
	read *,ici
	write (8,*) 'ici=',ici
c energie disponible dans le cdm proton-(anti)proton pour collider
c ou ebeam pour experience cible fixe
	read *,ebeam
c rapidite max. des photons
	read *,ysup
	write(8,*) 'ysup=',ysup
c rapidite min. des photons
	read *,yinf
	write(8,*) 'yinf=',yinf
c impulsion transverse max. des photons
	read *,gptsup
	write(8,*) 'gptsup=',gptsup
c impulsion transverse min. des photons
	read *,gptinf
	write(8,*) 'gptinf=',gptinf
c constantes de su(3)_couleur
	n = 3.d0
	cf = (n*n-1.d0)/(2.d0*n)
	gtr = dfloat(jnf)/2.d0
c ih1 et ih2 determine le choix des particules initiales, ih = 0:
c proton; ih = 1: antiproton, ih = 2: photon, ih = 3: pion, ih = 4: Be
c ih2 est la particule cible, ih1 est la particule incidente
	read *,ih1
	write(8,*) 'ih1=',ih1
	read *,ih2
	write(8,*) 'ih2=',ih2
c calcul de racine de s
	if (ici.eq.0) then
	  rs = ebeam
	  s = rs**2
	  write(8,*) 'rs =',rs
	else if (ici.eq.1) then
	  xmp = .93828d0
	  xmpi = .139567d0
	  if (ih1.le.1) then
	    xm_inci = xmp
	  else if (ih1.eq.3) then
	    xm_inci = xmpi
	  else
	    write (*,*) 'bad choice of incident particule'
	  endif
	  if (ih2.le.1) then
	    xm_cible = xmp
	  else if (ih2.eq.4) then
	    xm_cible = xmp
	  else
	    write (*,*) 'bad choice of target particule'
	  endif
	  rs = dsqrt(xm_inci**2+xm_cible**2+2.d0*xm_cible*ebeam)
	  s = rs**2
	  write(8,*) 'ebeam =',ebeam
	  write(8,*) 'rs =',rs
	endif
c ih3 et ih4 determine le choix des fonctions de fragmentation de
c parton en hadron:
c
c ih3,ih4 = 
c
c 150  : bkk, pi_0,  lo
c 151  : bkk, pi_0, nlo
c 
c 250  : kkp, pi_0,  lo
c 251  : kkp, pi_0, nlo
c
c 301  : owens,   gamma, nlo
c
c 401  : bouhris, gamma, nlo, set1
c 402  : bouhris, gamma, nlo, set2
c
c
	read *,ih3
	write(8,*) 'ih3=',ih3
	read *,ih4
	write(8,*) 'ih4=',ih4
c si ih3 = ih4, alors on met les facteur de symetrie
c c'est a dire on met a true la variable facteur_symetrie
	if (ih3.eq.ih4) then
	  facteur_symetrie = .true.
	else
	  facteur_symetrie = .false.
	endif
	write(8,*) 'facteur de symetrie=',facteur_symetrie
c coupures experimentale d'isolation: 0 pas de coupures; 1 coupures a
c la cdf (ca ne sert a rien dans cette version du programme)
	ptm_exp = 2.d0
	r_exp = .7d0
	icoup_exp = 0
c coupures theoriques en pt et r
	icoup = 0
	read *,ptm
	write(8,*) 'ptm=',ptm
	read *,r
	write(8,*) 'r=',r
c selection des sousprocessus
c partie directe
c	j0 = 1 : qi + qbi --> ph + ph
c	j0 = 2 : qi + g --> ph + ph
	read *,j_processd_min
	write(8,*) 'j_processd_min=',j_processd_min
	read *,j_processd_max
	write(8,*) 'j_processd_max=',j_processd_max
c partie un brem
c	qi + qk --> qi + ph
c	  j0 = 1 : d + u --> d + ph
c	  j0 = 2 : d + dp --> d + ph
c	  j0 = 3 : u + d --> u + ph
c	  j0 = 4 : u + up --> u + ph
c	qi + qbk --> qi + ph
c	  j0 = 5 : d + ub --> d + ph
c	  j0 = 6 : d + dpb --> d + ph
c	  j0 = 7 : u + db --> u + ph
c	  j0 = 8 : u + upb --> u + ph
c	qi + qbi --> qk + ph
c	  j0 = 9 : d + db --> u + ph
c	  j0 = 10 : d + db --> dp + ph
c	  j0 = 11 : u + ub --> d + ph
c	  j0 = 12 : u + ub --> up + ph
c	j0 = 13 : qi + qi --> qi + ph
c	j0 = 14 : qi + qbi --> qi + ph
c	j0 = 15 : qi + g --> qi + ph
c	j0 = 16 : qi + g --> g + ph
c	j0 = 17 : qi + qbi --> g + ph
c	j0 = 18 : g + g --> qi + ph
	read *,j_processo_min
	write(8,*) 'j_processo_min=',j_processo_min
	read *,j_processo_max
	write(8,*) 'j_processo_max=',j_processo_max
c partie deux brem	
c	j0 = 1 : qi + qk --> qi + qk
c	j0 = 2 : qi + qk --> qi + g
c	j0 = 3 : qi + qbk --> qi + qbk
c	j0 = 4 : qi + qbk --> qi + g
c	j0 = 5 : qi + qbi --> qk + qbk
c	j0 = 6 : qi + qbi --> qk + g
c	j0 = 7 : qi + g --> qi + qk
c	j0 = 8 : qi + g --> qi + qbk
c	j0 = 9 : qi + g --> qk + qbk
c	j0 = 10 : qi + qi --> qi + qi
c	j0 = 11 : qi + qi --> qi + g
c	j0 = 12 : qi + qbi --> qi + qbi
c	j0 = 13 : qi + qbi --> qi + g
c	j0 = 14 : qi + g --> qi + qi
c	j0 = 15 : qi + g --> qi + qbi
c	j0 = 16 : qi + qbi --> g + g
c	j0 = 17 : qi + g --> qi + g
c	j0 = 18 : qi + g --> g + g
c	j0 = 19 : g + g --> qi +qbi
c	j0 = 20 : g + g --> qi + g
c	j0 = 21 : g + g --> g + g
	read *,j_processt_min
	write(8,*) 'j_processt_min=',j_processt_min
	read *,j_processt_max
	write(8,*) 'j_processt_max=',j_processt_max
c
	read *,nlo
	write(8,*) 'nlo=',nlo
c
	read *,lo
	write(8,*) 'lo=',lo
c si integ est true il n'y a que l'integration (section efficace
c integree)
	read *,integ
	write(8,*) 'integ=',integ
c si intega est true il n'y a que l'integration (valeur absolue de la
c section efficace)
	read *,intega
	write(8,*) 'intega=',intega
c si gener est true alors on a la generation et l'integration	
	read *,gener
	write(8,*) 'gener=',gener
c inbevent est le nb d'evenements generes
 	read *,inbevent
	write(8,*) 'inbevent=',inbevent
c jtmx1 et jtmx2 sont le nb d'iteration dans bases pour la grille et
c l'integrale
  	read *,jtmx1
	write(8,*) 'jtmx1=',jtmx1
	read *,jtmx2
	write(8,*) 'jtmx2=',jtmx2
c jcall,jcall1 et jcall2 sont le nb d'appelles par iteration pour
c respectivement les quasi deux donne deux, les deux donne trois et les
c vrai deux donne deux
  	read *,jcall1
	write(8,*) 'jcall1=',jcall1
  	read *,jcall
	write(8,*) 'jcall=',jcall
  	read *,jcall2
	write(8,*) 'jcall2=',jcall2
c ixtry est le nb d'essai pour la rejection d'evenement dans spring
  	read *,ixtry
	write(8,*) 'ixtry=',ixtry
c accu est la precion en pour cent de bases
	read *,accu
	write(8,*) 'accu=',accu
c choix de l'echelle: ichoi_scale selectionne les familles d'echelle: 1
c (pt3+pt4)*n, 2 sqrt(pt3^2+pt4^2), 3 mgg (masse invariante gamma-gamma)
  	read *,ichoi_scale
	write(8,*) 'ichoi_scale=',ichoi_scale
c cm,cmu et cmf sont des facteurs de normalisation des differents choix
c pour respectivement les echelles de factorisation initial, de
c renormalisation et de factorisation final
	read *,cm
	write(8,*) 'cm=',cm
	read *,cmu
	write(8,*) 'cmu=',cmu
	read *,cmf
	write(8,*) 'cmf=',cmf
c dir, onef et twof servent a selectionner les parties directes, one
c brem et two brem
	read *,dir
	write (8,*) 'dir=',dir
	read *,onef
	write (8,*) 'onef=',onef
	read *,twof
	write (8,*) 'twof=',twof
c coupure entre les 2 photons
c iphocut=0 angle azymuthal > fimin, iphocut = 1 masse invariante
c photon-photon > mmin et < mmax, iphocut = 2 distance en pseudo-rapidite angle
c azymuthal > dmin
	read *,iphocut
	write(8,*) 'iphocut=',iphocut
	read *,fimin
	write(8,*) 'fimin=',fimin
	read *,mmin
	write(8,*) 'mmin=',mmin
	read *,mmax
	write(8,*) 'mmax=',mmax
	read *,dmin
	write(8,*) 'dmin=',dmin
c coupure d'isolement au niveau de bases
c r_isol est la taille du cone d'isolement
c etmax est le maximum d'energie transverse depose dans le cone
	read *,i_flag_iso
	write(8,*) 'i_flag_iso=',i_flag_iso
	read *,r_isol
	write(8,*) 'r_isol=',r_isol
	read *,etmax
	write(8,*) 'etmax=',etmax
c
c si distrib est vrai, on calcul les termes en distribution, si resom
c est vrai alors on calcul la resomation de ces distribution
	read *,distrib
** 	write(8,*) 'distrib=',distrib
	read *,resom
** 	write(8,*) 'resom=',resom
c valeur de q0
	read *,q0
** 	write(8,*) 'q0=',q0
	q0 = 2.d0
c valeur de bmax, en general 1/q0
	bmax = 1.d0/q0
c valeur de g1
	read *,g1
** 	write(8,*) 'g1=',g1
c valeur de g2
	read *,g2
** 	write(8,*) 'g2=',g2
	read *,g3
** 	write(8,*) 'g3=',g3
c	g3 = 0.d0
	gammae = .5772156649
	c = 2.d0*dexp(-gammae)
c valeur de c1
	read *,c1
** 	write(8,*) 'c1=',c1
	if (c1.eq.0.d0) then
	  c1 = c
	else
	  c1 = c1
	endif
c valeur de c2
	read *,c2
** 	write(8,*) 'c2=',c2
c pt_pair vrai on calcule dsigma/d qt, si pt_photon vrai, on calcul dsigma/dpt
	read *,pt_pair
** 	write(8,*) 'pt_pair=',pt_pair
	read *,pt_photon
** 	write(8,*) 'pt_photon=',pt_photon
c valeur des parametres q0,g1,g2,c1,c2,bmax pour gg --> photon photon
c on suppose dans un premier temps qu'ils sont egaux aux parametres de
c qqb --> photon photon
c valeur de g1 pour gg
	read *,g1gg
** 	write(8,*) 'g1gg=',g1gg
c valeur de g2 pour gg
	read *,g2gg
** 	write(8,*) 'g2gg=',g2gg
	q0gg = q0
	g1gg = g1
	g2gg = g2
	c1gg = c1
	c2gg = c2
	cgg = c
	bmaxgg = bmax
c 
	read *,name_experiment
	ilname = icharnum(name_experiment)	
** 	write (8,*) 'ilname =',ilname
** 	write (8,*) 
**      #	'name_experiment=',name_experiment(1:ilname)
c
	close (8)
	return
	end
	
