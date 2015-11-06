      subroutine bonn
c
c        Note: this version of bonn.f includes 'save'
c              statements in all subroutines where needed;
c              03/22/07.
c
c        bonn and its subroutines (ob...) compute the
c        relativistic one-boson-exchange nn-nn potential 
c        in momentum space; only bonn needs to be called.
c        this version of the code uses numerical integration
c        for the partial wave decomposition.
c        this package contains all subroutines needed.
c
c
c        this version of the code is published in "computational 
c        nuclear physics 2---nuclear reactions", Langanke, Maruhn, 
c        Koonin, eds. (springer, new york, 1993)
c
c
c        author:      r. machleidt
c                     department of physics
c                     university of idaho
c                     moscow, idaho 83844
c                     u. s. a.
c                     e-mail: machleid@uidaho.edu
c
c                     formerly:
c                     institut fuer theoretische kernphysik der
c                     universitaet bonn
c                     nussallee 14-16
c                     d - 5300  bonn, w. germany
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine:
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
c
c
c        xmev and ymev are the final and initial relative momenta,
c        respectively, in units of mev.
c        v is the potential in units of mev**(-2).
c        if heform=.true., v contains the 6 helicity matrix elements
c        associated with one j in the following order:
c        0v, 1v, 12v, 34v, 55v, 66v (helicity formalism).
c        if heform=.false., v contains the 6 lsj-state matrix elements
c        associated with one j in the following order:
c        0v, 1v, v++, v--, v+-, v-+ (lsj-formalism).
c        j is the total angular momentum. there is essentially no upper
c        limit for j.
c        sing, trip, and coup should in general be .true..
c        endep and label can be ignored.
c
c
c        this has been the end of the common-blocks containing
c        the arguments and values of this subroutine
c
c        specifications for these two common blocks
c
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
c
      data pi/3.141592653589793d0/
      character*4 mesong(12)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            '1+  ','2+  '/
      logical index
      data index/.false./
      logical indsig
      logical indmg(12)
      data indmg/12*.false./
c
      save
c
c
c
c
      inter=1
c
c
c
c
c        call subroutine obpar once and only once
c
c
      if (index) go to 50
      index=.true.
c
c
      call obpar
c
c
      iftgo=ift(inter)+1
      dwn=1.d0/wnn(inter)
      iman=imaa(inter)
      imen=imea(inter)
      imenn=imen
      indsig=indpar(inter)
      if (indsig) imenn=imen-2
c
c
c        prepare constant over-all factor
c
      fac=1.d0/(2.d0*pi)*dwn*dwn
c     --------------------------
c
c
c
c
c
c
c
c        prepare expressions depending on x and y
c        ----------------------------------------
c        ----------------------------------------
c
c
c
c
   50 x=xmev*dwn
      y=ymev*dwn
      indxy=.false.
      xx=x*x
      yy=y*y
      xy2=x*y*2.d0
      xxpyy=xx+yy
      ex=dsqrt(1.d0+xx)
      ey=dsqrt(1.d0+yy)
      eem12=(ex*ey-1.d0)*2.d0
c
c
c
c
      xy=xy2*0.5d0
      ee=ex*ey
      ree=dsqrt(ee)
      eem1=ee-1.d0
      eme=ex-ey
      emeh=eme*0.5d0
      emehq=emeh*emeh
      eep1=ee+1.d0
       epe=ex+ey
      xxyy=xx*yy
c
c
c
c
c        prepare over-all factor
c
c
      go to (70,71,72),iftgo
c
c        no additional factor
c
   70 fff=fac
      go to 90
c
c        minimal relativity
c
   71 fff=fac/ree
      go to 90
c
c        factor m/e*m/e
c
   72 fff=fac/ee
c
c
c
c
c
c
   90 do 93 iv=1,6
   93 v(iv)=0.d0
      do 95 il=iman,imen
      do 95 iv=1,32
   95 vj(iv,il)=0.d0
c
c
c
c
c        contributions of mesons
c        -----------------------
c        -----------------------
c
c
c
c
      do 1995 img=1,mge
      mg=mggo(img,inter)
      if (mg.eq.0) go to 2000
      if (mg.gt.7) go to 9000
      me=mgg(mg,inter)
      go to (100,200,9000,400,9000,9000,700),mg
c
c
c
c        0-  , pseudo-scalar coupling
c        ----------------------------
c
c
c
c
  100 mc=1
c
      ff=1.d0
      f(1)=eem1
      f(2)=-xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c        0-t , pseudo-vector coupling
c        ----------------------------
c
c
c
c
  200 mc=1
c
      ff=1.d0
      f(1)=eem1+emehq*(ee+3.d0)
      f(2)=-xy+emehq*xy
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=-eme-eme*(emehq+eem1)
      f(8)=-f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        0+  , scalar coupling
c        ---------------------
c
c
c
c
  400 mc=1
c
      ff=1.d0
      f(1)=-eep1
      f(2)=xy
      f(3)=f(1)
      f(4)=f(2)
      f(5)=f(2)
      f(6)=f(1)
      f(7)=epe
      f(8)=f(7)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        1-t , vector mesons
c        -------------------
c
c
c
c
c        vector-vector coupling
c
c
c
c
  700 mc=1
c
      ff=2.d0
      f(1)=eem1+ee
      f(2)=0.d0
      f(3)=ee
      f(4)=xy
      f(5)=xy2
      f(6)=1.d0
      f(7)=-ey
      f(8)=-ex
c
      call obstr(1,1,me)
c
c
c
c
c        tensor-tensor coupling
c
c
c
c
      mc=2
c
      ff=0.25d0
      f(1)=(3.d0*ee+1.d0)*xxpyy
      f(2)=-(6.d0*ee+2.d0-xxpyy)*xy
      f(3)=eem1*xxpyy+4.d0*xxyy
      f(4)=-(4.d0*ee+xxpyy)*xy
      f(5)=(4.d0-3.d0*xxpyy)*xy
      f(6)=6.d0*xxyy-(ee+3.d0)*xxpyy
      f(7)=(ex+3.d0*ey)*xx+eme*yy
      f(8)=(ey+3.d0*ex)*yy-eme*xx
c        factors for additional terms
      f(9)=-2.d0*xxyy
      f(10)=eep1*xy2
      f(11)=-epe*xy2
c
      call obstr(2,1,me)
c
c
c
c
c        vector-tensor coupling
c
c
c
c
      mc=3
c
      ff=1.d0
      f(1)=xxpyy
      f(2)=-xy2
      f(3)=-f(1)
      f(4)=-f(2)
      f(5)=6.d0*xy
      f(6)=3.d0*f(3)
      f(7)=(ex*yy+3.d0*ey*xx)
      f(8)=(ey*xx+3.d0*ex*yy)
c
      call obstr(1,1,me)
      go to 1995
c
c
c
c
c        this has been the end of the contributions of mesons
c        ----------------------------------------------------
c
c
c
c
c        errors and warnings
c        -------------------
c
c
c
c
 9000 if (indmg(mg)) go to 1995
      write (kwrite,19000) mesong(mg)
19000 format(////' warning in bonn: meson-group  ',a4,'  does not exi
     1st in this program.'/' contribution ignored. execution continued.'
     2////)
      indmg(mg)=.true.
c
c
c
c
 1995 continue
c
c
c
c
c        add up contributions of mesons
c        ------------------------------
c
c
c
c
 2000 do 2005 il=iman,imenn
      do 2005 iv=1,6
 2005 v(iv)=v(iv)+vj(iv,il)
c
      if (.not.indsig) go to 2100
c
c        case of different sigmas for t=1 and t=0
c
      mmod=mod(j,2)
      imxx=imenn+1+mmod
      imyy=imenn+2-mmod
      do 2025 iv=1,6
      if (iv.eq.2) go to 2023
      v(iv)=v(iv)+vj(iv,imxx)
      go to 2025
 2023 v(iv)=v(iv)+vj(iv,imyy)
 2025 continue
c
 2100 continue
c
c
c
c
c****    for cutting out the attractive short-range part in 1p1
c****    enable the statement below
c**** if(j.eq.1.and.v(1).lt.0.d0) v(1)=0.d0
c        this can be useful for more numerical stability in 1p1.
c
c
      return
      end
      subroutine obpar
c
c        obpar reads, writes and stores the parameters for the 
c        one-boson-exchange potential.
c
c
      implicit real*8 (a-h,o-z)
c
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
c
      dimension cc(5)
      integer name(3),nname(15)
      integer imga(3)
      integer cut,cutg,end,stars
      data cut/'cut '/,cutg/'cutg'/,end/'end '/,stars/'****'/
      integer mesong(12)
      data mesong/'0-  ','0-t ','0-st','0+  ','0+st',
     1            '1-  ','1-t ','1-tt','1-st','1-ss',
     2            '1+  ','2+  '/
      logical index
      data index/.false./
      logical zerocp,indcut
      data zerocp/.true./,indcut/.false./
      data uf/197.3286d0/
c
      save
c
c
c
c
10000 format (2a4,a2,15a4)
10002 format (/' jp  name      g**2      f/g       mass    isospin
     1   iprop'/9x,'cut type      c u t - o f f   p a r a m e t e r s')
10003 format (2a4,a2,5f10.4)
10004 format (1h ,2a4,a2,2f10.4,f9.2,1x,2(f7.1,3x))
10005 format (1h ,2a4,a2,f3.1,f11.1,f9.4,f14.4,f10.4)
10006 format (2a4,a2,3i3)
10007 format (1h ,2a4,a2,3i3)
10008 format (1h ,57(1h-))
10011 format (1h1  //' bonn:  one-boson-exchange nn-nn potential (numer.
     1 integ.)')
10015 format (' input-parameter-set:'/1h ,20(1h-))
10016 format (1h ,2a4,a2,15a4)
c
c
c
c
      if (index) go to 50
      index=.true.
c
      x=-1.d0
      y=-1.d0
c
c
c
c
c        maxima of certain indices related to the dimension as follows:
c        dimension c(mme,imee),ic(mice,imee),indc(mindce,imee),
c                  mgg(mge,3),mggo(mge,3),mesong(mge),vj(32,imee),
c                  ima(mee,mge,3)
c
      mge=12
      mee=5
      mme=10
      mice=10
      mindce=2
      imb=1
      ime=0
      imee=15
      imec=0
c        mme always ge mice, mindce
c
c        set all meson-parameters and indices to zero or .false.
c
      do 1 int=1,3
      imga(int)=0
      indpar(int)=.false.
      do 1 mgx=1,mge
      mgg(mgx,int)=0
    1 mggo(mgx,int)=0
c
c
      do 2 il=1,imee
      do 2 mm=1,mme
      if (mm.le.mindce) indc(mm,il)=.false.
      if (mm.le.mice) ic(mm,il)=0
    2 c(mm,il)=0.d0
      endep=.false.
c
c
c
c
c
c
c        reading and writing of first 4 records
c        --------------------------------------
c        --------------------------------------
c
c
c
c        write headline and read and write name of parameter set
c
   50 write (kwrite,10011)
      write (kwrite,10008)
      write (kwrite,10015)
      read  (kread, 10000) name,nname
      write (kwrite,10016) name,nname
      label=name(1)
c
c        read and write index-parameter concerning the factor of the
c        potential
c
      read  (kread, 10006) name,ift(inter)
      write (kwrite,10007) name,ift(inter)
      iftyp=ift(inter)
      if (iftyp.lt.0.or.iftyp.gt.2) go to 9003
c
c        read and write parameters for numerical integration
c
      read  (kread, 10006) name,mint(inter),maxt(inter)
      write (kwrite,10007) name,mint(inter),maxt(inter)
c
c        read and write mass of nucleon
c
      read  (kread, 10003) name,wn
      write (kwrite,10004) name,wn
      wnq=wn*wn
      dwn=1.d0/wn
      dwnq=dwn*dwn
      wnn(inter)=wn
c
c        write headline for meson parameters
c
      write (kwrite,10002)
      write (kwrite,10008)
c
c
c
c
c        read, write and store meson parameters
c        --------------------------------------
c        --------------------------------------
c
c
c
   61 read  (kread, 10003) name,cc
c
c        check if record just read contains cut-off parameters
c
      if (name(1).eq.cut.or.name(1).eq.cutg) go to 70
c
c        check if ****, indicating that there is a different sigma
c        for t=1 and t=0
c
      if (name(1).ne.stars) go to 69
      indpar(inter)=.true.
      write (kwrite,10004) name
      go to 61
c
c        check if end of mesons
c
   69 if (name(1).eq.end) go to 2000
c
c
c
c
c        write meson-parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
      indcut=.false.
c
      write (kwrite,10004) name,cc
c
c        check if coupling constants are zero
c
      if (cc(1).ne.0.d0) go to 62
      zerocp=.true.
      go to 61
c
   62 zerocp=.false.
c
c        find out number of meson-group mg
c
      do 63 mg=1,mge
      if (name(1).eq.mesong(mg)) go to 64
   63 continue
      go to 9000
c
c
c
c
c        store meson parameters, which are no cut-off parameters
c        -------------------------------------------------------
c
c
c
c
   64 ime=ime+1
      if (ime.gt.imee) go to 9011
      mgg(mg,inter)=mgg(mg,inter)+1
      m=mgg(mg,inter)
      if (m.gt.mee) go to 9001
      ima(m,mg,inter)=ime
      if (m.ne.1) go to 65
      imga(inter)=imga(inter)+1
      mggo(imga(inter),inter)=mg
   65 continue
c
c        store coupling constant g**2/4pi
      c(1,ime)=cc(1)
c        store coupling constant f*g/4pi
      c(3,ime)=cc(1)*cc(2)
c        store coupling constant f**2/4pi
      c(2,ime)=cc(2)*c(3,ime)
c        store meson mass squared in units of nucleon mass squared
      c(4,ime)=cc(3)*cc(3)*dwnq
c
c        get iso-spin
      icc=cc(4)
      if (icc.ne.0.and.icc.ne.1) go to 9004
c         store isospin as logical constant
      if (icc.eq.1) indc(1,ime)=.true.
c        store parameter for meson propagator (iprop)
      ic(1,ime)=cc(5)       
      if (ic(1,ime).ne.0) go to 9005
c
c        index values for further storing
      mi=4
      mm=5
      go to 61
c
c
c
c
c        write cut-off parameters
c        ------------------------
c
c
c
c
   70 write (kwrite,10005) name,cc
c
c
c        check if individual cut or general cut
c
      if (name(1).eq.cut) go to 73
c        case of general cut-off
      if (indcut) go to 90
      if (imec.ge.ime) go to 61
      imac=imec+1
      imec=ime
      if (imac.lt.imb) imac=imb
      go to 90
c        case of individuel cut-off
   73 imac=ime
      imec=ime
      if (zerocp) go to 61
c
   90 indcut=.true.
c
c        if cutoff type = 0, ignore cutoff
      if (cc(1).eq.0.d0) go to 61
c
c        save present values of indices
      mix=mi
      mmx=mm
c
c        start loop of mesons, which present cut-off refers to
c
      do 1095 im=imac,imec
      mi=mix
      mm=mmx
c
c
c
c
c        store cut-off parameters
c        ------------------------
c
c
c
c
c        store type of cut-off
      ic(mi,im)=cc(1)
      ityp=ic(mi,im)
      if (ityp.lt.1.or.ityp.gt.2) go to 9002
c        store and test type of denominator of cut-off
      ic(mi+1,im)=cc(2)
      if (ic(mi+1,im).ne.0) go to 9006
c
c
c        cut-off of monopole/dipole type
c        *******************************
c
c
c        store and test exponent of cut-off
      ic(mi+2,im)=cc(3)
      if (ic(mi+2,im).lt.0) go to 9009
      if (ic(mi+2,im).gt.0) go to 101
c        exponent is zero, omit cut-off
      ic(mi,im)=0
      ic(mi+1,im)=0
      go to 1000
c        store cut-off mass for denominator
  101 c(mm+1,im)=cc(4)*cc(4)*dwnq
c        store numerator of cut-off
      c(mm,im)=c(mm+1,im)
      if (ityp.eq.2)     c(mm,im)=c(mm,im)-c(4,im)
      mi=mi+3
      mm=mm+2
c
c
c
c
c        end cut-offs
c        ************
c
c        test dimensions
 1000 if (mi.gt.mice.or.mm-1.gt.mme) go to 9010
c
c
 1095 continue
      go to 61
c
c
c
c
c        last record
c        -----------
c        -----------
c
c
c
c
c        write end mesons
 2000 imaa(inter)=imb
      imea(inter)=ime
      imb=ime+1
      write (kwrite,10004) name
      write (kwrite,10008)
      write (kwrite,10008)
c
c
c
c
      return
c
c
c
c        errors
c        ------
c        ------
c
c
c
c
 9000 write (kwrite,19000) name(1)
19000 format (/////' error in obpar:  meson-group   ',a4,'   does not
     1 exist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9001 write (kwrite,19001)
19001 format (/////' error in obpar: too many mesons within a meson-g
     1roup with respect to'/' the given dimensions. execution termina
     2ted.'////)
      go to 9999
c
c
 9002 write (kwrite,19002) cc(1)
19002 format (/////' error in obpar: cut-off type',f10.4,' does not e
     1xist in this program.'/' execution terminated.'////)
      go to 9999
c
c
 9003 write (kwrite,19003) iftyp
19003 format (/////' error in obpar: factor type has the non-permissi
     1ble value',i4,' .'/' execution terminated.'////)
      go to 9999
c
c
 9004 write (kwrite,19004) cc(4)
19004 format (/////' error in obpar: isospin has the non-permissible
     1value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9005 write (kwrite,19005) cc(5)
19005 format (/////' error in obpar:     iprop has the non-permissibl
     1e value',f10.4,'  .'/' execution terminated.'////)
      go to 9999
c
c
 9006 write (kwrite,19006) cc(2)
19006 format (/////' error in obpar: the index for the denominator of
     1 the cut-off has the'/' non-permissible value',f10.4,' . execut
     2ion terminated.'////)
      go to 9999
c
c
 9009 write (kwrite,19009)
19009 format (/////' error in obpar: the exponent of the cut-off is l
     1ess than zero.'/' execution terminated.'////)
      go to 9999
c
c
 9010 write (kwrite,19010)
19010 format (/////' error in obpar: too many cut-off parameters with
     1 respect to the given'/' dimensions. execution terminated.'////)
      go to 9999
c
c
 9011 write (kwrite,19011)
19011 format (/////' error in obpar:  too many mesons with respect to
     1 the dimensions given'/' in this program. execution terminated.'
     2////)
      go to 9999
c
c
 9999 stop
      end
      subroutine obstr (icase,max,mex)
c
c        obstr computes the structure of one-boson-exchanges
c
c
      implicit real*8 (a-h,o-z)
c
c
c        common blocks
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c     further specifications
c
      dimension vv(32)
      dimension tt(2,3)
      data jj/-1/
      logical index
      data index/.false./
      logical indiso
c
      save
c
c
c
c
      if (index) go to 50
      index=.true.
c
c
      tt(1,1)=1.d0
      tt(2,1)=-3.d0
c
      do 1 ii=2,3
      do 1 i=1,2
    1 tt(i,ii)=1.d0
c
c
c
c
c
   50 do 1095 m=max,mex
      im=ima(m,mg,inter)
c
c
      if (mc.ne.1) go to 60
c
c
c
c
c        call integrals
c        --------------
c
c
c
c
      call obai
c
c
c
c
   60 continue         
c
      if (c(mc,im).eq.0.d0) go to 1095
c
c
c
c
c        nn-nn helicity amplitudes
c        -------------------------
c
c
c        vv(1), ..., vv(6) contain in the following order:
c        0v, 1v, 12v, 34v, 55v, 66v.
c
c
c        basic structure
c
c
  100 ive=6
c
      vv(1)=f(1)*ai(1,m)+f(2)*ai(2,m)
      vv(2)=f(3)*ai(1,m)+f(4)*ai(3,m)
      vv(3)=f(5)*ai(1,m)+f(6)*ai(2,m)
      vv(4)=f(4)*ai(1,m)+f(3)*ai(3,m)
      vv(5)=f(7)*ai(4,m)
      vv(6)=f(8)*ai(4,m)
c
c
      go to (1000,120),icase
c
c
c        additional terms for the case of tensor-tensor coupling
c
c
  120 vv(1)=vv(1)+f(9)*ai(5,m)
      vv(2)=vv(2)+f(10)*ai(2,m)+f(9)*ai(6,m)
      vv(3)=vv(3)+f(10)*ai(5,m)
      vv(4)=vv(4)+f(9)*ai(2,m)+f(10)*ai(6,m)
         e1=f(11)*ai(7,m)
      vv(5)=vv(5)+e1
      vv(6)=vv(6)+e1
      go to 1000
c
c
c
c
c        set certain cases to zero 
c
 1000 if (j.ne.0) go to 1021
      vv(2)=0.d0
      vv(4)=0.d0
      vv(5)=0.d0
      vv(6)=0.d0
c
 1021 if (.not.sing) vv(1)=0.d0
      if (.not.trip) vv(2)=0.d0
      if (coup) go to 1030
      do 1025 iv=3,6
 1025 vv(iv)=0.d0
c
 1030 if (heform) go to 1040
c
c
c        transformation into lsj-formalism
c        (if requested)
      if (j.eq.jj) go to 1035
      jj=j
      aj=dfloat(j)
      aj1=dfloat(j+1)
      d2j1=1.d0/dfloat(2*j+1)
      arjj1=dsqrt(aj*aj1)
c
 1035 v3=vv(3)
      v4=vv(4)
      v5=vv(5)
      v6=vv(6)
      v34=arjj1*(v3-v4)
      v56=arjj1*(v5+v6)
      vv(3)=d2j1*(aj1*v3+aj*v4-v56)
      vv(4)=d2j1*(aj*v3+aj1*v4+v56)
      vv(5)=d2j1*(v34+aj1*v5-aj*v6)
      vv(6)=d2j1*(v34-aj*v5+aj1*v6)
c        
c        after transformation into lsj formalism,
c        vv(3), ..., vv(6) contain:
c        v++, v--, v+-, v-+.
c
c
c
c
c        multiply with factors
c        ---------------------
c
c
c
c
 1040 is=mod(j,2)+1
      it=mod(is,2)+1
      indiso=indc(1,im)
c        get coupling constant
      cmc=c(mc,im)
      fc=fff*ff*cmc
      do 1045 iv=1,ive
c
c        multiply with coupling-constant and factors fff and ff
c
      vv(iv)=vv(iv)*fc
c
c        multiply with isospin factor
c
      if (.not.indiso) go to 1045
      if (iv.eq.2) go to 1043
      vv(iv)=vv(iv)*tt(is,inter)
      go to 1045
 1043 vv(iv)=vv(iv)*tt(it,inter)
c
c
c        add up in case of several couplings for one meson and store
 1045 vj(iv,im)=vj(iv,im)+vv(iv)
c
c
 1095 continue
c
c
      return
      end
      subroutine obai
c
c        obai performs the numerical integration over angle theta
c        necessary for the partial wave decomposition.
c
c
      implicit real*8 (a-h,o-z)
c
      common /cstate/ j,heform,sing,trip,coup,endep,label
      logical heform,sing,trip,coup,endep
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c        further specifications
      dimension gi(5)
c
      dimension pj(5,96)
      real*4 axy2,aomq,am
      data nnt/-1/,jj/-1/
      logical indj
      logical index
      data index/.false./
c
      save
c
c
c
c
      if (index) go to 50
      index=.true.
      min=mint(inter)
      max=maxt(inter)
c
      igeint=5
c
      wn=wnn(inter)
      dwn=1.d0/wn
      wnq=wn*wn
c
c
c
c
   50 if (j.eq.jj) go to 70
      jj=j
      indj=.false.
c
c
      aj=dfloat(j)
      aj1=dfloat(j+1)
      dj1=1.d0/aj1
      ajdj1=aj*dj1
      aaj=dsqrt(ajdj1)
c
c
c
c
c        find out appropriate number of gauss-points, nt
c        -----------------------------------------------
c
c
c
c        c4 is the meson mass squared
c
   70 c4=c(4,im)
c
c
c        compute am
c
c
      axy2=xy2
      aomq=xxpyy+c4
      am=axy2/aomq
c
c
c
c        compute number of gausspoints (nt)
c        necessary for sufficient numerical accuracy of the integrals
c
c
      if (am.gt.0.999) go to 94
c
c
      if (am.gt.0.85) am=am**(-log(1.-am)-0.9)
c
c
      nt=float(min)/(1.-am)+0.9
c
c
      if (nt.gt.max) nt=max
      go to 95
c
c
   94 nt=max
c
c
   95 nt=nt+j
c
c        compute nt, which is suitable for gset
c
      if (nt.le.16) go to 98
      if (nt.gt.24) go to 96
      nt=4*(nt/4)
      go to 98
   96 if (nt.gt.48) go to 97
      nt=8*(nt/8)
      go to 98
   97 nt=16*(nt/16)
      if (nt.gt.96) nt=96
c
   98 if (nt.eq.nnt.and.indj) go to 100
c
c
c
c
c        call gauss-points
c        -----------------
c
c
c
c
      call gset (-1.d0,1.d0,nt,ct,wt)
      nnt=nt
c
c
c
c
c        call legendre-polynominals if necessary
c        -----------------------------------
c
c
c
c
      indxy=.false.
      indj=.true.
      do 99 i=1,nt
      t=ct(i)
      call legp (pj(1,i),pj(3,i),t,j)
      pj(2,i)=pj(1,i)*t
      pj(4,i)=pj(2,i)*t
   99 pj(5,i)=pj(3,i)*t
c
c
c
c
c        call integrand
c        --------------
c
c
c
c
  100 call obaa
c
c
c
c
c        prepare for integration
c
c
c
c
      do 2001 ig=1,igeint
 2001 gi(ig)=0.d0
c
c
c
c
c        integration-loop of cosine theta
c        --------------------------------
c
c
c
c
      do 2005 i=1,nt
      do 2005 ig=1,igeint
 2005 gi(ig)=gi(ig)+pj(ig,i)*aa(i)
c
c
c
      if (j.ne.0) go to 2010
      gi(3)=0.d0
      gi(5)=0.d0
c
c
c
c
c        combinations of integrals
c        -------------------------
c
c
c
c
 2010 ai(1,m)=gi(1)
c
      ai(2,m)=gi(2)
      ai(3,m)= ajdj1*gi(2)+dj1*gi(3)
      gi23m  =gi(2)-gi(3)
      ai(4,m)=aaj*gi23m
c
c
      ai(5,m)=gi(4)
      ai(6,m)= ajdj1*gi(4)+dj1*gi(5)
      gi45m  =gi(4)-gi(5)
      ai(7,m)=aaj*gi45m
c
c
c
c
      return
      end
      subroutine obaa
c
c        obaa computes the propagator and the cutoff of one         
c        ob-exchange as a function of angle theta;
c        thus, it provides part of the integrand used in the
c        integration performed by obai
c
c
      implicit real*8 (a-h,o-z)
c
c
c        common block for all ob-subroutines
c
      common /cob/    vj(32,25),c(10,15),fff,ff,f(52),aa(96),ai(19,5),
     1                wnn(3),wdd(3),x,xx,y,yy,xy2,xxpyy,ex,ey,eem12,
     2                ez1,ez2,ct(96),wt(96),
     3                ic(10,15),ift(3),mint(3),maxt(3),nt,
     4                mge,mgg(12,3),mggo(12,3),ima(5,12,3),
     5                imaa(3),imea(3),ime,im,mc,m,mg,inter,ide,idde,
     6                indc(2,15),indpar(3),indxy
c
c         specifications for this common block
c
      logical indc,indxy,indpar
c
c
c
c        further specifications
      real*8 deltaq(96)
c
      save
c
c
c
c
c        transferred three-momentum squared
c        ----------------------------------
c
c
c
c
   60 if (indxy) go to 1000
      indxy=.true.
      do 65 i=1,nt
      xy2t=xy2*ct(i)
c
c
   65 deltaq(i)=xxpyy-xy2t
c     --------------------
c
c
c
c
c        meson propagator
c        ----------------
c        ----------------
c
c
c
c
c        c4 is the meson mass squared
c
 1000 c4=c(4,im)
c
c                                   
      do 1011 i=1,nt
 1011 aa(i)=wt(i)/(c4+deltaq(i))
c
c
c
c
c        cut-offs
c        --------
c        --------
c
c
c
c
      mi=4
      mm=5
c
c
  999 ityp=ic(mi,im)
      if (ityp.eq.0) go to 2000
c
c
c
c
c        cut-off of monopole/dipole type
c        *******************************
c
c
c
c
  100 c5=c(mm,im)
      c6=c(mm+1,im)
      nexp=ic(mi+2,im)
c
      do 105 i=1,nt
c
      aaa=c5/(c6+deltaq(i))
c     ---------------------
c
      do 105 ii=1,nexp
  105 aa(i)=aa(i)*aaa
c
c
      mi=mi+3
      mm=mm+2
      go to 999
c
c
c
c
 2000 return
      end
      subroutine legp (pj,pjm1,x,j)
c
c
c        subroutine legp   computes the legendre polynominals
c
      real*8 pj,pjm1,x,a,b
c
c
c
c        compute legendre polynominal for j equals zero
c
c
c
      pjm1=1.d0
      if (j.gt.0) go to 1
      pj=1.d0
      return
c
c
c
c        compute legendre polynominals for j equals one
c
c
c
    1 pj=x
      if (j.eq.1) return
c
c
c
c        compute legendre polynominal for j greater or equal two
c
c
c
      do 2 i=2,j
      a=x*pj
      b=a-pjm1
      pjm1=pj
    2 pj=-b/dfloat(i)+b+a
c
c
      return
      end
      subroutine gset(ax,bx,n,z,w)
c
c
c        this code has been obtained from the CERN computer library
c
c
      implicit real*8 (a-h,o-z)
c
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
      dimension     a(273),x(273),ktab(96)
      dimension z(2),w(2)
c
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c
c**** n=2
      data x(1)/0.577350269189626  d0/, a(1)/1.000000000000000  d0/
c**** n=3
      data x(2)/0.774596669241483  d0/, a(2)/0.555555555555556  d0/
      data x(3)/0.000000000000000  d0/, a(3)/0.888888888888889  d0/
c**** n=4
      data x(4)/0.861136311594053  d0/, a(4)/0.347854845137454  d0/
      data x(5)/0.339981043584856  d0/, a(5)/0.652145154862546  d0/
c**** n=5
      data x(6)/0.906179845938664  d0/, a(6)/0.236926885056189  d0/
      data x(7)/0.538469310105683  d0/, a(7)/0.478628670499366  d0/
      data x(8)/0.000000000000000  d0/, a(8)/0.568888888888889  d0/
c**** n=6
      data x(9)/0.932469514203152  d0/, a(9)/0.171324492379170  d0/
      data x(10)/0.661209386466265 d0/, a(10)/0.360761573048139 d0/
      data x(11)/0.238619186083197 d0/, a(11)/0.467913934572691 d0/
c**** n=7
      data x(12)/0.949107912342759 d0/, a(12)/0.129484966168870 d0/
      data x(13)/0.741531185599394 d0/, a(13)/0.279705391489277 d0/
      data x(14)/0.405845151377397 d0/, a(14)/0.381830050505119 d0/
      data x(15)/0.000000000000000 d0/, a(15)/0.417959183673469 d0/
c**** n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c**** n=9
      data x(20)/0.968160239507626 d0/, a(20)/0.081274388361574 d0/
      data x(21)/0.836031107326636 d0/, a(21)/0.180648160694857 d0/
      data x(22)/0.613371432700590 d0/, a(22)/0.260610696402935 d0/
      data x(23)/0.324253423403809 d0/, a(23)/0.312347077040003 d0/
      data x(24)/0.000000000000000 d0/, a(24)/0.330239355001260 d0/
c**** n=10
      data x(25)/0.973906528517172 d0/, a(25)/0.066671344308688 d0/
      data x(26)/0.865063366688985 d0/, a(26)/0.149451349150581 d0/
      data x(27)/0.679409568299024 d0/, a(27)/0.219086362515982 d0/
      data x(28)/0.433395394129247 d0/, a(28)/0.269266719309996 d0/
      data x(29)/0.148874338981631 d0/, a(29)/0.295524224714753 d0/
c**** n=11
      data x(30)/0.978228658146057 d0/, a(30)/0.055668567116174 d0/
      data x(31)/0.887062599768095 d0/, a(31)/0.125580369464905 d0/
      data x(32)/0.730152005574049 d0/, a(32)/0.186290210927734 d0/
      data x(33)/0.519096129206812 d0/, a(33)/0.233193764591990 d0/
      data x(34)/0.269543155952345 d0/, a(34)/0.262804544510247 d0/
      data x(35)/0.000000000000000 d0/, a(35)/0.272925086777901 d0/
c**** n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c**** n=13
      data x(42)/0.984183054718588 d0/, a(42)/0.040484004765316 d0/
      data x(43)/0.917598399222978 d0/, a(43)/0.092121499837728 d0/
      data x(44)/0.801578090733310 d0/, a(44)/0.138873510219787 d0/
      data x(45)/0.642349339440340 d0/, a(45)/0.178145980761946 d0/
      data x(46)/0.448492751036447 d0/, a(46)/0.207816047536889 d0/
      data x(47)/0.230458315955135 d0/, a(47)/0.226283180262897 d0/
      data x(48)/0.000000000000000 d0/, a(48)/0.232551553230874 d0/
c**** n=14
      data x(49)/0.986283808696812 d0/, a(49)/0.035119460331752 d0/
      data x(50)/0.928434883663574 d0/, a(50)/0.080158087159760 d0/
      data x(51)/0.827201315069765 d0/, a(51)/0.121518570687903 d0/
      data x(52)/0.687292904811685 d0/, a(52)/0.157203167158194 d0/
      data x(53)/0.515248636358154 d0/, a(53)/0.185538397477938 d0/
      data x(54)/0.319112368927890 d0/, a(54)/0.205198463721296 d0/
      data x(55)/0.108054948707344 d0/, a(55)/0.215263853463158 d0/
c**** n=15
      data x(56)/0.987992518020485 d0/, a(56)/0.030753241996117 d0/
      data x(57)/0.937273392400706 d0/, a(57)/0.070366047488108 d0/
      data x(58)/0.848206583410427 d0/, a(58)/0.107159220467172 d0/
      data x(59)/0.724417731360170 d0/, a(59)/0.139570677926154 d0/
      data x(60)/0.570972172608539 d0/, a(60)/0.166269205816994 d0/
      data x(61)/0.394151347077563 d0/, a(61)/0.186161000015562 d0/
      data x(62)/0.201194093997435 d0/, a(62)/0.198431485327111 d0/
      data x(63)/0.000000000000000 d0/, a(63)/0.202578241925561 d0/
c**** n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/
c**** n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/
c**** n=24
      data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
      data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
      data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
      data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
      data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
      data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
      data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
      data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
      data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
      data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
      data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
      data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=32
      data x(94)/0.997263861849481 d0/, a(94)/0.007018610009470 d0/
      data x(95)/0.985611511545268 d0/, a(95)/0.016274394730905 d0/
      data x(96)/0.964762255587506 d0/, a(96)/0.025392065309262 d0/
      data x(97)/0.934906075937739 d0/, a(97)/0.034273862913021 d0/
      data x(98)/0.896321155766052 d0/, a(98)/0.042835898022226 d0/
      data x(99)/0.849367613732569 d0/, a(99)/0.050998059262376 d0/
      data x(100)/0.794483795967942d0/, a(100)/0.058684093478535d0/
      data x(101)/0.732182118740289d0/, a(101)/0.065822222776361d0/
      data x(102)/0.663044266930215d0/, a(102)/0.072345794108848d0/
      data x(103)/0.587715757240762d0/, a(103)/0.078193895787070d0/
      data x(104)/0.506899908932229d0/, a(104)/0.083311924226946d0/
      data x(105)/0.421351276130635d0/, a(105)/0.087652093004403d0/
      data x(106)/0.331868602282127d0/, a(106)/0.091173878695763d0/
      data x(107)/0.239287362252137d0/, a(107)/0.093844399080804d0/
      data x(108)/0.144471961582796d0/, a(108)/0.095638720079274d0/
      data x(109)/0.048307665687738d0/, a(109)/0.096540088514727d0/
c**** n=40
      data x(110)/0.998237709710559d0/, a(110)/0.004521277098533d0/
      data x(111)/0.990726238699457d0/, a(111)/0.010498284531152d0/
      data x(112)/0.977259949983774d0/, a(112)/0.016421058381907d0/
      data x(113)/0.957916819213791d0/, a(113)/0.022245849194166d0/
      data x(114)/0.932812808278676d0/, a(114)/0.027937006980023d0/
      data x(115)/0.902098806968874d0/, a(115)/0.033460195282547d0/
      data x(116)/0.865959503212259d0/, a(116)/0.038782167974472d0/
      data x(117)/0.824612230833311d0/, a(117)/0.043870908185673d0/
      data x(118)/0.778305651426519d0/, a(118)/0.048695807635072d0/
      data x(119)/0.727318255189927d0/, a(119)/0.053227846983936d0/
      data x(120)/0.671956684614179d0/, a(120)/0.057439769099391d0/
      data x(121)/0.612553889667980d0/, a(121)/0.061306242492928d0/
      data x(122)/0.549467125095128d0/, a(122)/0.064804013456601d0/
      data x(123)/0.483075801686178d0/, a(123)/0.067912045815233d0/
      data x(124)/0.413779204371605d0/, a(124)/0.070611647391286d0/
      data x(125)/0.341994090825758d0/, a(125)/0.072886582395804d0/
      data x(126)/0.268152185007253d0/, a(126)/0.074723169057968d0/
      data x(127)/0.192697580701371d0/, a(127)/0.076110361900626d0/
      data x(128)/0.116084070675255d0/, a(128)/0.077039818164247d0/
      data x(129)/0.038772417506050d0/, a(129)/0.077505947978424d0/
c**** n=48
      data x(130)/0.998771007252426d0/, a(130)/0.003153346052305d0/
      data x(131)/0.993530172266350d0/, a(131)/0.007327553901276d0/
      data x(132)/0.984124583722826d0/, a(132)/0.011477234579234d0/
      data x(133)/0.970591592546247d0/, a(133)/0.015579315722943d0/
      data x(134)/0.952987703160430d0/, a(134)/0.019616160457355d0/
      data x(135)/0.931386690706554d0/, a(135)/0.023570760839324d0/
      data x(136)/0.905879136715569d0/, a(136)/0.027426509708356d0/
      data x(137)/0.876572020274247d0/, a(137)/0.031167227832798d0/
      data x(138)/0.843588261624393d0/, a(138)/0.034777222564770d0/
      data x(139)/0.807066204029442d0/, a(139)/0.038241351065830d0/
      data x(140)/0.767159032515740d0/, a(140)/0.041545082943464d0/
      data x(141)/0.724034130923814d0/, a(141)/0.044674560856694d0/
      data x(142)/0.677872379632663d0/, a(142)/0.047616658492490d0/
      data x(143)/0.628867396776513d0/, a(143)/0.050359035553854d0/
      data x(144)/0.577224726083972d0/, a(144)/0.052890189485193d0/
      data x(145)/0.523160974722233d0/, a(145)/0.055199503699984d0/
      data x(146)/0.466902904750958d0/, a(146)/0.057277292100403d0/
      data x(147)/0.408686481990716d0/, a(147)/0.059114839698395d0/
      data x(148)/0.348755886292160d0/, a(148)/0.060704439165893d0/
      data x(149)/0.287362487355455d0/, a(149)/0.062039423159892d0/
      data x(150)/0.224763790394689d0/, a(150)/0.063114192286254d0/
      data x(151)/0.161222356068891d0/, a(151)/0.063924238584648d0/
      data x(152)/0.097004699209462d0/, a(152)/0.064466164435950d0/
      data x(153)/0.032380170962869d0/, a(153)/0.064737696812683d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/
c**** n=80
      data x(186)/0.999553822651630d0/, a(186)/0.001144950003186d0/
      data x(187)/0.997649864398237d0/, a(187)/0.002663533589512d0/
      data x(188)/0.994227540965688d0/, a(188)/0.004180313124694d0/
      data x(189)/0.989291302499755d0/, a(189)/0.005690922451403d0/
      data x(190)/0.982848572738629d0/, a(190)/0.007192904768117d0/
      data x(191)/0.974909140585727d0/, a(191)/0.008683945269260d0/
      data x(192)/0.965485089043799d0/, a(192)/0.010161766041103d0/
      data x(193)/0.954590766343634d0/, a(193)/0.011624114120797d0/
      data x(194)/0.942242761309872d0/, a(194)/0.013068761592401d0/
      data x(195)/0.928459877172445d0/, a(195)/0.014493508040509d0/
      data x(196)/0.913263102571757d0/, a(196)/0.015896183583725d0/
      data x(197)/0.896675579438770d0/, a(197)/0.017274652056269d0/
      data x(198)/0.878722567678213d0/, a(198)/0.018626814208299d0/
      data x(199)/0.859431406663111d0/, a(199)/0.019950610878141d0/
      data x(200)/0.838831473580255d0/, a(200)/0.021244026115782d0/
      data x(201)/0.816954138681463d0/, a(201)/0.022505090246332d0/
      data x(202)/0.793832717504605d0/, a(202)/0.023731882865930d0/
      data x(203)/0.769502420135041d0/, a(203)/0.024922535764115d0/
      data x(204)/0.744000297583597d0/, a(204)/0.026075235767565d0/
      data x(205)/0.717365185362099d0/, a(205)/0.027188227500486d0/
      data x(206)/0.689637644342027d0/, a(206)/0.028259816057276d0/
      data x(207)/0.660859898986119d0/, a(207)/0.029288369583267d0/
      data x(208)/0.631075773046871d0/, a(208)/0.030272321759557d0/
      data x(209)/0.600330622829751d0/, a(209)/0.031210174188114d0/
      data x(210)/0.568671268122709d0/, a(210)/0.032100498673487d0/
      data x(211)/0.536145920897131d0/, a(211)/0.032941939397645d0/
      data x(212)/0.502804111888784d0/, a(212)/0.033733214984611d0/
      data x(213)/0.468696615170544d0/, a(213)/0.034473120451753d0/
      data x(214)/0.433875370831756d0/, a(214)/0.035160529044747d0/
      data x(215)/0.398393405881969d0/, a(215)/0.035794393953416d0/
      data x(216)/0.362304753499487d0/, a(216)/0.036373749905835d0/
      data x(217)/0.325664370747701d0/, a(217)/0.036897714638276d0/
      data x(218)/0.288528054884511d0/, a(218)/0.037365490238730d0/
      data x(219)/0.250952358392272d0/, a(219)/0.037776364362001d0/
      data x(220)/0.212994502857666d0/, a(220)/0.038129711314477d0/
      data x(221)/0.174712291832646d0/, a(221)/0.038424993006959d0/
      data x(222)/0.136164022809143d0/, a(222)/0.038661759774076d0/
      data x(223)/0.097408398441584d0/, a(223)/0.038839651059051d0/
      data x(224)/0.058504437152420d0/, a(224)/0.038958395962769d0/
      data x(225)/0.019511383256793d0/, a(225)/0.039017813656306d0/
c**** n=96
      data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
      data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
      data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
      data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
      data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
      data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
      data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
      data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
      data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
      data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
      data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
      data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
      data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
      data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
      data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
      data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
      data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
      data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
      data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
      data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
      data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
      data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
      data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
      data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
      data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
      data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
      data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
      data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
      data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
      data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
      data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
      data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
      data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
      data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
      data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/
      data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/
      data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/
      data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/
      data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/
      data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/
      data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/
      data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/
      data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/
      data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
      data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
      data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
      data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
      data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
c
c
c-----test n
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c
  100 zn=n
      write(kwrite,200) zn
  200 format(/////' error in gset. n has the non-permissible value',
     1e11.3/' execution terminated.')
      stop
      end
