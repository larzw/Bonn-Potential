      subroutine phases (pot,vv,s,u,a,b,aa,qq,eq)
c
c****  in this version, the on-shell K-matrix is written out.
c
c
c        phases computes the r-matrix (or 'k-matrix') and  
c        the phase-shifts of nucleon-nucleon scattering
c        for a given nn-potential pot.
c        the calculations are performed in momentum space.
c        this package contains all subroutines needed, except for
c        subroutine gset that is contained in the package bonn.
c        furthermore, a potential subroutine is required, e.g. bonn.
c        an 'external' statement in the calling program has to specify
c        the name of the potential subroutine to be applied that will
c        replace 'pot' in this program.
c
c
c        this version of the code is to be published in "computational 
c        nuclear physics", vol. II, koonin, langanke, maruhn, eds.
c        (springer, heidelberg).
c
c
c        author: r. machleidt
c                department of physics
c                university of idaho
c                moscow, idaho 83843
c                u. s. a.
c                e-mail: machleid@phys.uidaho.edu
c
c                formerly:
c                institut fuer theoretische kernphysik bonn
c                nussallee 14-16
c                d-5300 bonn, w. germany
c
c
      implicit real*8 (a-h,o-z)
      external pot
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments of the potential subroutine pot being called in this
c        program
c
      common /cpot/   v(6),xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
c        xmev and ymev are the final and initial relative momenta,
c        respectively, in units of mev.
c        v is the potential in units of mev**(-2).
c        if heform=.true., v must contain the 6 helicity matrix
c        elements associated with one j in the following order:
c        0v, 1v, 12v, 34v, 55v, 66v (helicity formalism).
c        if heform=.false., v must contain the six lsj-state matrix
c        elements associated with one j in the following order:
c        0v, 1v, v++, v--, v+-, v-+ (lsj formalism).
c        j is the total angular momentum. there is essentially
c        no upper limit for j.
c
c        specifications for these arguments
      logical heform,sing,trip,coup,endep
c
c        further specifications
c
c        all matrices are stored columnwise in vector arrays.
c        see main program calling phases, cph, for the dimensions
c        of the following arrays
      real*8 vv(1),s(1),u(1),a(1),b(1),aa(1),qq(1),eq(1)
      real*8 q(97)
      real*8 delta(5)
      data delta/5*0.d0/
      real*8 rb(2)
      real*8 r(6)
      data uf/197.3286d0/
      data pih/1.570796326794897d0/
      real*4 ops
      data ops/1.e-15/
      data nalt/-1/
c        the following dimension for at most 40 elabs
      data memax/41/
      real*8 elab(41),q0(40),q0q(40),eq0(40)
c     
      integer nj(20)
c
      character*4 name(3),nname(15)
      character*1 state(2),state3
      character*1 multi(4)
      data multi/'1',3*'3'/
      integer ldel(4)
      data ldel/0,0,-1,1/
      character*1 spd(50)
      data spd/'s','p','d','f','g','h','i','k','l','m','n',
     1'o','q','r','t','u','v','w','x','y','z',29*' '/
      character*1 chars,chart
      data chars/'s'/,chart/'t'/
      character*1 lab1,lab2
      character*4 blanks
      data blanks/'    '/
      logical indbrn
      logical indrma
      logical indqua
      logical indwrt
      logical indpts
      logical inderg
      data indbrn/.false./
      data indrma/.false./
      data indqua/.false./
      data indwrt/.false./
      data indpts/.false./
      data inderg/.false./
      logical indj
c
c
c
c
10000 format (2a4,a2,20i3)
10001 format (1h ,2a4,a2,20i3)
10002 format (2a4,a2,6f10.4)
10003 format (1h ,2a4,a2,6f10.4)
10004 format (//' input-parameters for phases'/1h ,27(1h-))
10005 format (//' transformed gauss points and weights for c =',f8.2,
     1'  and n =',i3//' points')
10006 format (7x,4f15.4)
10007 format (/' weights')
10008 format (2a4,a2,15a4)
10009 format (1h ,2a4,a2,15a4)
10013 format (///' error in phases. matrix inversion error index =',
     1 i5///)
10014 format (1h ,2a4,a2,6f10.4)
10016 format (/' low energy parameters    a',a1,' =',
     1 f10.4,'    r',a1,' =',f10.4)
10050 format (//' elab(mev)',5x,4(2a1,i2,10x),' e',i2/)
10051 format (1h ,f8.2,5e14.6)
10052 format (//' elab(mev)',5x,'1s0',25x,'3p0'/)
10053 format (1h1//' p h a s e - s h i f t s (radians)'/1h ,33(1h-))
10054 format (1h1//' p h a s e - s h i f t s (degrees)'/1h ,33(1h-))
10055 format (1h ,f8.2,5f14.6)
10110 format (1h ,a4,2x,2a1,'-',a1,i1,3d20.12)
10111 format (i3)
10112 format (' elab (mev) ',f10.4)
10113 format (4d20.12)
c
c
c
c
c        read and write input parameters for this program
c        ------------------------------------------------
c
c
      write (kwrite,10004)
c
c        read and write comment
      read  (kread, 10008) name,nname
      write (kwrite,10009) name,nname
c
c        read and write the range for the total angular momentum j
c        for which phase shifts are to be calculated. there are
c        essentially no limitations for j.
      read  (kread ,10000) name,jb,je
      write (kwrite,10001) name,jb,je
c
c        the born approximation will be used for j.ge.jborn.
      read  (kread ,10000) name,jborn
      write (kwrite,10001) name,jborn
      jb1=jb+1
      je1=je+1
      jee1=je1
      if (jee1.gt.20) jee1=20
c
c        read and write the number of gauss points to be used for the
c        matrix inversion; this is j-dependent.
      read  (kread ,10000) name,(nj(j1),j1=1,jee1)
      write (kwrite,10001) name,(nj(j1),j1=1,jee1)
c
c        ihef=0: lsj formalism is used (heform=.false.),
c        ihef.ne.0: helicity formalism is used (heform=.true.).
      read  (kread ,10000) name,ihef
      write (kwrite,10001) name,ihef
c
c        set ising=itrip=icoup=1, if you want that for each j
c        all possible states are considered.
      read  (kread ,10000) name,ising,itrip,icoup
      write (kwrite,10001) name,ising,itrip,icoup
c
c        iprop=1: non-relativistic propagator in scattering equation;
c        iprop=2: relativistic propagator.
      read  (kread ,10000) name,iprop
      write (kwrite,10001) name,iprop
c
c        iphrel=1: non-relativistic phase-relation;
c        iphrel=2: relativistic phase-relation.
      read  (kread ,10000) name,iphrel
      write (kwrite,10001) name,iphrel
c
c        c is the factor involved in the transformation of the
c        gauss points.
      read  (kread ,10002) name,c
      write (kwrite,10003) name,c
c
c        wn is the nucleon mass.
      read  (kread ,10002) name,wn
      write (kwrite,10003) name,wn
c
      do 1 k=1,memax
c        read and write the elab's for which phase shifts are to be
c        calculated. any number of elab's between 1 and 40 is permissible.
c        the last elab has to be zero and is understood as the end
c        of the list of elab's.
      read  (kread ,10002) name,elab(k)
      write (kwrite,10014) name,elab(k)
      if (elab(k).eq.0.d0) go to 2
    1 continue
    2 melab=k-1
c
c        irma=0:    the r-matrix is not written.
c        irma.ne.0: the r-matrix is written in terms of lsj states
c                   to unit kpunch,
c                   iqua=0: only the half off-shell r-matrix is written,
c                   iqua.ne.0: the (quadratic) fully off-shell r-matrix 
c                              is calculated and written.
      read  (kread ,10000) name,irma,iqua
      write (kwrite,10001) name,irma,iqua
c
c        ipoint=0: the transformed gauss-points and -weights are
c                  not printed;
c        ipoint.ne.0: ... are printed.
      read  (kread ,10000) name,ipoint
      write (kwrite,10001) name,ipoint
c
c        ideg=0: phase-shifts are printed in radians;
c        ideg=1: phase-shifts are printed in degrees.
      read  (kread ,10000) name,ideg
      write (kwrite,10001) name,ideg
c
c
c        prepare constants
c        -----------------
c
c
      sing=.false.
      trip=.false.
      coup=.false.
      if (ising.ne.0) sing=.true.
      if (itrip.ne.0) trip=.true.
      if (icoup.ne.0) coup=.true.
      heform=.false.
      if (ihef.ne.0) heform=.true.
      if (irma.ne.0) indrma=.true.
      if (irma.ne.0.and.iqua.ne.0) indqua=.true.
      if (ipoint.ne.0) indpts=.true.
c
      if(jb.le.1.and.elab(1).le.2.d0.and.elab(2).le.2.d0.and.melab.ge.2)
     1inderg=.true.
      wnh=wn*0.5d0
      wnq=wn*wn
      wp=wn*pih
      rd=90.d0/pih
      iideg=ideg+1
c**** label=blanks
c
c
c
c        prepare energies and on-shell momenta
c
      do 10 k=1,melab
      q0q(k)=wnh*elab(k)
      q0(k)=dsqrt(q0q(k))
   10 eq0(k)=dsqrt(q0q(k)+wnq)
c
c
c
c
c        loop of total angular momentum j
c        --------------------------------
c        --------------------------------
c
c
c
c
      do 2000 j1=jb1,je1
c
c
      indj=.false.
      j2=j1+1
      j=j1-1
      aj=dfloat(j)
      aj1=dfloat(j+1)
      a2j1=dfloat(2*j+1)
      d2j1=1.d0/a2j1
      aaj=dsqrt(aj*aj1)
c
c
      if (j.ge.jborn) indbrn=.true.
c
c
      if (j1.le.20) go to 105
      if (nalt.ge.1) go to 300
      n=16
      go to 115
c
c        number of gausspoints for this j
  105 if (nj(j1).eq.0) nj(j1)=nj(j)
      n=nj(j1)
      if (n.eq.nalt) go to 300
c
c
c        get gauss points and weights
c
  115 call gset (0.d0,1.d0,n,u,s)
c
      nalt=n
      n1=n+1
      n2=2*n1
      n3=3*n1
      n4=4*n1
      nx=n1
      if (indqua) nx=nx*nx
      nx2=2*nx
      nx2mn=nx2-n1
      nx2pn=nx2+n1
      nx4=4*nx
      nx4mn=nx4-n1
c
c        transform gauss points and weights
c
      do 201 i=1,n
      xx=pih*u(i)
c
c        transformed gauss point
      q(i)=dtan(xx)*c
      qq(i)=q(i)*q(i)
      eq(i)=dsqrt(qq(i)+wnq)
c
c        transformed gauss weight
      dc=1.d0/dcos(xx)
  201 s(i)=pih*c*dc*dc*s(i)
c
      if (.not.indqua) go to 205
c
      write (kpunch,10111) n
      write (kpunch,10113) (q(i),i=1,n)
c
  205 if (.not.indpts) go to 300
c
c        write gauss points and weights
c
      write (kwrite,10005) c,n
      write (kwrite,10006) (q(i),i=1,n)
c
      write (kwrite,10007)
      write (kwrite,10006) (s(i),i=1,n)
c
c
c
c
c
c
c        loop of elabs
c        -------------
c        -------------
c
c
c
c
  300 do 1000 k=1,melab
c
      q(n1)=q0(k)
c
c
      if (indrma) write (kpunch,10112) elab(k)
c
c
      if (indbrn.and..not.indqua) go to 500
c
c
c        check if right potential matrix does already exist
      if (indj) go to 500
      indj=.true.
c
c
c        compute potential matrix
c        ------------------------
c
c
      iii=0
      do 401 ix=1,n
c
c
      xmev=q(ix)
c
c
      do 401 iy=ix,n
c
c
      ymev=q(iy)
c
c
      call pot
c
c
      iaa=iii*6
      iii=iii+1
      do 401 iv=1,6
  401 aa(iv+iaa)=v(iv)
c
c
c        compute potential vector
c        ------------------------
c
c
  500 if (indbrn.and..not.indrma) go to 510
c
      ymev=q0(k)
c
c
      do 501 ix=1,n
c
c
      xmev=q(ix)
c
c
      call pot
c
c
      do 501 iv=1,6
      ivv=ix+(iv-1)*n
  501 vv(ivv)=v(iv)
c
c
c        compute potential element
c        -------------------------
c
c
  510 xmev=q0(k)
      ymev=q0(k)
c
c
      call pot
c
c
c
c        compute factor for the phase relation
c
      go to (601,602),iphrel
  601 wpq0=-wp*q0(k)
      go to 605
  602 wpq0=-pih*eq0(k)*q0(k)
  605 continue
c
c
      if (indbrn.and..not.indqua) go to 700
c
c
c        compute propagator
c        ------------------
c
c
      uq0=0.d0
      do 620 i=1,n
      sdq=s(i)/(qq(i)-q0q(k))
c
c        calculate propagator of lippmann-schwinger equation
c
      go to (621,622),iprop
  621 u(i)=sdq*qq(i)*wn
      go to 620
  622 u(i)=sdq*qq(i)*(eq(i)+eq0(k))*0.5d0
c
  620 uq0=uq0+sdq
c
c
      go to (631,632),iprop
  631 uq0=-uq0*q0q(k)*wn
      go to 700
  632 uq0=-uq0*q0q(k)*eq0(k)
c
c
c
c
c
c        build up matrix to be inverted
c        ------------------------------
c
c
  700 ni=0
      nii=0
      nv=n1
      mv=1
      if (indqua) mv=mv*n1
      ib=0
      eins=1.d0
c
c
      if (.not.sing) go to 720
      iv=1
      go to 770
c
c
  720 if (.not.trip.or.j.eq.0) go to 730
      iv=2
      go to 770
c
c
  730 if (.not.coup) go to 900
      iv=3
      if (j.eq.0) go to 770
      nv=n2
      mv=2
      if (indqua) mv=mv*n1
      go to 770
c
  740 if (j.eq.0) go to 800
      iv=4
      ib=n3
      ni=n1
      nii=n1
      go to 770
c
  750 iv=5
      ivi=6
      ib=n2
      ni=0
      nii=n1
      eins=0.d0
      go to 770
c
  760 iv=6
      ivi=5
      ib=n1
      ni=n1
      nii=0
c
c
c
c
  770 iii=0
      if (iv.le.4) ivi=iv
      igg=(iv-1)*n
      i1=(nii+n)*nv
      i2=(nii-1)*nv
c
c
      if (indbrn.and..not.indrma) go to 785
c
c
      do 780 i=1,n
      i3=i*nv
      i4=ni+i
c
c
      do 781 ii=i,n
      iaa=iii*6
      iii=iii+1
      i5=i2+i3+ni+ii
      i6=ivi+iaa
      i7=i2+i4+ii*nv
      i8=iv+iaa
      if (i.eq.ii) go to 782
c
c        matrix a
      a(i7)=aa(i8)*u(ii)
      a(i5)=aa(i6)*u(i)
      if (.not.indqua) go to 781
c
c        matrix b
      b(i7)=aa(i8)
      b(i5)=aa(i6)
      go to 781
c        diagonal element
  782 a(i7)=aa(i8)*u(i)+eins
      if (.not.indqua) go to 781
      b(i7)=aa(i8)
  781 continue
c
c        last column
      i9=i1+i4
      i10=i+igg
      a(i9)=vv(i10)*uq0
c        last row
      i11=i2+i3+ni+n1
      ivv=i+(ivi-1)*n
      a(i11)=vv(ivv)*u(i)
      if (.not.indqua) go to 783
      b(i9)=vv(i10)
      b(i11)=vv(ivv)
      go to 780
c
c        vector b
  783 b(ib+i)=vv(i+igg)
c
  780 continue
c
c
c        last element
      i12=i1+ni+n1
      a(i12)=v(iv)*uq0+eins
      if (.not.indqua) go to 785
      b(i12)=v(iv)
      go to 790
  785 b(ib+n1)=v(iv)
c
c
c
c
  790 go to (800,800,740,750,760,800),iv
c
c
c
c
c        invert matrix
c        -------------
c
c
c
c
  800 if (indbrn) go to 801
      call dgelg (b,a,nv,mv,ops,ier)
c
c
      if (ier.ne.0) write(kwrite,10013) ier
c
c
c
c
c        compute phase shifts
c        --------------------
c
c
c
c
  801 if (iv.gt.2.and.j.ne.0) go to 820
c
c        uncoupled cases
c
      delta(iv)=datan(b(nx)*wpq0)
c
c        prepare for effective range
      if (inderg.and.j.eq.0.and.iv.eq.1.and.k.le.2) 
     1rb(k)=q0(k)/(b(nx)*wpq0)
c
c
      if (.not.indrma) go to 810
c
c
c
c        write r-matrix
c
c
      state(1)=multi(iv)
      ispd=j1+ldel(iv)
      if (j.eq.0.and.iv.eq.3) ispd=2
      state(2)=spd(ispd)
      state3=state(2)
      if (indqua) write (kpunch,10110) label,state,state3,j
      do 805 i=n1,n1
      if (indqua) go to 804
c
c        write half off-shell r-matrix
      write (kpunch,10110) label,state,state3,j,q(i),q0(k),b(i)
      go to 805
c
c        write fully off-shell r-matrix (lower triangle)
  804 i1=(i-1)*n1
      write (kpunch,10113) (b(i1+ii),ii=i,n1)
  805 continue
c
c
c
c
  810 go to (720,730,900),iv
c
c
c        coupled cases
c
  820 if (heform) go to 822
c
c        calculate phase shifts from lsj-state r-matrix elements
c
      r0=b(nx2) 
      r1=b(nx2mn)-b(nx4)
      r2=b(nx2mn)+b(nx4)
      rr=-2.d0*r0/r1
c        epsilon
      delta(5)=datan(rr)/2.d0
      rr=r1*dsqrt(1.d0+rr*rr)
c        prepare for effective range
      if (inderg.and.j.eq.1.and.k.le.2) 
     1rb(k)=2.d0*q0(k)/(wpq0*(r2-rr))
c        delta minus
      delta(3)=datan((r2-rr)*wpq0*0.5d0)
c        delta plus
      delta(4)=datan((r2+rr)*wpq0*0.5d0)
      go to 824
c
c
c        calculate phase shifts from helicity-state r-matrix elements
c
  822 r0=b(nx2)
      r1=b(nx2mn)-b(nx4)
      r2=b(nx2mn)+b(nx4)
      rr=-2.d0*(aaj*r1+r0)/(r1-4.d0*aaj*r0)
c        epsilon
      delta(5)=datan(rr)/2.d0
      rr=(r1-4.d0*aaj*r0)*dsqrt(1.d0+rr*rr)*d2j1
c        prepare for effective range
      if (inderg.and.j.eq.1.and.k.le.2) 
     1rb(k)=2.d0*q0(k)/(wpq0*(r2-rr))
c        delta minus
      delta(3)=datan((r2-rr)*wpq0*0.5d0)
c        delta plus
      delta(4)=datan((r2+rr)*wpq0*0.5d0)
c
c        so far the delta(..) have been the blatt-biedenharn phase-
c        shifts, transform the delta(..) now into bar-phase-shifts
c        according to stapp et al., phys. rev. 105 (1957) 302.
  824 if (delta(5).eq.0.d0) go to 829
      pp=delta(3)+delta(4)
      pm=delta(3)-delta(4)
      d52=2.d0*delta(5)
      pm2=dsin(d52)*dsin(pm)
      delta(5)=0.5d0*dasin(pm2)
      pm1=dtan(2.d0*delta(5))/dtan(d52)
      pm1=dasin(pm1)
      delta(3)=0.5d0*(pp+pm1)
      delta(4)=0.5d0*(pp-pm1)
  829 continue
c
c
      if (.not.indrma) go to 900
c
c
c
      na=1
      if (indqua) na=n1
      do 835 i=1,na
      i1=(i-1)*n2
      do 835 ii=i,n1
      i3=i1+ii
      i4=nx2pn+i3
      i5=nx2+i3
      i6=n1+i3
      r(3)=b(i3)
      r(4)=b(i4)
      r(5)=b(i5)
      r(6)=b(i6)
      if (.not.heform) go to 837
c
c        in case of heform=.true., transform into lsj-form
      r34=(r(3)-r(4))*aaj
      r56=(r(5)+r(6))*aaj
      b(i6)=(aj1*r(3)+aj*r(4)-r56)*d2j1
      b(i3)=(aj*r(3)+aj1*r(4)+r56)*d2j1
      b(i5)=(r34+aj1*r(5)-aj*r(6))*d2j1
      b(i4)=(r34-aj*r(5)+aj1*r(6))*d2j1
      go to 835
c
c        in case of heform=.false., reorganize
  837 b(i6)=r(3)
      b(i3)=r(4)
      b(i5)=r(5)
      b(i4)=r(6)
  835 continue
c
c        write r-matrix
      ivx=0
      do 840 iv1=3,4
      do 840 iv2=3,4
      ivx=ivx+1
      state(1)=multi(iv1)
      state(2)=spd(j1+ldel(iv1))
      state3=spd(j1+ldel(iv2))
      go to (831,832,833,834),ivx
  831 ny=0
      go to 836
  832 ny=nx2pn
      go to 836
  833 ny=nx2
      go to 836
  834 ny=n1
c
c
  836 if (indqua) write (kpunch,10110) label,state,state3,j
      do 839 i=n1,n1
      if (indqua) go to 838
c
c        write half off-shell r-matrix
      ii1=ny+i
      write (kpunch,10110) label,state,state3,j,q(i),q0(k),b(ii1)
      go to 839
c
c        write fully off-shell r-matrix (lower triangle)
  838 i1=(i-1)*n2+ny
      write (kpunch,10113) (b(i1+ii),ii=i,n1)
  839 continue
  840 continue
c
c
c
c
  900 continue
c
c
c        write phase-shifts
c        ------------------
c
      if (indwrt) go to 921
      indwrt=.true.
      go to (931,932),iideg
  931 write (kwrite,10053)
      go to 933
  932 write (kwrite,10054)
  933 continue
c
  921 if (k.ne.1) go to 923
      if (j.ne.0) go to 922
      write (kwrite,10052)
      go to 923
  922 write (kwrite,10050) multi(1),spd(j1),j,
     1                     multi(2),spd(j1),j,
     2                     multi(2),spd(j),j,
     3                     multi(2),spd(j2),j,
     4                     j
  923 if (iideg.eq.1) go to 926
      do 925 iv=1,5
  925 delta(iv)=delta(iv)*rd
      write (kwrite,10055) elab(k),delta
      go to 1000
  926 write (kwrite,10051) elab(k),delta
c
c
c
c
 1000 continue
c        this has been the end of the elab loop
c
c
c
c
c        calculate and write low energy parameters
c        -----------------------------------------
c
c
      if (.not.inderg) go to 2000
      if (j.gt.1) go to 2000
      rb2=4.d0/wn*(rb(1)-rb(2))/(elab(1)-elab(2))*uf
      rb1=wn/4.*elab(1)*rb2/uf-rb(1)
      rb1=1./rb1*uf
      if (j.ne.0) go to 1090
      lab1=chars
      lab2=chars
      go to 1091
 1090 lab1=chart
      lab2=chart
 1091 write (kwrite,10016) lab1,rb1,lab2,rb2
c
c
c
c
 2000 continue
c        this has been the end of the j loop
c
c
      stop
      end
c*************************************************************
c
c name:      dgelg
c
c from:      programmbibliothek rhrz bonn, germany;   02/02/81
c            (free soft-ware)
c language:  fortran iv (fortran-77 compatible)
c
c purpose:
c
c to solve a general system of simultaneous linear equations.
c
c usage:   call dgelg(r,a,m,n,eps,ier)
c
c parameters:
c
c r:       double precision m by n right hand side matrix
c          (destroyed). on return r contains the solutions
c          of the equations.
c
c a:       double precision m by m coefficient matrix
c          (destroyed).
c
c m:       the number of equations in the system.
c
c n:       the number of right hand side vectors.
c
c eps:     single precision input constant which is used as
c          relative tolerance for test on loss of
c          significance.
c
c ier:     resulting error parameter coded as follows
c           ier=0  - no error,
c           ier=-1 - no result because of m less than 1 or
c                   pivot element at any elimination step
c                   equal to 0,
c           ier=k  - warning due to possible loss of signifi-
c                   cance indicated at elimination step k+1,
c                   where pivot element was less than or
c                   equal to the internal tolerance eps times
c                   absolutely greatest element of matrix a.
c
c remarks: (1) input matrices r and a are assumed to be stored
c              columnwise in m*n resp. m*m successive storage
c              locations. on return solution matrix r is stored
c              columnwise too.
c          (2) the procedure gives results if the number of equations m
c              is greater than 0 and pivot elements at all elimination
c              steps are different from 0. however warning ier=k - if
c              given indicates possible loss of significance. in case
c              of a well scaled matrix a and appropriate tolerance eps,
c              ier=k may be interpreted that matrix a has the rank k.
c              no warning is given in case m=1.
c
c method:
c
c solution is done by means of gauss-elimination with
c complete pivoting.
c
c
c author:         ibm, ssp iii
c
c**********************************************************************
      subroutine dgelg(r,a,m,n,eps,ier)
c
c
      dimension a(1),r(1)
      real*8 r,a,piv,tb,tol,pivi
c
c
c
c
      if(m)23,23,1
c
c     search for greatest element in matrix a
    1 ier=0
      piv=0.d0
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb=dabs(a(l))
      if(tb-piv)3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
c     a(i) is pivot element. piv contains the absolute value of a(i).
c
c
c     start elimination loop
      lst=1
      do 17 k=1,m
c
c     test on singularity
      if(piv)23,23,4
    4 if(ier)7,5,7
    5 if(piv-tol)6,6,7
    6 ier=k-1
    7 pivi=1.d0/a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
c     i+k is row-index, j+k column-index of pivot element
c
c     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
c
c     is elimination terminated
      if(k-m)9,18,18
c
c     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
c
c     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
c
c     save column interchange information
      a(lst)=j
c
c     element reduction and next pivot search
      piv=0.d0
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb=dabs(a(l))
      if(tb-piv)15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
c     end of elimination loop
c
c
c     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+.5d0
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
c
c
c     error return
   23 ier=-1
      return
      end
