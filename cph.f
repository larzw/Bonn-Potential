c        call phases
      program cph
      implicit real*8 (a-h,o-z)
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c        the parameter n has to be larger than or equal to the number
c        of gauss points (n) used in phases for the matrix inversion.
      parameter (n=64)
      parameter (n6=6*n,na=2*(n+1)*2*(n+1),naa=6*n/2*(n+1))
      dimension vv(n6),s(n),u(n),a(na),b(na),aa(naa),qq(n),eq(n)
      external bonn
      open (unit=5,file='dphbonnbpv.d')
      open (unit=6,file='phbonnbpv.d')
      open (unit=7,file='rma.d')
      kread=5
      kwrite=6
      kpunch=7
      call phases (bonn,vv,s,u,a,b,aa,qq,eq)
      end
