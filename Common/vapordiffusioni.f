      subroutine vapordiffusioni(nz,z,dt,D,p,p0,Tn,Tnp1,rhof,T0,D0)
C***********************************************************************
C   vapordiffusion: diffusion and sublimation of water vapor in 
C                   sub-surface on irregular grid, with spatially
C                   variable diffusivity and finite void space 
C   Eqn: ((p/T)*(1-f) + k/mu*rhof)_t = (D*(1-f)*(p/T)_z)_z
C        with filling fraction f = rhof/rhoice
C   upper BC (z=0): prescribed partial pressure, permeable
C   lower BC (z=L): no flux of water vapor
C
C   D = diffusion constant without porosity (D12/tortuosity) [m^2/s]
C   p = partial pressure of H2O [Pa]
C   p0 = partial pressure of H2O at surface [Pa]
C   Tn = temperature profile at time step n  [K]
C   Tnp1 = temperature profile at time step n+1  [K]
C   rhof = density of free ice in space not filled by regolith [kg/m^3]
C   T0 = temperature at surface [K]
C   D0 = diffusion constant at surface (analogous to D) [m^2/s]
C
C   Grid: surface is at z=0
C         T(1) is at z(1); T(i) is at z(i); 
C         2*z(1)=z(2)-z(1) and arbitrary otherwise
C
C   Written by Norbert Schorghofer, ~2003
C***********************************************************************

      implicit none
      integer NMAX
      parameter (NMAX=1000)

      real*8 p(NMAX),Tn(NMAX),Tnp1(NMAX),rhof(NMAX)
      integer nz, j
      real*8 dt, z(NMAX), p0, T0, D(NMAX), kmu, rhs(NMAX), pfrost
      real*8 D0, a(NMAX), icedensity, a0, c, r(NMAX)
!      real*8 gamma, fp(NMAX), fT(NMAX)
      real*8 hp,hm,h2,cA0,cB0,cB1,cB2,cC0, rescale
      parameter (kmu = 8314.5/18.015, rescale=1.)
      parameter (icedensity = 927.*rescale) ! max. density of ice
      parameter (c=0.5)    ! 0<c<=0.5,  c=0.5 is most accurate
!      parameter (gamma=kmu/0.4)    !with adsorption
!      parameter (gamma=0.)    !no adsorption
      real*8 psv
      external psv

      do j=1,nz
         a(j) = D(j)*(1.-rhof(j)/icedensity)  ! constriction
         r(j) = p(j)/Tn(j)
      enddo
      a0 = D0*1.  ! phi0=1.
      rhs(1) = a(1)*(8./3.*p0/T0-4.*r(1)+4./3.*r(2)) +
     &     (-4./3.*a0+a(1)+1./3.*a(2))*
     &     (-4./3.*p0/T0+r(1)+1./3.*r(2))
      rhs(1) = rhs(1)/(z(2)-z(1))**2
      do j=2,nz-1
         hp=z(j+1)-z(j)
         hm=z(j)-z(j-1)
         h2=z(j+1)-z(j-1)
C        general corner-free
         cB0=-(2.*c+(1.-2.*c)*hp/hm)/(hm*hp)
         cA0=2.*(c-1.)/(hp*h2)
         cC0=(-1.+(1.-2.*c)*hp/hm)/(hm*h2)
         ! cA1=2.*(1.-c)/(hp*h2) = -cA0
         cB1=2.*c/(hp*h2)
         cB2=(1.+(1.-2.*c)*hp/hm)/(hm*h2)
         ! cC2=(1.+(2.*c-1.)*hp/hm)/(hm*h2) = -cC0
         rhs(j) = 
     &        a(j+1)*cA0* (-r(j+1) + r(j)) +
     &        a(j)* (cB1*r(j+1) + cB0*r(j) + cB2*r(j-1)) +
     &        a(j-1)*cC0* (r(j) - r(j-1))
      enddo
      !rhs(nz) = 2.*a(nz)*(r(nz-1)-r(nz))/(z(nz)-z(nz-1))**2
      !rhs(nz) = 0.
      rhs(nz) = 2.*a(nz)*(r(nz-1)-r(nz))/(z(nz)-z(nz-1))**2/nz  ! trick

!      call adsorption(p,Tn,fp,fT,nz)  ! slow, comment out when not needed
!      print *,rhs(1),gamma*fp(1)*p(1),gamma*fT(1)*Tn(1)

      do j=1,nz
         rhs(j) = dt*rhs(j) + r(j)*(1.-rhof(j)/icedensity) + 
     &        kmu*rhof(j) 
!     &        + gamma*fp(j)*p(j) - gamma*fT(j)*(Tnp1(j)-Tn(j))
      enddo

      do j=1,nz
         !p(j)=rhs(j)*Tnp1(j)/(1.+Tnp1(j)*gamma*fp(j))
         p(j)=rhs(j)*Tnp1(j)
         pfrost=psv(Tnp1(j))
         if (p(j)>pfrost) then
            p(j) = pfrost
            rhof(j) = (rhs(j)-pfrost/Tnp1(j))/
!            rhof(j) = (rhs(j)-pfrost/Tnp1(j)-gamma*fp(j)*pfrost)/
     &           (kmu-pfrost/(Tnp1(j)*icedensity))
            if (rhof(j)>icedensity) rhof(j)=icedensity
         else
            rhof(j) = 0.
         endif
         if (p(j)<0..or.rhof(j)<0..or.rhof(j)>icedensity) then
            print *,'unphysical result (numerical instability)',
     &           j,p(j),rhof(j)
C           stop
         endif
      enddo

      end

