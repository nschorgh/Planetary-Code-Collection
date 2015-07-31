      subroutine conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhoc,
     &     Fgeotherm,Fsurf)
C***********************************************************************
C   conductionT:  program to calculate the diffusion of temperature 
C                 into the ground with prescribed surface temperature 
C                 and variable thermal properties on irregular grid
C   Crank-Nicholson scheme, flux conservative
C
C   Eqn: rhoc*T_t = (k*T_z)_z 
C   BC (z=0): T=T(t)
C   BC (z=L): heat flux = Fgeotherm
C
C   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
C   rhoc = rho*c  heat capacity per volume [J m^-3]  VECTOR
C   ti and rhoc are not allowed to vary in the layers immediately adjacent
C               to the surface or the bottom
C   T = vertical temperature profile [K]
C   Tsurf, Tsurfp1 = surface temperatures at times n and n+1
C
C   Grid: surface is at z=0
C         T(1) is at z(1); ...; T(i) is at z(i)
C         k(i) is midway between z(i-1) and z(i)
C         rhoc(i) is midway between z(i-1) and z(i)
C***********************************************************************

      implicit none
      integer NMAX
      parameter (NMAX=1000)

      integer nz, i
      real*8 z(NMAX), dt, T(NMAX), Tsurf, Tsurfp1, ti(NMAX), rhoc(NMAX)
      real*8 alpha(NMAX), k(NMAX), gamma(NMAX), buf, Fgeotherm, Fsurf
      real*8 a(NMAX), b(NMAX), c(NMAX), r(NMAX)
C     set some constants
      do i=1,nz
         k(i) = ti(i)**2/rhoc(i) ! thermal conductivity
      enddo
      alpha(1) = k(2)*dt/rhoc(1)/(z(2)-z(1))/z(2) 
      gamma(1) = k(1)*dt/rhoc(1)/z(1)/z(2) 
      do i=2,nz-1
         buf=dt/(z(i+1)-z(i-1))
         alpha(i) = k(i+1)*buf*2./(rhoc(i)+rhoc(i+1))/(z(i+1)-z(i))
         gamma(i) = k(i)*buf*2./(rhoc(i)+rhoc(i+1))/(z(i)-z(i-1))
      enddo
      buf=dt/(z(nz)-z(nz-1))**2
      gamma(nz) = k(nz)*buf/(rhoc(nz)+rhoc(nz)) ! assumes rhoc(nz+1)=rhoc(nz)

C     elements of tridiagonal matrix
      do i=1,nz
         a(i) = -gamma(i)   !  a(1) is not used
         b(i) = 1. + alpha(i) + gamma(i)
         c(i) = -alpha(i)   !  c(nz) is not used
      enddo
      b(nz) = 1. + gamma(nz)

C     Set RHS         
      r(1)= alpha(1)*T(2) + (1.-alpha(1)-gamma(1))*T(1)
     &     + gamma(1)*(Tsurf+Tsurfp1)
      do i=2,nz-1
         r(i) = gamma(i)*T(i-1) + (1.-alpha(i)-gamma(i))*T(i)
     &        + alpha(i)*T(i+1)
      enddo
      r(nz) = gamma(nz)*T(nz-1) + (1.-gamma(nz))*T(nz) +
     &     dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1)) ! assumes rhoc(nz+1)=rhoc(nz)

C     Solve for T at n+1
      call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion

      Fsurf = -k(1)*(T(1)-Tsurfp1)/z(1) ! heat flux into surface
      end
