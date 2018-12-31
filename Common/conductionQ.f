      subroutine conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhoc,emiss,Tsurf,
     &     Fgeotherm,Fsurf)
C***********************************************************************
C   conductionQ:  program to calculate the diffusion of temperature 
C                 into the ground and thermal emission at the surface 
C                 with variable thermal properties on irregular grid
C   Crank-Nicolson scheme, flux conservative
C                          uses Samar's radiation formula
C   Eqn: rhoc*T_t = (k*T_z)_z 
C   BC (z=0): Q(t) + kT_z = em*sig*T^4
C   BC (z=L): heat flux = Fgeotherm
C
C   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
C   rhoc = rho*c  VECTOR where rho=density [kg m^-3] and 
C                              c=specific heat [J K^-1 kg^-1]
C   ti and rhoc are not allowed to vary in the layers immediately adjacent
C               to the surface or the bottom
C   T = vertical temperature profile [K] (output)
C   Qn,Qnp1 = net solar insolation at time steps n and n+1 [Watts m^-2]
C   emiss = emissivity
C   Tsurf = surface temperature [K]  (output)
C   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
C   Fsurf = heat flux at surface [W/m^2]  (output)   
C
C   Grid: surface is at z=0
C         z(0)=0, z(2)=3*z(1), i.e., the top layer has half the width
C         T(1) is at z(1); ...; T(i) is at z(i)
C         k(i) is midway between z(i-1) and z(i)
C         rhoc(i) is midway between z(i-1) and z(i)
C     
C   originally written by Samar Khatiwala, 2001
C   extended to variable thermal properties
C         and irregular grid by Norbert Schorghofer
C***********************************************************************

      implicit none
      integer NMAX
      real*8 sigSB
      parameter (NMAX=1000)
      parameter (sigSB=5.6704d-8)

      integer nz, i
      real*8 z(NMAX), dt, Qn, Qnp1, T(NMAX), ti(NMAX),rhoc(NMAX)
      real*8 emiss, Tsurf, Fgeotherm, Fsurf
      real*8 k(NMAX), k1, alpha(NMAX), gamma(NMAX), Tr
      real*8 a(NMAX), b(NMAX), c(NMAX), r(NMAX)
      real*8 arad, brad, ann, annp1, bn, buf, dz, beta


C     set some constants
      do i=1,nz
         k(i) = ti(i)**2/rhoc(i) ! thermal conductivity
      enddo
      dz=2.*z(1)
      beta = dt/rhoc(1)/(2.*dz**2)   ! assumes rhoc(0)=rhoc(1)
      alpha(1) = beta*k(2)
      gamma(1) = beta*k(1)
      do i=2,nz-1
         buf = dt/(z(i+1)-z(i-1))
         alpha(i) = 2.*k(i+1)*buf/(rhoc(i)+rhoc(i+1))/(z(i+1)-z(i))
         gamma(i) = 2.*k(i)*buf/(rhoc(i)+rhoc(i+1))/(z(i)-z(i-1))
      enddo
      buf=dt/(z(nz)-z(nz-1))**2
      gamma(nz) = k(nz)*buf/(2*rhoc(nz)) ! assumes rhoc(nz+1)=rhoc(nz)

      k1=k(1)/dz

C     elements of tridiagonal matrix
      do i=1,nz
         a(i) = -gamma(i)   !  a(1) is not used
         b(i) = 1. + alpha(i) + gamma(i) !  b(1) has to be reset at every timestep
         c(i) = -alpha(i)   !  c(nz) is not used
      enddo
      b(nz) = 1. + gamma(nz)
      
C     Emission
      Tr = Tsurf                !   'reference' temperature
      arad = -3.*emiss*sigSB*Tr**4
      brad = 2.*emiss*sigSB*Tr**3
      ann = (Qn-arad)/(k1+brad)
      annp1 = (Qnp1-arad)/(k1+brad)
      bn = (k1-brad)/(k1+brad)
      b(1) = 1. + alpha(1) + gamma(1) - gamma(1)*bn

C     Set RHS         
      r(1) = gamma(1)*(annp1+ann) + 
     &     (1.-alpha(1)-gamma(1)+gamma(1)*bn)*T(1) + alpha(1)*T(2)
      do i=2,nz-1
         r(i) = gamma(i)*T(i-1) + (1.-alpha(i)-gamma(i))*T(i)
     &        + alpha(i)*T(i+1)
      enddo
      r(nz) = gamma(nz)*T(nz-1) + (1.-gamma(nz))*T(nz)
     &     + dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1)) ! assumes rhoc(nz+1)=rhoc(nz)

C     Solve for T at n+1
      call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion
      
      Tsurf = 0.5*(annp1 + bn*T(1) + T(1)) ! (T0+T1)/2
      Fsurf = - k(1)*(T(1)-Tsurf)/z(1) ! heat flux into surface
      end
