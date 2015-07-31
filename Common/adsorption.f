      subroutine adsorption(p,T,fp,fT,nz)
C***********************************************************************
C   adsorption: calculates the amount adsorbed and returns the 
C               derivatives with respect temperature and pressure
C
C   p = partial pressure of H2O [Pa]
C   T = temperature [K]
C   f = mass density of adsorbed H20 [kg/m^3]
C   fp = partial derivative of f with respect to pressure
C   fT = partial derivative of T with respect to temperature
C   nz = number of grid points
C***********************************************************************
      implicit none
      integer NMAX
      parameter (NMAX=1000)

      real*8 p(NMAX), T(NMAX), f, fp(NMAX), fT(NMAX)
      integer nz, j

c     Jakosky et al., Icarus 130, 87-95 (1997) based on Zent & Quinn
      real*8 K, K0, e, nu, vm
      parameter (K0=1.57e-8, e=2573.9, nu=0.48) 
      parameter (vm=1500.*1.e5*2.84e-7) ! J/mol   
!     (density of regolith x specific surface area x monolayer mass density)
!           (kg/m^3)       x        (m^2/kg)       x   (kg/m^2/monolayer)

c     Zent et al., basalt
!      real*8 rhor, gamma, B0, delta
!      parameter (rhor=2000., B0=0.51)  ! rhor = density of regolith
!      parameter (gamma=8.364e-16, delta=-2679.8)

      do j=1,nz
c        Jakosky et al. '97
         K = K0*exp(e/T(j))
         f = vm*(K*p(j)/(1.+K*p(j)))**nu
         fp(j) = f/p(j)*nu/(1.+K*p(j)) 
         fT(j) = -f*nu*e/T(j)**2/(1.+K*p(j))
c         Zent et al.
!         f = 1.e-5*rhor*(gamma*p(j))**B0*exp(-delta/T(j))
!         fp(j) = B0*f/p(j)      ! fails for p=0
!         fT(j) = delta*f/T(j)**2
      enddo
      end




      subroutine adsorption2(p,T,f,nz,psurf,Tsurf,fsurf)
C***********************************************************************
C   adsorption2: returns the amount of adsorbed H2O
C
C   (use exactly the same isotherms as in subroutine adsorption)
C
C   p = partial pressure of H2O [Pa]
C   T = temperature [K]
C   f = mass density of adsorbed H20 [kg/m^3]
C   nz = number of grid points
C   ?surf = quantities at surface
C***********************************************************************
      implicit none
      integer NMAX
      parameter (NMAX=1000)

      real*8 p(NMAX), T(NMAX), f(NMAX), psurf, Tsurf, fsurf
      integer nz, j

c     Jakosky et al., Icarus 130, 87-95 (1997) based on Zent & Quinn
      real*8 K, K0, e, nu, vm
      parameter (K0=1.57e-8, e=2573.9, nu=0.48) 
      parameter (vm=1500.*1.e5*2.84e-7) ! J/mol

c     Zent et al., basalt
!      real*8 rhor, gamma, B0, delta
!      parameter (rhor=2000., B0=0.51)
!      parameter (gamma=8.364e-16, delta=-2679.8)

      K = K0*exp(e/Tsurf)
      fsurf = vm*(K*psurf/(1.+K*psurf))**nu
      do j=1,nz
c        Jakosky et al. '97
         K = K0*exp(e/T(j))
         f(j) = vm*(K*p(j)/(1.+K*p(j)))**nu
c        Zent et al.
!        f(j) = 1.e-5*rhor*(gamma*p(j))**B0*exp(-delta/T(j))
      enddo
      end


