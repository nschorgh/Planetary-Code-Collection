subroutine conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhoc,Fgeotherm,Fsurf)
!***********************************************************************
!   conductionT:  program to calculate the diffusion of temperature 
!                 into the ground with prescribed surface temperature 
!                 and variable thermal properties on irregular grid
!   Crank-Nicholson scheme, flux conservative
!
!   Eqn: rhoc*T_t = (k*T_z)_z 
!   BC (z=0): T=T(t)
!   BC (z=L): heat flux = Fgeotherm
!
!   nz = number of grid points
!   dt = time step
!   T = vertical temperature profile [K]  (in- and output)
!   Tsurf, Tsurfp1 = surface temperatures at times n and n+1  
!   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
!   rhoc = rho*c  heat capacity per volume [J m^-3 K^-1]  VECTOR
!   ti and rhoc are not allowed to vary in the layers immediately
!               adjacent to the surface or the bottom
!   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
!   Fsurf = heat flux at surface [W/m^2]  (output)
!
!   Grid: surface is at z=0
!         T(1) is at z(1); ...; T(i) is at z(i)
!         k(i) is midway between z(i-1) and z(i)
!         rhoc(i) is midway between z(i-1) and z(i)
!***********************************************************************
  implicit none

  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), dt, Tsurf, Tsurfp1, ti(nz), rhoc(nz)
  real(8), intent(IN) :: Fgeotherm
  real(8), intent(INOUT) :: T(nz)
  real(8), intent(OUT) :: Fsurf
  integer i
  real(8) alpha(nz), k(nz), gamma(nz), buf
  real(8) a(nz), b(nz), c(nz), r(nz)
  
  ! set some constants
  k(:) = ti(:)**2 / rhoc(:) ! thermal conductivity
  alpha(1) = dt * k(2) / rhoc(1) / (z(2)-z(1)) / z(2)
  gamma(1) = dt * k(1) / rhoc(1) / z(1) / z(2)
  do i=2,nz-1
     buf = 2.*dt / (z(i+1)-z(i-1)) / (rhoc(i)+rhoc(i+1))
     alpha(i) = k(i+1) * buf/ (z(i+1)-z(i))
     gamma(i) = k(i) * buf / (z(i)-z(i-1))
  enddo
  gamma(nz) = dt * k(nz) / (2.*rhoc(nz)) / (z(nz)-z(nz-1))**2
  
  ! elements of tridiagonal matrix
  a(:) = -gamma(:)   !  a(1) is not used
  b(:) = 1. + alpha(:) + gamma(:)
  c(:) = -alpha(:)   !  c(nz) is not used
  a(nz) = -2*gamma(nz)
  b(nz) = 1. + 2*gamma(nz)
  
  ! Set RHS         
  r(1)= alpha(1)*T(2) + (1.-alpha(1)-gamma(1))*T(1) + gamma(1)*(Tsurf+Tsurfp1)
  do concurrent (i=2:nz-1)
     r(i) = gamma(i)*T(i-1) + (1.-alpha(i)-gamma(i))*T(i) + alpha(i)*T(i+1)
  enddo
  r(nz) = 2.*gamma(nz)*T(nz-1) + (1.-2.*gamma(nz))*T(nz) + &
       &     2.*dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1))

  ! Solve for T at n+1
  call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion
  
  Fsurf = -k(1)*(T(1)-Tsurfp1)/z(1) ! heat flux into surface
  
end subroutine conductionT
