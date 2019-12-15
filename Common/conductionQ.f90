subroutine conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhoc,emiss,Tsurf, &
     &     Fgeotherm,Fsurf)
!***********************************************************************
!   conductionQ:  program to calculate the diffusion of temperature 
!                 into the ground and thermal emission at the surface 
!                 with variable thermal properties on irregular grid
!   Crank-Nicolson scheme, flux conservative
!                          uses Samar's radiation formula
!   Eqn: rhoc*T_t = (k*T_z)_z 
!   BC (z=0): Q(t) + kT_z = em*sig*T^4
!   BC (z=L): heat flux = Fgeotherm
!
!   nz = number of grid points
!   dt = time step
!   Qn,Qnp1 = net solar insolation at time steps n and n+1 [W/m^2]
!   T = vertical temperature profile [K]  (in- and output)  
!   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
!   rhoc = rho*c  VECTOR where rho=density [kg/m^3] and 
!                              c=specific heat [J K^-1 kg^-1]
!   ti and rhoc are not allowed to vary in the layers immediately 
!               adjacent to the surface or the bottom
!   emiss = emissivity
!   Tsurf = surface temperature [K]  (in- and output)
!   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
!   Fsurf = heat flux at surface [W/m^2]  (output)
!
!   Grid: surface is at z=0
!         z(0)=0, z(2)=3*z(1), i.e., the top layer has half the width
!         T(1) is at z(1); ...; T(i) is at z(i)
!         k(i), rhoc(i), ti(i) are midway between z(i-1) and z(i)
!     
!   originally written by Samar Khatiwala, 2001
!   extended to variable thermal properties
!         and irregular grid by Norbert Schorghofer
!   added predictor-corrector 9/2019  
!***********************************************************************

  implicit none
  real*8, parameter :: sigSB=5.6704d-8
  
  integer, intent(IN) :: nz
  real*8, intent(IN) :: z(nz), dt, Qn, Qnp1, ti(nz),rhoc(nz)
  real*8, intent(IN) :: emiss, Fgeotherm
  real*8, intent(INOUT) :: T(nz), Tsurf
  real*8, intent(OUT) :: Fsurf
  integer i, iter
  real*8 k(nz), k1, alpha(nz), gamma(nz), Tr
  real*8 a(nz), b(nz), c(nz), r(nz), Told(nz)
  real*8 arad, brad, ann, annp1, bn, buf, dz, beta
  
  ! set some constants
  k(:) = ti(:)**2/rhoc(:) ! thermal conductivity
  dz = 2.*z(1)
  beta = dt/rhoc(1)/(2.*dz**2)   ! assumes rhoc(0)=rhoc(1)
  alpha(1) = beta*k(2)
  gamma(1) = beta*k(1)
  do i=2,nz-1
     buf = dt/(z(i+1)-z(i-1))
     alpha(i) = 2.*k(i+1)*buf/(rhoc(i)+rhoc(i+1))/(z(i+1)-z(i))
     gamma(i) = 2.*k(i)*buf/(rhoc(i)+rhoc(i+1))/(z(i)-z(i-1))
  enddo
  buf = dt/(z(nz)-z(nz-1))**2
  gamma(nz) = k(nz)*buf/(2*rhoc(nz)) ! assumes rhoc(nz+1)=rhoc(nz)
  
  k1 = k(1)/dz
  
  ! elements of tridiagonal matrix
  a(:) = -gamma(:)   !  a(1) is not used
  b(:) = 1. + alpha(:) + gamma(:) !  b(1) has to be reset at every timestep
  c(:) = -alpha(:)   !  c(nz) is not used
  b(nz) = 1. + gamma(nz)

  Tr = Tsurf            ! 'reference' temperature
  iter = 0
  Told(:) = T(:)
30 continue  ! in rare case of iterative correction

  ! Emission
  arad = -3.*emiss*sigSB*Tr**4
  brad = 2.*emiss*sigSB*Tr**3
  ann = (Qn-arad)/(k1+brad)
  annp1 = (Qnp1-arad)/(k1+brad)
  bn = (k1-brad)/(k1+brad)
  b(1) = 1. + alpha(1) + gamma(1) - gamma(1)*bn
  
  ! Set RHS         
  r(1) = gamma(1)*(annp1+ann) + &
       &     (1.-alpha(1)-gamma(1)+gamma(1)*bn)*T(1) + alpha(1)*T(2)
  do concurrent (i=2:nz-1)
     r(i) = gamma(i)*T(i-1) + (1.-alpha(i)-gamma(i))*T(i)+ alpha(i)*T(i+1)
  enddo
  r(nz) = gamma(nz)*T(nz-1) + (1.-gamma(nz))*T(nz) &
       &     + dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1)) ! assumes rhoc(nz+1)=rhoc(nz)

  ! Solve for T at n+1
  call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion
  
  Tsurf = 0.5*(annp1 + bn*T(1) + T(1)) ! (T0+T1)/2

  ! iterative predictor-corrector
  if ((Tsurf > 1.2*Tr .or. Tsurf < 0.8*Tr) .and. iter<10) then  ! linearization error expected
     ! redo until Tr is within 20% of new surface temperature
     ! (under most circumstances, the 20% threshold is never exceeded)
     iter = iter+1
     Tr = sqrt(Tr*Tsurf)  ! linearize around an intermediate temperature
     T(:) = Told(:)
     goto 30
  endif
  !if (iter>=5) print *,'consider taking shorter time steps',iter
  !if (iter>=10) print *,'warning: too many iterations'
  
  Fsurf = - k(1)*(T(1)-Tsurf)/z(1) ! heat flux into surface
end subroutine conductionQ
