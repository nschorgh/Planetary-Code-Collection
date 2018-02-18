module conductionQ
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
!   This version of conductionQ pre-computes coefficients and
!         stores them privately
!
!   orginally written by Samar Khatiwala, 2001
!   extended to variable thermal properties
!         and irregular grid by Norbert Schorghofer
!***********************************************************************

  implicit none
  real*8, private :: k1, r_geo
  real*8, private, allocatable :: alpha(:), gamma(:), a(:), b(:), c(:)

  
contains
  subroutine conductionQ2_init(nz,z,dt,ti,rhoc,Fgeotherm)
!***********************************************************************
!   calculate coefficients alpha(:), gamma(:), a(:), b(:), c(:)
!   has no public output
!
!   dt = time step [s]
!   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
!   rhoc = rho*c  VECTOR where rho=density [kg m^-3] and 
!                              c=specific heat [J K^-1 kg^-1]
!   ti and rhoc are not allowed to vary in the layers immediately adjacent
!         to the surface or the bottom
!
!   Grid: surface is at z=0
!         T(1) is at z(1); ...; T(i) is at z(i)
!         k(i) is midway between z(i-1) and z(i); same for rhoc(i)
!
!   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]    
!***********************************************************************
    
    implicit none
    integer, intent(IN) :: nz
    real*8, intent(IN) :: z(nz),dt,ti(nz),rhoc(nz),Fgeotherm
    integer i
    real*8 k(nz), dz, beta, buf

    if (.not. allocated(alpha)) then
       allocate(alpha(nz), gamma(nz), a(nz), b(nz), c(nz))
    endif
    
    ! set some constants
    dz=2.*z(1)
    k(:) = ti(:)**2/rhoc(:) ! thermal conductivity
    k1=k(1)/dz   ! = k(1)/(2*z(1))

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

    a(:) = -gamma(:)   !  a(1) is not used
    b(2:nz-1) = 1. + alpha(2:nz-1) + gamma(2:nz-1) ! b(1) will be set later
    c(:) = -alpha(:)   !  c(nz) is not used
    b(nz) = 1. + gamma(nz)
    
    r_geo = dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1)) ! assumes rhoc(nz+1)=rhoc(nz)
  end subroutine conductionQ2_init


  
  subroutine conductionQ2(nz,Qn,Qnp1,T,emiss,Tsurf,Fsurf)
!***********************************************************************
!   Qn,Qnp1 = net solar insolation at time steps n and n+1 [Watts/m^2]
!   T = vertical temperature profile [K]  (in- and output)    
!   emiss = emissivity
!   Tsurf = surface temperature [K]  (in- and output)
!   Fsurf = heat flux into surface [W/m^2]  (output)    
!***********************************************************************

    implicit none
    real*8, parameter :: sigSB=5.6704d-8
    
    integer, intent(IN) :: nz
    real*8, intent(IN) :: Qn, Qnp1, emiss
    real*8, intent(INOUT) :: T(nz), Tsurf
    real*8, intent(OUT) :: Fsurf
    integer i
    real*8 Tr, r(nz)
    real*8 arad, brad, ann, annp1, bn

!   Emission
    Tr = Tsurf                !   'reference' temperature
    arad = -3.*emiss*sigSB*Tr**4
    brad = 2.*emiss*sigSB*Tr**3
    ann = (Qn-arad)/(k1+brad)
    annp1 = (Qnp1-arad)/(k1+brad)
    bn = (k1-brad)/(k1+brad)
    b(1) = 1. + alpha(1) + gamma(1) - gamma(1)*bn
    
!   Set RHS         
    r(1) = gamma(1)*(annp1+ann) + &
         &     (1.-alpha(1)-gamma(1)+gamma(1)*bn)*T(1) + alpha(1)*T(2)
    forall(i=2:nz-1)
       r(i) = gamma(i)*T(i-1) + (1.-alpha(i)-gamma(i))*T(i) + alpha(i)*T(i+1)
    end forall
    r(nz) = gamma(nz)*T(nz-1) + (1.-gamma(nz))*T(nz) + r_geo
    
!   Solve for T at n+1
    call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion
    
    Tsurf = 0.5*(annp1 + bn*T(1) + T(1)) ! (T0+T1)/2
    Fsurf = -2*k1*(T(1)-Tsurf)
    
  end subroutine conductionQ2

  ! note that alpha(1) and gamma(1) for conductionT are not the same as for conductionQ
end module conductionQ
