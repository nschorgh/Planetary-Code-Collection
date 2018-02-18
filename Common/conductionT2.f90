module conductionT
!***********************************************************************
!   conductionT:  program to calculate the diffusion of temperature 
!                 into the ground with prescribed surface temperature 
!                 and variable thermal properties on irregular grid
!   Crank-Nicolson scheme, flux conservative
!
!   Eqn: rhoc*T_t = (k*T_z)_z 
!   BC (z=0): T=T(t)
!   BC (z=L): heat flux = Fgeotherm
!
!   This version of conductionQ pre-computes coefficients and
!         stores them privately
!***********************************************************************

  implicit none
  real*8, private :: k1, r_geo
  real*8, private, allocatable :: alpha(:), gamma(:), a(:), b(:), c(:)


contains
  subroutine conductionT2_init(nz,z,dt,ti,rhoc,Fgeotherm)
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
!         k(i) is midway between z(i-1) and z(i)
!         rhoc(i) is midway between z(i-1) and z(i)
!
!   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]    
!***********************************************************************

    implicit none
    integer, intent(IN) :: nz
    real*8, intent(IN) :: z(nz), dt, ti(nz), rhoc(nz), Fgeotherm
    integer i
    real(8) k(nz), buf
    
    if (.not. allocated(alpha)) then
       allocate(alpha(nz), gamma(nz), a(nz), b(nz), c(nz))
    endif
       
    !   set some constants
    k(:) = ti(:)**2/rhoc(:) ! thermal conductivity
    alpha(1) = k(2)*dt/rhoc(1)/(z(2)-z(1))/z(2) 
    gamma(1) = k(1)*dt/rhoc(1)/z(1)/z(2) 
    do i=2,nz-1
       buf=dt/(z(i+1)-z(i-1))
       alpha(i) = k(i+1)*buf*2./(rhoc(i)+rhoc(i+1))/(z(i+1)-z(i))
       gamma(i) = k(i)*buf*2./(rhoc(i)+rhoc(i+1))/(z(i)-z(i-1))
    enddo
    buf=dt/(z(nz)-z(nz-1))**2
    gamma(nz) = k(nz)*buf/(rhoc(nz)+rhoc(nz)) ! assumes rhoc(nz+1)=rhoc(nz)

!   elements of tridiagonal matrix
    a(:) = -gamma(:)   !  a(1) is not used
    b(:) = 1. + alpha(:) + gamma(:)
    c(:) = -alpha(:)   !  c(nz) is not used
    b(nz) = 1. + gamma(nz)

    k1 =  k(1)/z(1)
    r_geo = dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1)) ! assumes rhoc(nz+1)=rhoc(nz)
  end subroutine conductionT2_init

  
  subroutine conductionT2(nz,Tsurf,Tsurfp1,T,Fsurf)
!***********************************************************************
!   Tsurf, Tsurfp1 = surface temperatures at times n and n+1
!   T = vertical temperature profile [K]
!   Fsurf = heat flux into surface [W/m^2] 
!***********************************************************************

    implicit none
    integer, intent(IN) :: nz
    real*8, intent(IN) :: Tsurf, Tsurfp1
    real*8, intent(INOUT) :: T(nz)
    real*8, intent(OUT) :: Fsurf
    integer i
    real*8 r(nz)

!   Set RHS         
    r(1)= alpha(1)*T(2) + (1.-alpha(1)-gamma(1))*T(1) + gamma(1)*(Tsurf+Tsurfp1)
    forall(i=2:nz-1)
       r(i) = gamma(i)*T(i-1) + (1.-alpha(i)-gamma(i))*T(i) + alpha(i)*T(i+1)
    end forall
    r(nz) = gamma(nz)*T(nz-1) + (1.-gamma(nz))*T(nz) + r_geo
      
!   Solve for T at n+1
    call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion
      
    Fsurf = -k1*(T(1)-Tsurfp1) ! heat flux into surface
  end subroutine conductionT2

end module conductionT
