subroutine setgrid(nz,z,zmax,zfac)
  ! construct regularly or geometrically spaced 1D grid
  ! z(n)-z(1) = 2*z(1)*(zfac**(n-1)-1)/(zfac-1)
  ! choice of z(1) and z(2) is compatible with conductionQ
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: zfac, zmax
  real(8), intent(OUT) :: z(nz)
  integer i
  real(8) dz
  
  dz = zmax/nz
  do i=1,nz
     z(i) = (i-0.5)*dz
  enddo
  
  if (zfac>1.) then
     dz = zmax/(3.+2.*zfac*(zfac**(nz-2)-1.)/(zfac-1.))
     z(1) = dz
     z(2) = 3*z(1)
     do i=3,nz
        z(i) = (1.+zfac)*z(i-1) - zfac*z(i-2)
        ! z(i) = z(i-1) + zfac*(z(i-1)-z(i-2)) ! equivalent
     enddo
  endif
  
end subroutine setgrid



real(8) function smartzfac(nz_max,zmax,nz_delta,delta)
  ! output can be used as input to setgrid
  ! produces zfac with desired number of points within depth delta
  implicit none
  integer, intent(IN) :: nz_max, nz_delta
  real(8), intent(IN) :: zmax, delta
  integer j
  real(8) f, g, gprime
  
  if (nz_max<nz_delta .or. nz_delta<=1) then
     stop 'inappropriate input to smartzfac'
  endif
  f = 1.05
  do j=1,7  ! Newton iteration
     !print *,j,f
     g = (f-3+2*f**(nz_max-1))/(f-3+2*f**(nz_delta-1))-zmax/delta
     gprime=(1+2*(nz_max-1)*f**(nz_max-2))/(f-3+2*f**(nz_delta-1)) - &
          &        (f-3+2*f**(nz_max-1))*(1+2*(nz_delta-1)*f**(nz_delta-2))/ & 
          &        (f-3+2*f**(nz_delta-1))**2
     f = f-g/gprime
  enddo
  smartzfac = f
  if (smartzfac<1. .or. smartzfac>2.) then
     print *,'zfac=',smartzfac
     stop 'unwanted result in smartzfac'
  endif
end function smartzfac


      
!-grid-dependent utility functions

function colint(y,z,nz,i1,i2)
  ! column integrates y
  implicit none
  integer, intent(IN) :: nz, i1, i2
  real(8), intent(IN) :: y(nz), z(nz)
  real(8) colint
  integer i
  real(8) dz(nz)
  
  dz(1) = (z(2)-0.)/2
  do i=2,nz-1
     dz(i) = (z(i+1)-z(i-1))/2.
  enddo
  dz(nz) = z(nz)-z(nz-1)
  colint = sum(y(i1:i2)*dz(i1:i2))
end function colint


 
subroutine dzvector(nz,z,dz) 
  ! matches function colint
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz)
  real(8), intent(OUT) :: dz(nz)
  integer i
  
  dz(1) = (z(2)-0.)/2
  do i=2,nz-1
     dz(i) = (z(i+1)-z(i-1))/2.
  enddo
  dz(nz) = z(nz)-z(nz-1)
end subroutine dzvector



subroutine heatflux_from_temperature(nz,z,T,k,H)
  ! calculates heat flux from temperature profile
  ! like k, the heat flux H is defined mid-point
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), T(nz), k(nz)
  real(8), intent(OUT) :: H(nz)
  integer j
  
  ! H(1) = -k(1)*(T(1)-Tsurf)/z(1)
  H(1) = 0. ! to avoid ill-defined value
  do j=2,nz
     H(j) = -k(j)*(T(j)-T(j-1))/(z(j)-z(j-1))
  enddo
end subroutine heatflux_from_temperature
