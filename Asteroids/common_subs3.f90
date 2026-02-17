!***********************************************************************
! a collection of commonly-used functions
!***********************************************************************

pure function flux2T(Q,albedo,emiss)
  ! convert incoming flux Q to equilibrium temperature
  implicit none
  real(8) flux2T
  real(8), intent(IN) :: Q, emiss, albedo
  real(8), parameter :: sigSB = 5.6704d-8

  flux2T = ((1-albedo)*Q/sigSB/emiss)**0.25
end function flux2T


pure function a2Torb(semia)
  ! returns orbital period in Earth days
  implicit none
  real(8) a2Torb
  real(8), parameter :: pi=3.1415926535897932
  real(8), intent(IN) :: semia  ! semimajor axis [AU]
  real(8) T  ! orbital period [sec]
  
  T = sqrt(4*pi**2/(6.674e-11*1.989e30)*(semia*149.598e9)**3)
  a2Torb = T/86400.
end function a2Torb


pure integer function sols_per_year(semia,solarDay)
  ! number of solar days per orbit, rounded
  implicit none
  real(8), parameter :: pi=3.1415926535897932
  real(8), intent(IN) :: semia, solarDay
  real(8) T  ! orbital period [sec]

  ! the numbers match those in subroutine generalorbit
  T = sqrt(4*pi**2/(6.674e-11*1.989e30)*(semia*149.598e9)**3)
  sols_per_year = nint(T/solarDay)
end function sols_per_year


pure function faintsun(t)
  ! solar constant of the past
  implicit none
  real(8) faintsun
  real(8), intent(IN) :: t   ! time before present [years]
  ! Gough, D. O. (1981), Sol. Phys., 74, 21-34
  faintsun = 1./(1+0.4*abs(t)/4.57e9)
end function faintsun

