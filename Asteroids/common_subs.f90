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


pure function interp1(x0,x,y0,y,xi,nz)
  ! linear interpolation
  ! x0<x(1)<x(2)<...<x(nz)
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: x(nz),y(nz),xi,x0,y0
  integer i
  real(8) interp1, yi
  
  yi = -9.
  do i=1,nz-1
     if (x(i)>x(i+1)) error stop 'incorrect direction'
     if (x(i)<xi) then
        yi = (y(i)*(x(i+1)-xi)+y(i+1)*(xi-x(i)))/(x(i+1)-x(i))
     endif
  enddo
  if (xi<x(1)) yi = (y0*x(1)-y(1)*x0)/(x(1)-x0)
  
  interp1 = yi
end function interp1


pure function heatcapacity(T)
  ! specific heat capacity of silicates
  implicit none
  real(8), intent(IN) :: T  ! [K]
  real(8) c, heatcapacity   ! [J/(kg K)]
  
  ! Ledlow et al., ApJ 348, 640 (1992), <350K
  !c = 0.1812 + 0.1191*(T/300.-1) + 0.0176*(T/300.-1)**2 + &
  !     0.2721*(T/300.-1)**3 + 0.1869*(T/300.-1)**4
  !heatcapacity = c*1000*4.184  ! cal/(g K) -> J/(kg K)

  ! Winter & Saari, ApJ 156, 1135 (1969), 20K<T<500K
  c = -0.034*T**0.5 + 0.008*T - 0.0002*T**1.5
  heatcapacity = c*1000   ! J/(g K) -> J/(kg K)

  ! Hayne et al., JGR 122, 2371 (2017)
  !c = 8.9093E-09*T**4 -1.2340E-05*T**3 +2.36160E-03*T**2 + 2.7431*T -3.6125
  !heatcapacity = c
end function heatcapacity


pure function vapordiffusivity(diam,porosity,T)
  ! vapor diffusivity of porous material [m^2/s]
  ! diam = rms grain diameter
  implicit none
  real(8) vapordiffusivity
  real(8), parameter :: pi = 3.1415926535897932
  real(8), parameter :: Ru = 8314.46
  real(8), intent(IN) :: diam, porosity, T
  real(8) vbar, r
  real(8), parameter :: tau = 2.  ! tortuosity

  r = diam/2.
  vbar = sqrt(8*Ru*T/(pi*18))
  ! for 0<=porosity<=0.5
  vapordiffusivity = pi/(8+pi)*porosity/(1-porosity)*vbar*r/tau
end function vapordiffusivity


pure function faintsun(t)
  ! solar constant of the past
  implicit none
  real(8) faintsun
  real(8), intent(IN) :: t   ! time before present [years]
  ! Gough, D. O. (1981), Sol. Phys., 74, 21–34
  faintsun = 1./(1+0.4*abs(t)/4.57e9)
end function faintsun


pure function gettype(zdepth,nz,z)
  implicit none
  integer gettype
  integer, intent(IN) :: nz
  real(8), intent(IN) :: zdepth, z(nz)
  integer j
  gettype = -9 
  do j=1,nz
     if (z(j)>zdepth) then
        gettype = j  
        exit
     endif
  enddo
end function gettype

