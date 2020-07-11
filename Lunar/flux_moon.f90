pure function flux_moon(R,decl,latitude,HA,albedo0)
!**********************************************************************
!   flux_moon: calculates absorbed solar flux without atmosphere;
!              albedo depends on incidence angle  
!     R: distance from sun (AU)
!     decl: planetocentric solar declination (radians)
!     latitude: (radians)
!     HA: hour angle (radians from noon, clockwise)
!     albedo0: albedo at zero incidence angle  
!**********************************************************************
  implicit none
  real(8) flux_moon
  real(8), parameter :: So=1365.  ! solar constant
  real(8), parameter :: pi=3.1415926535897931, d2r=pi/180.
  real(8), intent(IN) :: R, decl, latitude, HA, albedo0
  real(8) c1, s1, sinbeta, cosbeta, azSun, buf
  real(8) incid, albedo
  real(8), parameter :: a=0.06, b=0.25  ! Keihm, Icarus 60, 568 (1984)
  
  c1 = cos(latitude)*cos(decl)
  s1 = sin(latitude)*sin(decl)
  ! beta = 90 minus incidence angle for horizontal surface
  ! beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1
  
  cosbeta = sqrt(1-sinbeta**2)
  ! ha -> az
  buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta)
  ! buf can be NaN if cosbeta = 0
  if (buf>+1.) buf=1.0; if (buf<-1.) buf=-1.0; ! roundoff
  azSun = acos(buf)
  if (sin(HA)>=0) azSun=2*pi-azSun

  if (sinbeta>0.) then
     flux_moon = sinbeta*So/(R**2)
     incid = acos(sinbeta)
     !if (incid/=incid) error stop 'NaN'
     ! parametrization from Keihm, Icarus 60, 568 (1984) 
     albedo = albedo0 + a*(incid*4./pi)**3 + b*(incid*2./pi)**8
     flux_moon = (1-albedo)*flux_moon
  else
     flux_moon = 0.
  end if
  
end function flux_moon
