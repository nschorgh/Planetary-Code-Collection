elemental function flux_mars(R,decl,latitude,HA,albedo, &
     &   fracir,fracdust,surfaceSlope,azFac,emax,viewfactor)
!***********************************************************************
!   flux:  program to calculate solar insolation
!     identical to function flux, but with distant horizon
!
!     R: distance from sun (AU)
!     decl: planetocentric solar declination (radians)
!     latitude: (radians)
!     HA: hour angle (radians from noon, clockwise)
!     fracir: fraction of absorption
!     fracdust: fraction of scattering
!     surfaceSlope: >0, (radians) 
!     azFac: azimuth of gradient (radians east of north)
!     emax: maximum horizon elevation in direction of azimuth (radians)
!     viewfactor: fraction of unobstructed sky (0...1)  
!***********************************************************************
  implicit none
  real(8) flux_mars
  real(8), parameter :: pi=3.1415926535897931, So=1365., d2r=pi/180.
  real(8), intent(IN) :: R,decl,latitude,HA,albedo,fracIR,fracDust
  real(8), intent(IN) :: surfaceSlope,azFac,emax,viewfactor
  real(8) c1,s1,solarAttenuation,Q
  real(8) sinbeta,sinbetaNoon,cosbeta,sintheta,azSun,buf

  c1=cos(latitude)*cos(decl)
  s1=sin(latitude)*sin(decl)
! beta = 90 minus incidence angle for horizontal surface
! beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1
  sinbetaNoon = c1 + s1
  sinbetaNoon = max(sinbetaNoon,0.d0)  ! horizon
      
  cosbeta = sqrt(1-sinbeta**2)
  buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta)
  if (buf>+1.) buf=+1.d0  ! roundoff
  if (buf<-1.) buf=-1.d0  ! roundoff
  azSun = acos(buf)
  if (sin(HA)>=0) azSun=2*pi-azSun
  
! theta = 90 minus incidence angle for sloped surface
  if (cosbeta==0.) sintheta = cos(surfaceSlope)*sinbeta ! new line
  sintheta = cos(surfaceSlope)*sinbeta - &   ! note the sign 
       &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)

  sintheta = max(sintheta,0.d0)  ! horizon
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  if (sinbeta<sin(emax)) sintheta=0. 

! net flux
  Q = (1.-albedo)*sintheta
! contributions from atmosphere
  solarAttenuation = (1.- fracIR - fracDust)**(1./max(sinbeta,0.04d0))
  Q = Q*solarAttenuation + viewfactor*sinbetaNoon*fracIR 
  if (sinbeta>0.d0) then
     Q = Q + (1.-albedo)*0.5*viewfactor*fracDust
  endif
  flux_mars=Q*So/R**2
  
end function flux_mars

