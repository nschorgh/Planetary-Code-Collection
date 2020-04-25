pure function flux_noatm(R,decl,latitude,HA,surfaceSlope,azFac)
!**********************************************************************
!   flux_noatm: calculates incoming solar flux without atmosphere
!     R: distance from sun (AU)
!     decl: planetocentric solar declination (radians)
!     latitude: (radians)
!     HA: hour angle (radians from noon, clockwise)
!     surfaceSlope: >0, (radians) 
!     azFac: azimuth of topographic gradient (radians east of north)
!**********************************************************************
  implicit none
  real(8) flux_noatm
  real(8), parameter :: So=1365.  ! solar constant [W/m^2]
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), intent(IN) :: R,decl,latitude,HA,surfaceSlope,azFac
  real(8) c1,s1,sinbeta,cosbeta,sintheta,azSun,buf
  
  c1=cos(latitude)*cos(decl)
  s1=sin(latitude)*sin(decl)
  ! beta = 90 minus incidence angle for horizontal surface
  ! beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1
  
  cosbeta = sqrt(1-sinbeta**2)
  ! ha -> az (option 1)
  !azSun=asin(-cos(decl)*sin(ha)/cosbeta)
  ! ha -> az (option 2)
  buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta)
  ! buf can be NaN if cosbeta = 0
  if (buf>+1.) buf=1.d0; if (buf<-1.) buf=-1.d0; ! roundoff
  azSun = acos(buf)
  if (sin(HA)>=0) azSun=2*pi-azSun
  ! ha -> az (option 3)  without beta
  !azSun=sin(latitude)*cos(decl)*cos(ha)-cos(latitude)*sin(decl)
  !azSun=atan(sin(ha)*cos(decl)/azSun)

  ! theta = 90 minus incidence angle for sloped surface
  sintheta = cos(surfaceSlope)*sinbeta + &
       &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)
  if (cosbeta==0.) sintheta = cos(surfaceSlope)*sinbeta
  sintheta = max(sintheta,0.d0)  ! horizon
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  
  flux_noatm = sintheta*So/(R**2)
  
  !write(6,'(99(1x,f6.2))') decl/d2r,HA/d2r,flux_noatm, &
  !     &     asin(sintheta)/d2r,asin(sinbeta)/d2r,azSun/d2r,buf
end function flux_noatm



pure function flux_wshad(R,sinbeta,azSun,surfaceSlope,azFac,emax)
!**********************************************************************
!   flux_wshad: calculates incoming solar flux without atmosphere
!     R: distance from sun (AU)
!     sinbeta: sin(altitude) 
!     azSun: azimuth of Sun (radians east of north)
!     surfaceSlope: >=0, (radians) 
!     azFac: azimuth of topographic gradient (radians east of north)
!     emax: elevation of horizon in direction of azimuth (radians)
!**********************************************************************
  implicit none
  real(8) flux_wshad
  real(8), parameter :: So=1365.  ! solar constant [W/m^2]
  real(8), intent(IN) :: R,azSun,sinbeta,surfaceSlope,azFac,emax
  real(8) cosbeta,sintheta
  
  cosbeta = sqrt(1.-sinbeta**2)

!-incidence angle
  ! theta = 90 minus incidence angle for sloped surface
  sintheta = cos(surfaceSlope)*sinbeta - &
       &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)
  if (cosbeta==0.) sintheta = cos(surfaceSlope)*sinbeta ! does not use azimuths

!-shadowing
  if (sintheta<0.) sintheta=0.  ! self-shadowing
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  if (sinbeta<sin(emax)) sintheta=0.  ! shadowing from distant horizon

!-intensity
  flux_wshad = sintheta*So/(R**2)
end function flux_wshad

