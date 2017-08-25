elemental real(8) function flux_wshad(R,sinbeta,azSun,surfaceSlope,azFac,smax)
!***********************************************************************
!   flux:  program to calculate incoming solar flux without atmosphere
!     R: distance from sun (AU)
!     sinbeta: sin(altitude) 
!     azSun: azimuth of Sun (radians east of north)
!     surfaceSlope: >=0, (radians) 
!     azFac: azimuth of topographic gradient (radians east of north)
!     smax: elevation of horizon in direction of azimuth (radians)
!***********************************************************************
  implicit none
  real(8), parameter :: So=1365.  ! solar constant
  real(8), intent(IN) :: R,azSun,sinbeta,surfaceSlope,azFac,smax
  real(8) cosbeta,sintheta
  
  cosbeta = sqrt(1.-sinbeta**2)

!-incidence angle
  ! theta = 90 minus incidence angle for sloped surface
  sintheta = cos(surfaceSlope)*sinbeta - &
       &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)
  if (cosbeta==0.) sintheta = cos(surfaceSlope)*sinbeta ! does not use azimuths

!-shadowing
  sintheta = max(sintheta,0.d0)  ! self-shadowing
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  if (sinbeta<sin(smax)) sintheta=0.  ! shadowing

!-intensity
  flux_wshad = sintheta*So/(R**2)
end function flux_wshad

