elemental subroutine flux_mars2(R,decl,latitude,HA, &
     &   surfaceSlope,azFac,emax,Q,Qscat,Qlw)
!*************************************************************************
!     This subroutine for solar insolation at Mars is based on the
!        function flux, but returns several irradiances in W/m^2
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
!     Q: direct incoming short-wavelength irradiance
!     Qscat: scattered (isotropic) short-wavelength irradiance from atmosphere
!     Qlw: (isotropic) long-wavelength irradiance from atmosphere
!*************************************************************************
  implicit none
  real(8), parameter :: pi=3.1415926535897931, So=1365., d2r=pi/180.
  real(8), intent(IN) :: R,decl,latitude,HA,surfaceSlope,azFac,emax
  real(8), intent(OUT) :: Q,Qscat,Qlw
  real(8) c1,s1,solarAttenuation,Q0
  real(8) sinbeta,sinbetaNoon,cosbeta,sintheta,azSun,buf
  real(8), parameter :: fracIR=0.04, fracDust=0.02
  
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
  if (cosbeta==0.) sintheta = cos(surfaceSlope)*sinbeta
  sintheta = cos(surfaceSlope)*sinbeta - &  
       &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)

  sintheta = max(sintheta,0.d0)  ! local horizon
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  if (sinbeta<sin(emax)) sintheta=0. ! distant horizon

! fluxes and contributions from atmosphere
  Q0 = So/R**2  ! solar constant at Mars
  solarAttenuation = (1.- fracIR - fracDust)**(1./max(sinbeta,0.04d0))
  Q = Q0*sintheta*solarAttenuation
  ! Qlw = fracIR*max(sinbetaNoon*Q0, sigSB*Tbase**4)  in polar region
  Qlw = sinbetaNoon*fracIR*Q0
  if (sinbeta>0.d0) then
     Qscat = 0.5*fracDust*Q0
  else
     Qscat = 0.
  endif
  ! absorbed flux = (1-albedo)*(Q+Qscat*viewfactor) + emiss*Qlw*viewfactor
end subroutine flux_mars2



