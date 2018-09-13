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



elemental subroutine flux_mars2(R,decl,latitude,HA, &
     &   surfaceSlope,azFac,emax,Q,Qscat,Qlw)
!***********************************************************************
!     This subroutine for solar insolation at Mars is based on the
!        function flux_mars, but returns several quantities  
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
!     Q: direct incoming short-wavelength flux
!     Qscat: scattered (isotropic) short-wavelength flux from atmosphere
!     Qlw: (isotropic) long-wavelength flux from atmosphere
!***********************************************************************
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



