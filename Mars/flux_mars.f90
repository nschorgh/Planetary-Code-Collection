pure function flux_mars77(R,decl,latitude,HA,albedo,fracir,fracdust)
  !***********************************************************************
  !   flux_mars77: calculates insolation at Mars
  !     flat surface only; also works in polar regions
  !
  !     R: distance from sun [AU]
  !     decl: planetocentric solar declination [radians]
  !     latitude: [radians]
  !     HA: hour angle [radians from noon, clockwise]
  !     fracir: fraction of absorption
  !     fracdust: fraction of scattering
  !***********************************************************************
  implicit none
  real*8 flux_mars77
  real*8, intent(IN) :: R,decl,latitude,HA,albedo,fracIR,fracDust
  real*8, parameter :: pi=3.1415926535897931, So=1365., d2r=pi/180.
  real*8, parameter :: sigSB=5.6704d-8
  real*8 c1,s1,Qo,solarAttenuation,Q
  real*8 sinbeta,sinbetaNoon,cosbeta,azSun,buf,Fout

  c1=cos(latitude)*cos(decl)
  s1=sin(latitude)*sin(decl)
  ! beta = 90 minus incidence angle for horizontal surface
  ! beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1
  sinbetaNoon = c1 + s1
  sinbetaNoon = max(sinbetaNoon,0.d0)  ! horizon
  
  cosbeta = sqrt(1-sinbeta**2)
  ! ha -> az (option 1)
  ! azSun=asin(-cos(decl)*sin(ha)/cosbeta)
  ! ha -> az (option 2)
  buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta)
  if (buf>+1.) buf=1.d0   ! roundoff
  if (buf<-1.) buf=-1.d0  ! roundoff
  azSun = acos(buf)
  if (sin(HA)>=0) azSun=2*pi-azSun
  ! ha -> az (option 3)  without beta
  ! azSun=sin(latitude)*cos(decl)*cos(ha)-cos(latitude)*sin(decl)
  ! azSun=atan(sin(ha)*cos(decl)/azSun)
  ! print *,asin(sinbeta)/d2r,azSun/d2r
  
  Qo = So/(R**2)
  ! atmospheric contributions are based on Kieffer et al. (1977), JGR
  if (sinbeta>0.d0) then
     solarAttenuation = (1.-albedo)* &
          &        (1.- fracIR - fracDust)**(1./max(sinbeta,0.04)) 
     Q = Qo*(sinbeta*solarAttenuation + 0.5*fracDust)
  else
     Q = 0.
  endif
  ! net flux: short-wavelength insolation + IR
  Fout = 1.*sigSB*150**4  ! matters only in polar region
  Q = Q + fracIR*max(Qo*sinbetaNoon,Fout)
  flux_mars77 = Q
end function flux_mars77



pure subroutine flux_mars2(R,decl,latitude,HA,fracIR,fracDust, &
     &   surfaceSlope,azFac,emax,Q,Qscat,Qlw)
!***********************************************************************
!     This subroutine for solar insolation at Mars is based on the
!       function flux, but returns several irradiances [W/m^2]
!
!     R: distance from sun [AU]
!     decl: planetocentric solar declination [radians]
!     latitude: [radians]
!     HA: hour angle [radians from noon, clockwise]
!     fracir: fraction of absorption
!     fracdust: fraction of scattering
!     surfaceSlope: >0, [radians]
!     azFac: azimuth of gradient (radians east of north)
!     emax: maximum horizon elevation in direction of azimuth [radians]
!     Q: direct incoming short-wavelength irradiance [W/m^2]
!     Qscat: diffuse short-wavelength irradiance from atmosphere [W/m^2]
!     Qlw: diffuse long-wavelength irradiance from atmosphere [W/m^2]
!***********************************************************************
  implicit none
  real(8), parameter :: pi=3.1415926535897931, So=1365., d2r=pi/180.
  real(8), intent(IN) :: R,decl,latitude,HA,surfaceSlope,azFac,emax
  real(8), intent(IN) :: fracIR,fracDust
  real(8), intent(OUT) :: Q,Qscat,Qlw
  real(8) c1,s1,solarAttenuation,Qo
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
  sintheta = cos(surfaceSlope)*sinbeta - &  
       &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)

  if (sintheta<0.) sintheta=0.  ! local horizon
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  if (sinbeta<sin(emax)) sintheta=0. ! distant horizon

! fluxes and contributions from atmosphere
  Qo = So/R**2  ! solar constant at Mars
  solarAttenuation = (1.- fracIR - fracDust)**(1./max(sinbeta,0.04d0))
  Q = Qo*sintheta*solarAttenuation
  ! Fout = 1.*sigSB*150**4
  ! Qlw = fracIR*max(Qo*sinbetaNoon,Fout)  in polar region
  Qlw = fracIR*Qo*sinbetaNoon
  if (sinbeta>0.d0) then
     Qscat = 0.5*fracDust*Qo
  else
     Qscat = 0.
  endif
  ! absorbed flux = (1-albedo)*(Q+Qscat*viewfactor) + emiss*Qlw*viewfactor
end subroutine flux_mars2



