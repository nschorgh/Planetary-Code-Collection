pure function flux_mars77(R,decl,latitude,HA,albedo,fracir,fracscat)
!***********************************************************************
! flux_mars77: calculates insolation (incoming solar radiation) on Mars
!     flat surface only; also works in polar regions
!
!     R: distance from sun [AU]
!     decl: planetocentric solar declination [radians]
!     latitude: [radians]
!     HA: hour angle [radians from noon, clockwise]
!     fracir: fraction of absorption
!     fracscat: fraction of scattering
!***********************************************************************
  implicit none
  real(8) flux_mars77
  real(8), intent(IN) :: R,decl,latitude,HA,albedo,fracIR,fracScat
  real(8), parameter :: So=1365. ! [W/m^2]
  real(8), parameter :: sigSB=5.6704d-8
  real(8) c1,s1,Qo,solarAttenuation,Q
  real(8) sinbeta,sinbetaNoon,Qdir,Qscat,Qlw,Fout

  c1 = cos(latitude)*cos(decl)
  s1 = sin(latitude)*sin(decl)
  ! beta = 90 minus incidence angle for horizontal surface
  ! beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1
  sinbetaNoon = c1 + s1
  sinbetaNoon = max(sinbetaNoon,0.d0)  ! horizon

  Qo = So/(R**2)  ! solar constant at Mars
  
  ! short-wavelength irradiance
  if (sinbeta>0.d0) then
     solarAttenuation = (1.- fracIR - fracScat)**(1./max(sinbeta,0.04))
     Qdir = Qo*sinbeta*solarAttenuation
     Qscat = 0.5*Qo*fracScat
     Q = (1.-albedo)*(Qdir + Qscat)
  else
     Q = 0.
  endif

  ! atmospheric IR contribution based on Kieffer et al. (1977), JGR 82, 4249
  Fout = 1.*sigSB*150**4  ! matters only in polar region
  Qlw = fracIR*max(Qo*sinbetaNoon,Fout)

  ! net flux = short-wavelength + long-wavelength irradiance
  flux_mars77 = Q + Qlw
end function flux_mars77



pure subroutine flux_mars2(R,decl,latitude,HA,fracIR,fracScat, &
     &   SlopeAngle,azFac,emax,Qdir,Qscat,Qlw)
!*****************************************************************************
! flux_mars2: Insolation (incoming solar radiation) at Mars on sloped surface;
!             returns several irradiances
!
! INPUTS:
!     R: distance from sun [AU]
!     decl: planetocentric solar declination [radians]
!     latitude: [radians]
!     HA: hour angle (radians from noon, clockwise)
!     fracIR: fraction of absorption
!     fracScat: fraction of scattering
!     SlopeAngle: >0, [radians]
!     azFac: azimuth of topographic gradient (radians east of north)
!            azFac=0 is south-facing  
!     emax: maximum horizon elevation in direction of azimuth [radians]
! OUTPUTS:
!     Qdir: direct incoming short-wavelength irradiance [W/m^2]
!     Qscat: diffuse short-wavelength irradiance from atmosphere [W/m^2]
!     Qlw: diffuse long-wavelength irradiance from atmosphere [W/m^2]
!*****************************************************************************
  implicit none
  real(8), parameter :: pi=3.1415926535897932, So=1365.
  real(8), parameter :: sigSB=5.6704d-8
  real(8), intent(IN) :: R,decl,latitude,HA,SlopeAngle,azFac,emax
  real(8), intent(IN) :: fracIR,fracScat
  real(8), intent(OUT) :: Qdir,Qscat,Qlw
  real(8) c1,s1,solarAttenuation,Qo
  real(8) sinbeta,sinbetaNoon,cosbeta,sintheta,azSun,buf,Fout
  
  c1 = cos(latitude)*cos(decl)
  s1 = sin(latitude)*sin(decl)
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
  if (sin(HA)>=0) azSun = 2*pi-azSun
  
! theta = 90 minus incidence angle for sloped surface
  sintheta = cos(SlopeAngle)*sinbeta - &  
       &     sin(SlopeAngle)*cosbeta*cos(azSun-azFac)
  if (cosbeta==0) sintheta=cos(SlopeAngle)*sinbeta  ! zenith, where azSun=NaN
  
  if (sintheta<0.) sintheta=0.  ! local horizon (self-shadowing)
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  if (sinbeta<sin(emax)) sintheta=0. ! distant horizon

! fluxes and contributions from atmosphere
  Qo = So/R**2  ! solar constant at Mars
  solarAttenuation = (1.- fracIR - fracScat)**(1./max(sinbeta,0.04d0))
  Qdir = Qo*sintheta*solarAttenuation
  Fout = 1.*sigSB*150**4  ! matters only in polar region
  Qlw = fracIR*max(Qo*sinbetaNoon,Fout) 
  if (sinbeta>0.d0) then
     Qscat = 0.5*fracScat*Qo
  else
     Qscat = 0.
  endif
  
! For a horizontal and unobstructed surface
!   absorbed flux = (1-albedo)*(Qdir+Qscat) + emiss*Qlw
!
! For a tilted surface
!   absorbed flux = (1-albedo)*(Qdir+Qscat*viewfactor) + emiss*Qlw*viewfactor
!   then add irradiance from land in field of view  
!   in the case of a planar slope, viewfactor = cos(SlopeAngle/2.)**2
end subroutine flux_mars2
