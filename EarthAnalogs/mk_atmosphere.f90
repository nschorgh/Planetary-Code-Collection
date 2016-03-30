
real(8) function mk_atmosphere(Z,I0,D0)
  ! returns direct and indirect irradiance on Mauna Kea in W/m^2
  ! Nunez (1980), J. Biogeogr. 7, 173-186; paper full of typos
  implicit none
  real(8), intent(IN) :: Z  ! solar zenith angle (radians)
  real(8), intent(OUT) :: I0  ! clear-sky direct irradiance (unitless fraction) = transmittance
  real(8), intent(OUT) :: D0  ! clear-sky diffuse irradiance (unitless fraction)

  real(8) psi_wa  ! water vapor absorption
  real(8) psi_ws  ! water vapor scattering
  real(8) psi_rs  ! Rayleight scattering
  real(8) psi_da  ! dust absorption
  real(8) psi_ds  ! dust scattering
  real(8) w    ! precipitable water vapour (cm)
  real(8) m, mprime    ! optical air mass 
  real(8) p0   ! total pressure (Pa)

! Relative air mass
  ! m = 1/cos(Z)   ! simplest approximation
  ! Kasten (1966), Arch. Meteor. Geophys. Bioklimatol. B14, 206-223
  m = 1./(cos(Z) + 0.15*(93.885 - Z)**(-1.253)) 
  ! m = 1/(cos(Z) + 0.025*exp(-11*cos(Z)))
  if (m<0.) m=1e38

  p0 = 610   ! pressure on Mauna Kea 
  mprime = m*p0/1013.

! Water vapor
  w = 0.16   ! for Mauna Kea (cm)
  ! www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/mk-water-vapour-statistics

  psi_wa = 1 - 0.077*(w*m)**0.30  ! McDonald (1960), J. Meteor. 17, 319-328
  psi_ws = 1 - 0.025*w*m 
  if (psi_wa<0.) psi_wa=0.
  if (psi_ws<0.) psi_ws=0.

! Rayleigh scattering
  !psi_rs = 0.972 - 0.08262*m + 0.0933*m**2 - 0.00095*m**3 + 0.0000437*m**4 ! wrong
  psi_rs = exp(-0.08*mprime)  ! 8% at sea level according to Fig 3-3 in Bird & Hulstrom (1981)

! Aerosols 
  !psi_da = 0.94**(m/2)  ! Nunez (1980)
  !psi_ds = 0.94**(m/2)  ! Nunez (1980)
  ! aerosol optical depth on MK = 0.0084*(lambda/1um)**(-1.26), Buton et al. (2013), A&A 549, A8
  psi_ds = exp(-m*0.0084*0.5**(-1.26)) 
  psi_da = psi_ds   ! assumes single scattering albedo of 0.5

! Direct sunlight
  I0 = psi_wa*psi_da*psi_ws*psi_rs*psi_ds  ! multiply with solar constant to get flux

! Diffusive sunlight
  D0 = I0*cos(Z)*psi_wa*psi_da*(1-psi_ws*psi_rs*psi_ds)/2.
  ! D0 = D0*(sky view)  ! on inclined surface
  if (D0<=0.) D0=0.   ! because of -0

  mk_atmosphere = I0*cos(Z) + D0
  if (mk_atmosphere<=0.) mk_atmosphere=0.   ! because of -0

  !print *,Z,psi_wa,psi_ws,psi_rs,psi_ds
  !print *,Z,m,I0,D0,mk_atmosphere

  ! Sensible heat flux (depends on U)

  ! Long-wave downward radiation
  ! effective atmospheric emissivity 
  !wmbar = w/0.98  ! cm -> mbar
  !eatm = 0.62*(wmbar)**0.08     ! Staley & Jurica (1972). J. Appl. Meteor. 11, 349-356
  !eatm = 0.605 + 0.048*(wmbar)**0.5  ! Brunt's equation
  !Tatm = 273+4    ! temperature of the atmosphere at 2m height
  !I_LW = eatm*5.67e-8*Tatm**4   ! multiply this with surface emissivity
  !I_LW = I_LW*(sky view)

  ! annual mean cloud frequency 
  ! 0.26 ! according to http://evapotranspiration.geography.hawaii.edu/
end function mk_atmosphere


