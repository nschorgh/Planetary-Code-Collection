!**********************************************************************
! Regolith thermal properties for the Moon
! Options include those used in the "Diviner standard thermal model"
! described in Hayne et al., JGR 122, 2371 (2017)
!**********************************************************************


function solidconductivity(z)
  ! thermal conduction through solid
  implicit none
  real(8) solidconductivity
  real(8), intent(IN) :: z
  real(8) kc, k_s, k_d
  real(8), external :: twolayers
  
  ! Hayne et al. (2017) JGR
  !rho_s = 1100; rho_d = 1800
  k_s = 7.4e-4; k_d = 3.4e-3

  !kc = k_d - (k_d - k_s) * exp(-z/H)
  kc = twolayers(k_s, k_d, z)
  solidconductivity = kc 

  ! rho = (1-porosity)*rhosolid,  porosity = 1-rho/rhosolid
  ! 1-1100/2500 = 0.560 ! s
  ! 1-1800/2500 = 0.280 ! d
  ! porosity(z) = 0.28 + (0.56-0.28)*exp(-z/H)
end function solidconductivity



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
  !c = -0.034*T**0.5 + 0.008*T - 0.0002*T**1.5
  !heatcapacity = c*1000   ! J/(g K) -> J/(kg K)

  ! Hayne et al., JGR 122, 2371 (2017)
  c = 8.9093E-09*T**4 -1.2340E-05*T**3 +2.36160E-03*T**2 + 2.7431*T -3.6125
  heatcapacity = c
end function heatcapacity



pure function twolayers(s,d,z)
  ! s ... surface value (z=0)
  ! d ... value at great depth (z=Inf)
  ! z ... depth
  implicit none
  real(8) twolayers
  real(8), intent(IN) :: s, d, z
  real(8) H
  H = 0.07 ! Hayne et al. (2017) average
  !H = 0.06 ! Hayne et al. (2017) Fig. A2
  twolayers = d - (d-s) * exp(-z/H)
end function twolayers



function twolayerdensity(z)
  implicit none
  real(8) twolayerdensity
  real(8), intent(IN) :: z
  real(8), external :: twolayers
  real(8) rho_s, rho_d
  ! Hayne et al. (2017)
  rho_s = 1100; rho_d = 1800
  twolayerdensity = twolayers(rho_s,rho_d,z)
end function twolayerdensity
