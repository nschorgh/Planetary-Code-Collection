elemental function evap_vacuum_func(T)
  ! sublimation rate of H2O ice into vacuum
  implicit none
  real(8) evap_vacuum_func  ! [kg/m^2/s]
  real(8), intent(IN) :: T  ! [Kelvin]
  real(8), parameter :: mu=18.015, Ru=8314.46, pi=3.141592653589793
  real(8) psv
  
  ! eq. (7) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
  !psv = exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)
  !evap_vacuum_func = psv*sqrt(mu/(2*pi*R*T))

  ! first coefficient: add log(sqrt(mu/(2*pi*R))) 
  ! third coefficient: subtract 0.5
  evap_vacuum_func = exp(5.564214 -5723.265/T +3.03068*log(T) -0.00728332*T)

  if (T>273.16) then ! liquid, Hyland & Wexler (1983), Murphy & Koop (2005)
     psv = exp( -5800.2206/T +1.3914993 -0.48640239e-1*T +0.41764768e-4*T**2 &
          & -0.14452093e-7*T**3 + 6.5459673*log(T) )
     evap_vacuum_func = psv*sqrt(mu/(2*pi*Ru*T))
 endif
 
end function evap_vacuum_func



elemental function inv_evap_vacuum_func(PET)
  ! inverse of evap_vacuum_func for ice
  implicit none
  real(8) inv_evap_vacuum_func  ! [Kelvin]
  real(8), intent(IN) :: PET  ! [kg/m^2/s]
  real(8) T, dEdT, E
  real(8), parameter :: a1=5.564214, a2=-5723.265, a3=3.03068, a4=-0.00728332
  integer n
  T = (0.81004*log(PET) + 6049.7) / (21.762 - log(PET) ) ! approximate fit
  do n=1,4  ! refine with Newton iterations
     E = exp(a1 +a2/T +a3*log(T) +a4*T)
     dEdT = -a2/T**2 + a3/T + a4 ! times E
     T = T - (E-PET)/(E*dEdT)
  enddo
  inv_evap_vacuum_func = T
end function inv_evap_vacuum_func



real(8) function sublr_amorph(T)
  ! sublimation rate of amorphous or crystalline H2O ice
  ! in units of #molecules/m^2/s
  implicit none
  real(8), intent(IN) :: T  ! [Kelvin]
  real(8), parameter :: kB=1.38065e-23, m=18.015*1.66e-27
  real(8) E

  ! sublimation rate of amorphous H2O ice
  ! according to Sack and Baragiola, Phys. Rev. B 48, 9973 (1993)  
  E = 1.82e21*1e4*T**3.5*exp(-0.45*1.6022e-19/(kB*T)) ! [#/m^2/s]

  ! sublimation rate of crystalline H2O ice
  ! eq. (7) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
  ! psv = exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)
  ! E = psv/sqrt(2*pi*kB*T*m)
  ! first coefficient: subtract log(2*pi*kB*m)/2
  ! third coefficient: subtract 0.5
  !E = exp(64.3358 -5723.265/T +3.03068*log(T) -0.00728332*T)

  sublr_amorph = E
  !sublr_amorph = E*m ! [kg/m^2/s]
end function sublr_amorph



real(8) function alpha(T)
  ! condensation coefficient from Haynes et al. (1992)
  implicit none
  real(8), intent(IN) :: T
  !alpha = 1.06/(1.+ 1.*exp(-0.23/(0.00198588*T)))
  alpha = ( 1.06*(185-T) + 0.65*(T-20) ) / (185-20)
  if (alpha>1.) alpha=1.
  if (alpha<0.) alpha=0.
end function alpha



real(8) function padsr(x)
  ! relative vapor pressure of adsorbed H2O
  ! based on adsorption isotherms by Cadenhead & Stetter (1974)
  ! see Schorghofer & Aharonson (2014)
  implicit none
  real(8), intent(IN) :: x
  real(8) b, c
  real(8), parameter :: a=0.402, x0=2.45
  
  if (x<x0) then
     b = (exp(-a*x0)*(2+a*x0) - 2)/x0**3
     c = (3 - exp(-a*x0)*(3+a*x0))/x0**2
     padsr = b*x**3 + c*x**2
  else
     padsr = 1 - exp(-a*x)
  endif
  ! 0<= padsr <=1
end function padsr


