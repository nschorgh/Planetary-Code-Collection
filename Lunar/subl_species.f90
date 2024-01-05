elemental function evap_vacuum_species(T,species)
  ! sublimation rate of various ices into vacuum
  implicit none
  real(8) evap_vacuum_species  ! [kg/m^2/s]
  real(8), intent(IN) :: T  ! [Kelvin]
  character(*), intent(IN) :: species
  real(8), parameter :: pi=3.1415926535897932
  real(8), parameter :: Ru=8314.46  ! [J/kmol/K]
  real(8) psv, mu
  real(8) b0, b1, b2

  select case (species)
     case('CO2')
        mu = 44.01
  
        !real(8), parameter :: A=6.81228, B=1301.679, C=-3.494  ! NIST Webbook
        !psv = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation

        ! Fray & Schmitt (2009)
        !real(8), parameter :: A0=1.861e1, A1=-4.154e3, A2=1.041e5  ! 194.7<T<216.58
        !psv = exp(A0+A1/T+A2/T2**2)*1e5
        ! 40<T<194.7
        !real(8), parameter :: A0=1.476e1, A1=-2.571e3, A2=-7.781e4
        !real(8), parameter :: A3=4.325e6, A4=-1.207e8, A5=1.350e9
        !psv = exp(A0+A1/T+A2/T**2+A3/T**3+A4/T**4+A5/T**5)*1e5

        !block
        ! Bryson et al. (1974) after unit conversions
        !  real(8), parameter :: A=28.693, B=3270.9 ! 91-190K
        !  psv = exp(A-B/T)
        !end block
        
        ! new parametrization - Nov 2023
        b0=32.6; b1=3292; b2=-0.08
        psv = exp(b0-b1/T+b2*log(T))
        
     case('NH3') ! ammonia
        mu = 17.03
        b0=28.7; b1=3903
        psv = exp(b0-b1/T)
        
     case('SO2')
        mu = 64.06
        !b0=29.9; b1=4262
        b0=9; b1=3775; b2=3
        psv = exp(b0-b1/T+b2*log(T))
        
     case('HCN')
        mu = 27.0253
        b0=27.03; b1=4472
        psv = exp(b0-b1/T)
        
     case('CH3OH')
        mu = 32.0419
        if (T<157.4) then ! alpha-phase
           b0=15.94; b1=2453
        else ! beta-phase
           b0=15.02; b1=2308
        end if
        psv = exp(b0-b1/T)
        
     case default
        error stop ('evap_vacuum_species: no species matches')
        
  end select

  ! calculate sublimation rate from vapor pressure
  evap_vacuum_species = psv*sqrt(mu/(2*pi*Ru*T))
  
end function evap_vacuum_species
