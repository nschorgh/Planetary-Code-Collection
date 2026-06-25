elemental function psv_species(T,species)
  ! vapor pressures of various ices other than H2O
  ! most coefficients are from Schorghofer & Williams, Icarus 416, 116086 (2024)
  implicit none
  real(8) psv_species  ! [Pa]
  real(8), intent(IN) :: T  ! [Kelvin]
  character(*), intent(IN) :: species
  real(8) b0, b1, b2

  select case (species)
     case('CO2')
  
        !real(8), parameter :: A=6.81228, B=1301.679, C=-3.494  ! NIST Webbook
        !psv_species = 10**(A-(B/(T+C)))*1.e5   ! Antoine equation

        ! Fray & Schmitt (2009)
        !real(8), parameter :: A0=1.861e1, A1=-4.154e3, A2=1.041e5
        !psv = exp(A0+A1/T+A2/T2**2)*1e5
        ! 40<T<194.7
        !real(8), parameter :: A0=1.476e1, A1=-2.571e3, A2=-7.781e4
        !real(8), parameter :: A3=4.325e6, A4=-1.207e8, A5=1.350e9
        !psv_species = exp(A0+A1/T+A2/T**2+A3/T**3+A4/T**4+A5/T**5)*1e5

        !block
        ! Bryson et al. (1974) after unit conversions
        !  real(8), parameter :: A=28.693, B=3270.9 ! 91-190K
        !  psv_species = exp(A-B/T)
        !end block
        
        b0=32.61; b1=3291; b2=-0.7947
        psv_species = exp(b0-b1/T+b2*log(T))
        
     case('NH3') ! ammonia
        b0=28.71; b1=3903
        psv_species = exp(b0-b1/T)
        
     case('SO2')
        !b0=28.92; b1=4262
        b0=9.435; b1=3775; b2=3.225
        psv_species = exp(b0-b1/T+b2*log(T))
        
     case('HCN') ! hydrogen cyanide
        if (T>170.) then ! phase I
           b0=27.03; b1=4472.
        else ! phase II (Hudson & Gerakines, 2023)
           b0=26.50; b1=4568.
        end if
        ! Note: the two do not connect continuously at phase boundary
        psv_species = exp(b0-b1/T)
        
     case('CH3OH') ! methanol
        if (T<157.4) then ! alpha-phase
           b0=15.94; b1=2453
        else ! beta-phase
           b0=15.02; b1=2308
        end if
        psv_species = exp(b0-b1/T)
        
     case default
        error stop ('psv_species: no species matches')
        
  end select

end function psv_species
