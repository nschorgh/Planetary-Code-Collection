!************************************************************************
! thermophysical properties
!************************************************************************

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
  !c = 8.9093E-09*T**4 -1.2340E-05*T**3 +2.36160E-03*T**2 + 2.7431*T -3.6125
  !heatcapacity = c

  ! Biele et al., Int. J. Thermophys. 43, 144 (2022), eq. 24
  real(8) x
  x = log(T)
  c = exp((3.*x**3 -54.45*x**2 +306.8*x -376.6)/(x**2 -16.81*x +87.32))
  heatcapacity = c
end function heatcapacity


elemental function meanfreepathinsoil(diam,porosity)
  ! mean-free path projected on vertical, ell [m]
  implicit none
  real(8) meanfreepathinsoil
  real(8), intent(IN) :: diam, porosity

  meanfreepathinsoil = diam/sqrt(2.)*porosity
end function meanfreepathinsoil


subroutine assignthermalproperties3(nz,z,T,porosity,ti,rhocv,icefrac,zdepthT)
  ! assign thermal properties of soil
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), T(nz), porosity(nz)
  real(8), intent(OUT) :: ti(nz), rhocv(nz)
  real(8), intent(IN), optional :: icefrac(nz), zdepthT
  real(8), parameter :: rhodry = 2500  ! [kg/m^3]
  real(8), external :: heatcapacity

  integer j
  real(8) icedensity
  real(8) cice  ! heat capacity of pure ice
  real(8) kice  ! thermal conductivity of pure ice
  real(8) k(nz) ! thermal conductivity
  real(8) thIn, Tday, buf

  thIn = 20.
  !thIn = 100.
  do j=1,nz
     rhocv(j) = (1.-porosity(j)) * rhodry * heatcapacity(T(j))
     
     Tday = 150. ! estimated temperature when thIn was measured
     buf = (1.-porosity(j)) * rhodry * heatcapacity(Tday)
     k(j) = thIn**2/buf
  end do
  if (present(icefrac) .and. present(zdepthT)) then
     if (zdepthT>=0.) then
        do j=1,nz
           if (z(j)>zdepthT) then
              call iceproperties_species('H2O',T(j),icedensity,cice,kice)
              ! linear addition of conductivities in the spirit of Siegler et al. (2012)
              k(j) = k(j) + icefrac(j)*kice  ! icefrac <= porosity
              rhocv(j) = rhocv(j) + icedensity*cice*icefrac(j)
           endif
        end do
     end if
  end if
  
  ti(1:nz) = sqrt(k(1:nz)*rhocv(1:nz))
  print *,'min thIn=',minval(ti)
end subroutine assignthermalproperties3


subroutine iceproperties_species(species,T,icedensity,cice,kice)
  ! thermal properties of water ice and other ices
  implicit none
  character(*), intent(IN) :: species
  real(8), intent(IN) :: T
  real(8), intent(OUT) :: icedensity ! [kg/m^3]
  real(8), intent(OUT) :: cice  ! heat capacity of pure ice
  real(8), intent(OUT) :: kice  ! thermal conductivity of pure ice

  select case (species)
  case('H2O')
     icedensity = 933.  ! at 120K
     ! cice = 7.8*T
     cice = 7.49*T + 90.  ! Klinger (1981), Shulman (2004)
     !kice = 632./T+0.38-1.97e-3*T
     kice = 612./T  ! DOI:10.17632/ttzbgxs9fw.2

  case('HCN')  ! hydrogen cyanide
     icedensity = 1037.  ! Gerakines et al. (2022)
     ! 8.938 cal/K/mol at 120K (Giauque & Ruehrwein 1939)
     cice = 1384. ! 8.938*4184./27.02 J/K/kg
     kice = 3.  ! (assumed)

  case('CO2')  ! carbon dioxide
     icedensity = 1680.  ! approx. from Yu et al. (2023)
     cice = 951. ! 10*4184./44.  for more detail see Giauque & Egan (1937)
     kice = 0.7  ! approx. from Saiduzzaman et al. (2025), Fig. 8
     
  case default
     error stop ('iceproperties_species: no species matches')
     
  end select
  ! Note: icedensity should match number in subroutine icechanges3
end subroutine iceproperties_species
  

elemental function porosityprofile(z)
  implicit none
  real(8) porosityprofile
  real(8), intent(IN) :: z
  real(8) phi_s  ! surface porosity
  real(8) phi_d  ! porosity at great depth
  real(8) H ! depth-scale
  
  ! rho = (1-porosity)*rhosolid,  porosity = 1-rho/rhosolid
  ! twolayers = rho_d - (rho_d - rho_s) * exp(-z/H)
  
  ! values inspired by Fa, Earth and Space Science 7, e2019EA000801 (2020)
  H = 1.03
  phi_s = 0.745   ! surface porosity
  ! porosity at 5m depth is 0.323
  phi_d = 0.33  ! 0.326 at z->infinity
  porosityprofile = 0.33 + (0.66-0.33)*exp(-z/H)
end function porosityprofile
