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
  real(8), parameter :: rhodry = 2500  ! bulk density
  real(8), parameter :: icedensity = 933.  ! 120K  [kg/m^3]
  real(8), external :: heatcapacity

  integer j
  real(8) cice  ! heat capacity of pure ice
  real(8) kice  ! thermal conductivity of pure ice
  real(8) k(nz) ! thermal conductivity
  real(8) thIn

  thIn = 20.
  !thIn = 100.
  do j=1,nz
     rhocv(j) = (1.-porosity(j)) * rhodry * heatcapacity(T(j))
     k(j) = thIn**2/rhocv(j)
  end do
  if (present(icefrac) .and. present(zdepthT)) then
     do j=1,nz
        if (z(j)>zdepthT) then
           ! cice = 7.8*T(j) 
           cice = 7.49*T(j) + 90.  ! Klinger (1981), Shulman (2004)
           !kice = 632./T(j)+0.38-1.97e-3*T(j)
           kice = 612./T(j)  ! DOI:10.17632/ttzbgxs9fw.1
           k(j) = k(j) + icefrac(j)*kice  ! icefrac <= porosity
           rhocv(j) = rhocv(j) + icedensity*cice*icefrac(j)
        endif
     end do
  end if

  ti(1:nz) = sqrt(k(1:nz)*rhocv(1:nz))
  print *,'min thIn=',minval(ti)
end subroutine assignthermalproperties3
