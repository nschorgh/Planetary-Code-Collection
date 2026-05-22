module body
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.

  type orbitp
     real(8) semia    ! semimajor axis [AU]
     real(8) ecc      ! orbital eccentricity
     real(8) eps      ! axis tilt (obliquity) [radians]
     real(8) omega    ! Ls of perihelion relative to equinox [radians]
     real(8) solarDay ! length of solar day [seconds]
  end type orbitp

  real(8) emiss  ! IR emissivity of ice-free surface
  real(8) albedo
  integer nz     ! number of vertical grid points
  real(8) zfac   ! geometric increase in grid spacing
  real(8) zmax   ! domain depth

  type(orbitp) :: Orbit = &
       & orbitp(5.2, 0.07, 20.*d2r, 0., 12.*3600)  ! nominal Trojan
       !& orbitp(5.20, 0.044, 158.*d2r, 0., 8.702724*3600.)  ! Eurybates  
       !& orbitp(5.17, 0.057, 20.*d2r, 0., 11.5*3600.) ! Polymele
       !& orbitp(5.29, 0.024,  10.*d2r, 0., 445.683*3600.)  ! Leucus  
       !& orbitp(5.13, 0.013, 154.*d2r, 0., 13.48617*3600.)  ! Orus
       !& orbitp(5.22, 0.101, 20.*d2r, 0., 102.784*3600.) ! Patroclus
  
  parameter(albedo = 0.05)   ! nominal
  !parameter(albedo = 0.044)  ! Eurybates
  !parameter(albedo = 0.092)  ! Polymele  
  !parameter(albedo = 0.043)  ! Leucus
  !parameter(albedo = 0.040)  ! Orus
  ! none for Patroclus

  parameter(emiss = 0.90d0)
  
  parameter(nz=160, zfac=1.05d0, zmax=20.)  ! thIn=20
  !parameter(nz=100, zfac=1.05d0, zmax=1.)  ! without seasons
  real(8), parameter :: Tnominal = 120.  ! for initializations
  real(8), parameter :: diam = 100e-6  ! grain diameter [m]
  
  real(8), parameter :: dt = 0.01  ! [solar days]
  real(8), parameter :: Fgeotherm = 0.
  integer, parameter :: EQUILTIME = 5 ! [orbits]
end module body



subroutine outputmoduleparameters
  use body
  implicit none
  real(8), external :: sols_per_orbit
  
  print *,'Global parameters stored in modules'
  !print *,'  Ice bulk density',icedensity,'kg/m^3'
  print *,'  dt=',dt,'solar days'
  print *,'  Fgeotherm=',Fgeotherm,'W/m^2'
  print *,'  Albedo=',albedo
  print *,'  Emissivity=',emiss
  print *,'  Thermal model equilibration time',EQUILTIME,'orbits'
  print *,'  Semimajor axis',orbit%semia
  print *,'  Eccentricity',orbit%ecc
  print *,'  Obliquity',orbit%eps / d2r
  print *,'  Solar day',orbit%solarDay, &
       & 'Sols per orbit',sols_per_orbit( orbit%semia, orbit%solarDay )
  print *,'  Vertical grid: nz=',nz,' zfac=',zfac,'zmax=',zmax
  print *,'  Grain diameter diam=',diam
end subroutine outputmoduleparameters
