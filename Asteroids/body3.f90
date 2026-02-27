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
       & orbitp(5.2, 0., 20.*d2r, 0., 12.*3600)  ! nominal Trojan
       !& orbitp(5.20, 0.089, 158.*d2r, 0., 8.702724*3600.)  ! Eurybates  
       !& orbitp(semia=5.17, ecc=0.095, solarDay=11.5*3600.) ! Polymele
       !& orbitp(5.29, 0.064,  10.*d2r, 445.683*3600.)  ! Leucus  
       !& orbitp(5.13, 0.037, 154.*d2r, 13.48617*3600.)  ! Orus
       !& orbitp(semia=5.22, ecc=0.129, solarDay = 102.784*3600.) ! Patroclus
  !Orbit%eps = 0.; Orbit%omega=0.
  
  parameter(emiss = 0.90d0, albedo = 0.05)

  parameter(nz=160, zfac=1.05d0, zmax=20.)  ! thIn=20
  !parameter(nz=100, zfac=1.05d0, zmax=1.)  ! without seasons
  real(8), parameter :: Tnominal = 120.  ! for initializations
  real(8), parameter :: diam = 100e-6  ! grain diameter [m]
  
  real(8), parameter :: dt = 0.01  ! [solar days]
  real(8), parameter :: Fgeotherm = 0.
  integer, parameter :: EQUILTIME = 10 ! [orbits]
end module body



subroutine outputmoduleparameters
  use body
  !use allinterfaces, only : sols_per_orbit
  implicit none
  real(8), external :: sols_per_orbit
  
  print *,'Global parameters stored in modules'
  !print *,'  Ice bulk density',icedensity,'kg/m^3'
  print *,'  dt=',dt,'solar days'
  print *,'  Fgeotherm=',Fgeotherm,'W/m^2'
  print *,'  Emissivity of surface=',emiss
  print *,'  Thermal model equilibration time',EQUILTIME,'orbits'
  print *,'  Semimajor axis',orbit%semia
  print *,'  Eccentricity',orbit%ecc
  print *,'  Obliquity',orbit%eps
  print *,'  Solar day',orbit%solarDay, &
       & 'Sols per orbit',sols_per_orbit( orbit%semia, orbit%solarDay )
  print *,'  Vertical grid: nz=',nz,' zfac=',zfac,'zmax=',zmax
  print *,'  Grain diameter diam=',diam
end subroutine outputmoduleparameters
