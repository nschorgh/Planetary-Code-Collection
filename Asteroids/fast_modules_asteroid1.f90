!*****************************************************
! All modules for fast asteroid method
!*****************************************************


module constants
  ! miscellaneous parameters that are very constant
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer, parameter :: NMAX=1000
  real(8), parameter :: sigSB=5.6704e-8
  real(8), parameter :: kB=1.38065e-23
end module constants



module body
  implicit none

  real(8) semia  ! semimajor axis (AU)
  real(8) ecc  ! orbital eccentricity 
  real(8) Trot ! length of solar day in Earth days
  real(8) solarDay ! length of solar day in seconds
  real(8) emiss ! IR emissivity of dry surface
  real(8) solsperyear
  integer nz   ! number of vertical grid points
  real(8) zfac
  real(8) zmax  ! domain depth

  ! (1) Ceres
  parameter(semia = 2.76750591)
  !parameter(ecc = 0.075822766)  ! current
  parameter(ecc = 0.0)
  !parameter(ecc = 0.117)  ! proper
  parameter(Trot = 9.074170/24., solarDay = 9.076*3600.)
  parameter(solsperyear = 4446.)
  parameter(emiss = 0.95d0)
  parameter(nz=160, zfac=1.05d0, zmax=15.)  ! thIn=15


  real(8), parameter :: Tnominal = 140.   ! for Diff and Tinit
  real(8), parameter :: icedensity = 931.  ! 140K

  real(8), parameter :: dt = 0.01  ! in units of solar days
  real(8), parameter :: Fgeotherm = 0.
  integer, parameter :: EQUILTIME = 10 ! (orbits)
end module body



module allinterfaces
  ! interfaces from Fortran 90 subroutines and functions

  ! fast_subs_asteroid1.f90
  interface
     elemental function faintsun(t)
       implicit none
       real(8) faintsun
       real(8), intent(IN) :: t   ! time before present [years]
     end function faintsun
  end interface

  interface
     subroutine icelayer_asteroid(bigstep,NP,z,porosity,icefrac,Tinit, &
          & zdepthT,Tmean1,Tmean3,Tmin,Tmax,latitude,albedo,ecc,omega,eps,S0)
       use constants, only : d2r, NMAX
       use body, only : icedensity, Tnominal, nz
       implicit none
       integer, intent(IN) :: NP
       real(8), intent(IN) :: bigstep
       real(8), intent(IN) :: z(NMAX), porosity, icefrac
       real(8), intent(INOUT) :: zdepthT(NP), Tmean1(NP), Tmean3(NP)
       real(8), intent(OUT) :: Tmin(NP), Tmax(NP)
       real(8), intent(IN) :: latitude(NP), albedo(NP), ecc, omega, eps, S0
       logical, intent(IN) :: Tinit
     end subroutine icelayer_asteroid
  end interface

  interface
     subroutine ajsub_asteroid(latitude, albedo, z, ti, rhocv, ecc, omega, eps, &
          &     S0, typeT, avrho, Tinit, Tmean1, Tmean3, Tmin, Tmax)
       use constants
       use body, only : EQUILTIME, dt, solsperyear, Fgeotherm, semia, nz, emiss, solarDay
       implicit none
       real(8), intent(IN) :: latitude
       real(8), intent(IN) :: albedo, z(NMAX)
       real(8), intent(IN) :: ti(NMAX), rhocv(NMAX)
       real(8), intent(IN) :: ecc, omega, eps, S0
       integer, intent(IN) :: typeT
       real(8), intent(OUT) :: avrho
       logical, intent(IN) :: Tinit
       real(8), intent(INOUT) :: Tmean1, Tmean3
       real(8), intent(OUT) :: Tmin, Tmax
     end subroutine ajsub_asteroid
  end interface

  interface
     subroutine icechanges(nz,z,typeT,avrho,Deff,bigstep,zdepthT,porosity,icefrac)
       implicit none
       integer, intent(IN) :: nz, typeT
       real(8), intent(IN) :: z(nz), avrho
       real(8), intent(IN) :: Deff, bigstep, porosity, icefrac
       real(8), intent(INOUT) :: zdepthT
     end subroutine icechanges
  end interface

  interface
     pure function zint(y1,y2,z1,z2)
       implicit none
       real(8), intent(IN) :: y1,y2,z1,z2
       real(8) zint
     end function zint
  end interface

  interface ! moved to grids.f 
     pure function colint(y,z,nz,i1,i2)
       implicit none
       integer, intent(IN) :: nz, i1, i2
       real(8), intent(IN) :: y(nz),z(nz)
       real(8) colint
     end function colint
  end interface

  interface
     subroutine assignthermalproperties1(nz,z,Tnom,porosity,ti,rhocv,icefrac,zdepthT)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz), Tnom, porosity
       real(8), intent(OUT) :: ti(nz), rhocv(nz)
       real(8), intent(IN), optional :: icefrac, zdepthT
     end subroutine assignthermalproperties1
  end interface

  interface
     elemental function heatcapacity(T)
       implicit none
       real(8), intent(IN) :: T
       real(8) heatcapacity
     end function heatcapacity
  end interface

  interface
     elemental function conductivity(T)
       implicit none
       real(8), intent(IN) :: T
       real(8) conductivity
     end function conductivity
  end interface

  interface
     function vapordiffusivity(diam,porosity,T)
       implicit none
       real(8) vapordiffusivity
       real(8), intent(IN) :: diam,porosity,T
     end function vapordiffusivity
  end interface

  interface
     subroutine dzvector(nz,z,dz) 
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz)
       real(8), intent(OUT) :: dz(nz)
     end subroutine dzvector
  end interface

  interface
     integer function gettype(zdepth,nz,z)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: zdepth, z(nz)
     end function gettype
  end interface

  ! Common/*.f90
  interface
     elemental real(8) function flux_noatm(R,decl,latitude,HA,surfaceSlope,azFac)
       implicit none
       real(8), intent(IN) :: R,decl,latitude,HA,surfaceSlope,azFac
     end function flux_noatm
  end interface

end module allinterfaces


