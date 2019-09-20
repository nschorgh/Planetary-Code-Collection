!***********************************************************************
! All modules for fast asteroid method
!***********************************************************************


module constants
  ! miscellaneous parameters that are very constant
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer, parameter :: NMAX=1000
  real(8), parameter :: sigSB=5.6704e-8
  real(8), parameter :: kB=1.38065e-23
end module constants



module body
  implicit none

  real(8) semia  ! semimajor axis [AU]
  real(8) ecc  ! orbital eccentricity 
  real(8) Trot ! length of solar day in Earth days
  real(8) solarDay ! length of solar day [seconds]
  real(8) emiss ! IR emissivity of dry surface
  real(8) solsperyear
  integer nz   ! number of vertical grid points
  real(8) zfac
  real(8) zmax  ! domain depth

  ! (1) Ceres
  parameter(semia = 2.76750591)
  !parameter(ecc = 0.075822766)  ! current
  !parameter(ecc = 0.0)
  parameter(ecc = 0.117)  ! proper
  parameter(Trot = 9.074170/24., solarDay = 9.076*3600.)
  parameter(solsperyear = 4446.)
  parameter(emiss = 0.95d0)
  parameter(nz=160, zfac=1.05d0, zmax=15.)  ! thIn=15

  ! 133P/Elst-Pizarro
  !parameter(semia=3.163, ecc=0.159)  ! proper
  !parameter(Trot=3.471/24., solarDay=3.471*3600.)
  !parameter(solsperyear=14192)
  !parameter(emiss = 0.95d0)
  !parameter(nz=160, zfac=1.05d0, zmax=20.) 
  
  real(8), parameter :: Tnominal = 140.   ! for Diff and Tinit [K]
  real(8), parameter :: icedensity = 931.  ! 140K  [kg/m^3]

  real(8), parameter :: dt = 0.01  ! in units of solar days
  real(8), parameter :: Fgeotherm = 0.
  integer, parameter :: EQUILTIME = 20 ! (orbits)
end module body



module allinterfaces
  
  ! fast_subs_asteroid1.f90
  interface
     subroutine assignthermalproperties1(nz,z,Tnom,porosity,ti,rhocv,icefrac,zdepthT)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz), Tnom, porosity
       real(8), intent(OUT) :: ti(nz), rhocv(nz)
       real(8), intent(IN), optional :: icefrac, zdepthT
     end subroutine assignthermalproperties1
  end interface
  
  ! fast_subs_asteroid2.f90
  interface
     subroutine assignthermalproperties(nz,Tnom,porosity,ti,rhocv,porefill)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: Tnom, porosity(nz)
       real(8), intent(OUT) :: ti(nz), rhocv(nz)
       real(8), intent(IN), optional :: porefill(nz)
     end subroutine assignthermalproperties
  end interface

  ! Common/*.f90
  interface
     pure function flux_noatm(R,decl,latitude,HA,surfaceSlope,azFac)
       implicit none
       real(8) flux_noatm
       real(8), intent(IN) :: R,decl,latitude,HA,surfaceSlope,azFac
     end function flux_noatm
  end interface

  interface
     subroutine deriv1(z,nz,y,y0,yNp1,yp)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz),y(nz),y0,yNp1
       real*8, intent(OUT) :: yp(nz)
     end subroutine deriv1
  end interface

  interface
     subroutine deriv2_simple(z,nz,y,y0,yNp1,yp2)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz),y(nz),y0,yNp1
       real*8, intent(OUT) :: yp2(nz)
     end subroutine deriv2_simple
  end interface

  interface
     real(8) function deriv1_onesided(j,z,nz,y)
       implicit none
       integer, intent(IN) :: nz,j
       real(8), intent(IN) :: z(nz),y(nz)
     end function deriv1_onesided
  end interface

  ! Commmon/grids.f
  interface
     pure function colint(y,z,nz,i1,i2)
       implicit none
       integer, intent(IN) :: nz, i1, i2
       real(8), intent(IN) :: y(nz),z(nz)
       real(8) colint
     end function colint
  end interface

  interface
     pure subroutine dzvector(nz,z,dz) 
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz)
       real(8), intent(OUT) :: dz(nz)
     end subroutine dzvector
  end interface
  
end module allinterfaces


