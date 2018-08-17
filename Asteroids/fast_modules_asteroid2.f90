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
  
  real(8), parameter :: Tnominal = 140.   ! for Diff and Tinit
  real(8), parameter :: icedensity = 931.  ! 140K

  real(8), parameter :: dt = 0.01  ! in units of solar days
  real(8), parameter :: Fgeotherm = 0.
  integer, parameter :: EQUILTIME = 20 ! (orbits)
end module body



module allinterfaces
  ! interfaces from Fortran 90 subroutines and functions

  ! fast_subs_asteroid2.f90
  interface
     elemental function faintsun(t)
       implicit none
       real(8) faintsun
       real(8), intent(IN) :: t   ! time before present [years]
     end function faintsun
  end interface

  interface
     subroutine icelayer_asteroid(bigstep,NP,z,porosity,Tinit, &
          & zdepthP,sigma,Tmean1,Tmean3,Tmin,Tmax,latitude,albedo,ecc,omega,eps,S0)
       use constants, only : d2r, NMAX
       use body, only : icedensity, Tnominal, nz
       implicit none
       integer, intent(IN) :: NP
       real(8), intent(IN) :: bigstep
       real(8), intent(IN) :: z(NMAX), porosity(nz)
       real(8), intent(INOUT) :: sigma(nz,NP), zdepthP(NP), Tmean1(NP), Tmean3(NP)
       real(8), intent(OUT) :: Tmin(NP), Tmax(NP)
       real(8), intent(IN) :: latitude(NP), albedo(NP), ecc, omega, eps, S0
       logical, intent(IN) :: Tinit
     end subroutine icelayer_asteroid
  end interface

  interface
     subroutine ajsub_asteroid(latitude, albedo, z, ti, rhocv, ecc, omega, eps, &
          &     S0, typeP, Diff, Diff0, avrho, Tinit, ypp, Jp, Tmean1, Tmean3, Tmin, Tmax)
       use constants
       use body, only : EQUILTIME, dt, solsperyear, Fgeotherm, semia, nz, emiss, solarDay
       implicit none
       real(8), intent(IN) :: latitude
       real(8), intent(IN) :: albedo, z(NMAX)
       real(8), intent(IN) :: ti(NMAX), rhocv(NMAX)
       real(8), intent(IN) :: ecc, omega, eps, Diff(nz), Diff0, S0
       integer, intent(IN) :: typeP
       real(8), intent(OUT) :: avrho(nz)
       logical, intent(IN) :: Tinit
       real(8), intent(OUT) :: ypp(nz), Jp
       real(8), intent(INOUT) :: Tmean1, Tmean3
       real(8), intent(OUT) :: Tmin, Tmax
     end subroutine ajsub_asteroid
  end interface

  interface
     subroutine avmeth(nz, z, rhosatav, rhosatav0, rlow, typeP, Diff, Diff0, ypp, Jpump1)
       implicit none
       integer, intent(IN) :: nz, typeP
       real(8), intent(IN), dimension(nz) :: z, rhosatav, Diff
       real(8), intent(IN) :: rhosatav0, rlow, Diff0
       real(8), intent(OUT) :: ypp(nz), Jpump1
     end subroutine avmeth
  end interface

  interface
     subroutine icechanges(nz,z,typeP,avrho,ypp,Deff,bigstep,Jp,zdepthP,sigma)
       implicit none
       integer, intent(IN) :: nz, typeP
       real(8), intent(IN) :: z(nz), ypp(nz), avrho(nz)
       real(8), intent(IN) :: Deff, bigstep, Jp
       real(8), intent(INOUT) :: zdepthP, sigma(nz)
     end subroutine icechanges
  end interface

  interface
     pure function zint(y1,y2,z1,z2)
       implicit none
       real(8), intent(IN) :: y1,y2,z1,z2
       real(8) zint
     end function zint
  end interface

  interface
     subroutine compactoutput(unit,sigma,nz)
       implicit none
       integer, intent(IN) :: unit,nz
       real(8), intent(IN) :: sigma(nz)
     end subroutine compactoutput
  end interface

  interface
     subroutine assignthermalproperties(nz,Tnom,porosity,ti,rhocv,porefill)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: Tnom, porosity(nz)
       real(8), intent(OUT) :: ti(nz), rhocv(nz)
       real(8), intent(IN), optional :: porefill(nz)
     end subroutine assignthermalproperties
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
     elemental function constriction(porefill)
       implicit none
       real(8), intent(IN) :: porefill
       real(8) eta, constriction
     end function constriction
  end interface
  
  ! impactstirring.f90
  interface
     subroutine impactstirring(nz,z,dt,rho)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz), dt 
       real(8), intent(INOUT) :: rho(nz)
     end subroutine impactstirring
  end interface

  interface
     real(8) function pareto(idum,mean)
       implicit none
       integer, intent(INOUT) :: idum
       real(8), intent(IN) :: mean
       real(8), external :: ran2
     end function pareto
  end interface

  interface
     integer function gettype(zdepth,nz,z)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: zdepth, z(nz)
     end function gettype
  end interface

  ! Common/{*.f,*.f90}
  interface
     elemental real(8) function flux_noatm(R,decl,latitude,HA,surfaceSlope,azFac)
       implicit none
       real(8), intent(IN) :: R,decl,latitude,HA,surfaceSlope,azFac
     end function flux_noatm
  end interface

  interface
     pure function colint(y,z,nz,i1,i2)
       implicit none
       integer, intent(IN) :: nz, i1, i2
       real(8), intent(IN) :: y(nz),z(nz)
       real(8) colint
     end function colint
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
end module allinterfaces


