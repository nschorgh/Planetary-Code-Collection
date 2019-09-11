module constants
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8
  real(8), parameter :: So=1365.  ! solar constant [W/m^2]
end module constants


module body
  implicit none
  real(8) semia ! semimajor axis [AU]
  real(8) ecc   ! orbital eccentricity
  real(8) Trot  ! length of solar day in Earth days
  real(8) emiss ! IR emissivity
  real(8) albedo

  ! default
  !parameter(semia = 3.2, ecc = 0., Trot=10./24.)
  parameter(emiss = 0.96)
  !parameter(albedo = 0.05)

  ! Ceres
  parameter(semia = 2.76750591, ecc = 0.07582, Trot = 9.076/24.)  ! synodic
  parameter(albedo = 0.09)
  
  ! 133 P/Elst-Pizarro
  !parameter(semia = 3.1609, ecc=0.1617, Trot=3.471/24.)

  ! Trojan
  !parameter(semia = 5.2, ecc = 0.)
  
  ! Cybele
  !parameter (semia = 3.433; ecc=0.105, Trot=4.041/24.)
end module body


PROGRAM asteroid_thermal
  use constants, only : pi, d2r, So
  use body
  implicit none
  !integer n
  real(8) eps  ! obliquity [radians]
  real(8) omega  ! [radians]
  real(8) latitude, Qmean, Q4mean, Q0mean
  real(8) Tmean, thIn, Tmin, Tmax, Torb !, ti(13)
  real(8), external :: flux_noatm, flux2T, a2Torb

  ! Ceres
  eps = 4.*d2r; omega = 301.*d2r

  ! 133 P/Elst-Pizarro
  ! vernal equinox = 291.8+90=21.8
  ! eps = 75.*d2r; omega= (132.18-21.8)*d2r

  latitude = 0.*d2r

  print *,'a=',semia,'ecc=',ecc,'omega=',omega/d2r
  print *,'Latitude=',latitude/d2r,'obliquity=',eps/d2r
  write(*,'(1x,a7,1x,f5.3)') 'Albedo=',albedo

  ! orbital period (days)
  Torb = a2Torb(semia)
  print *,'Torb=',Torb,'days'

  print *,'Tss=',flux2T(flux_noatm(semia,0d0,0d0,0d0,0d0,0d0),albedo,emiss)
  print *,'Tss=',flux2T(So/semia**2,albedo,emiss)
  !print *, 'Tmax=',flux2T(So/(semia*(1-ecc))**2,albedo,emiss)

  ! Insolation - optional
  call insolonly(latitude,semia,omega,ecc,eps,Trot,Q0mean,Qmean,Q4mean)
  Qmean = (1-albedo)*Qmean;  Q4mean = (1-albedo)*Q4mean
  write(*,'(a,4(1x,f5.1))') 'Fluxes (W/m^2):',So/semia**2,Q0mean,(1-albedo)*Q0mean/pi,Qmean
  write(*,'(a,2(1x,f6.2))') 'End-member temperatures (K):', &
       & flux2T(Qmean,1d0,emiss),flux2T(Q4mean,1d0,emiss)

  ! Surface temperature with rotation and conduction
  !ti = (/ 10000., 2100., 2000., 1000., 500., 200., 100., 50., 25., 15., 10., 5., 3. /)
  !do n= 1,size(ti)
  !   thIn = ti(n)
  thIn = 15.
  call oneasteroid(latitude,omega,eps,thIn,Qmean,Tmean,Tmin,Tmax)
  print *, 'Mean insolation=',Qmean,'W/m^2'
  print *, 'Mean temperature=',Tmean,'K'
  print *, '#',thIn,Tmean,Tmin,Tmax
  !enddo
  print *

END PROGRAM asteroid_thermal

