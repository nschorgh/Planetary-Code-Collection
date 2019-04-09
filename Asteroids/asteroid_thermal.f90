module constants
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8
  real(8), parameter :: So=1365.  ! solar constant [W/m^2]
end module constants


module body
  implicit none
  real(8) a   ! semimajor axis [AU]
  real(8) ecc  ! orbital eccentricity
  real(8) Trot ! length of solar day in Earth days
  real(8) emiss ! IR emissivity

  ! default
  !parameter(a = 3.2, ecc = 0., Trot=10./24.)
  parameter(emiss=0.96)  

  ! Ceres
  parameter(a = 2.76750591, ecc = 0.07582, Trot = 9.076/24.)  ! synodic

  ! 133 P/Elst-Pizarro
  !parameter(a = 3.1609, ecc=0.1617, Trot=3.471/24.)

  ! P/2005 U1 (Read)
  !parameter(a = 3.1649, ecc=0.2528)

  ! 118401 (1999 RE70)
  !parameter(a = 3.1929, ecc=0.1933)

  ! Trojan
  !parameter(a = 5.2, ecc = 0.)
  
  ! Cybele
  !parameter (a = 3.433; ecc=0.105, Trot=4.041/24.)
end module body


program asteroid
  use constants, only : pi, d2r, So
  use body
  implicit none
  !integer n
  real(8) eps  ! obliquity [radians]
  real(8) omega  ! [radians]
  real(8) latitude, Qmean, Q4mean, Q0mean
  real(8) Tmean, thIn, Tmin, Tmax, Torb !, ti(13)
  real(8) albedo
  real(8), external :: flux_noatm, flux2T, a2Torb

  ! Ceres
  eps = 4.*d2r; omega = 301.*d2r
  ! albedo = 0.034  ! Bond albedo
  albedo = 0.09

  ! 133 P/Elst-Pizarro
  ! vernal equinox = 291.8+90=21.8
  ! eps = 75.*d2r; omega= (132.18-21.8)*d2r

  latitude = 0.*d2r

  print *,'a=',a,'ecc=',ecc,'omega=',omega/d2r
  print *,'Latitude=',latitude/d2r,'obliquity=',eps/d2r
  write(6,'(1x,a7,1x,f5.3)') 'Albedo=',albedo

  ! orbital period (days)
  Torb = a2Torb(a)
  print *,'Torb=',Torb,'days'

  write(6,*) 'Tss=',flux2T(flux_noatm(a,0d0,0d0,0d0,0d0,0d0),albedo,emiss)
  write(6,*) 'Tss=',flux2T(So/a**2,albedo,emiss)
  !write(6,*) 'Tmax=',flux2T(So/(a*(1-ecc))**2,albedo,emiss)

  ! Insolation - optional
  call insolonly(latitude,a,omega,ecc,eps,Trot,Q0mean,Qmean,Q4mean)
  Qmean = (1-albedo)*Qmean;  Q4mean = (1-albedo)*Q4mean
  write(6,'(a,4(1x,f5.1))') 'Fluxes (W/m^2):',So/a**2,Q0mean,(1-albedo)*Q0mean/pi,Qmean
  write(6,'(a,2(1x,f6.2))') 'End-member temperatures (K):', &
       & flux2T(Qmean,1d0,emiss),flux2T(Q4mean,1d0,emiss)

  ! Surface temperature with rotation and conduction
  !ti = (/ 10000., 2100., 2000., 1000., 500., 200., 100., 50., 25., 15., 10., 5., 3. /)
  !do n= 1,size(ti)
  !   thIn = ti(n)
  thIn = 15.
  call oneasteroid(latitude,omega,eps,albedo,thIn,Qmean,Tmean,Tmin,Tmax)
  write(6,*) 'Mean insolation=',Qmean,'W/m^2'
  write(6,*) 'Mean temperature=',Tmean,'K'
  write(6,*)  '#',thIn,Tmean,Tmin,Tmax
  !enddo
  write(6,*)

end program asteroid





