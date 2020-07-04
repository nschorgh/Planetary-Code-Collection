program oscidea1
  ! Lunar thermal ice pump
  ! written by Norbert Schorghofer, 2013
  use oscidea_params
  implicit none
  !real(8), parameter :: unitconv = m*365.24*86400 ! #/m^2/s -> kg/m^2/year
  real(8) Tm, Ta, Tmean, Esurf, Ebase, supply, spaceweather0, mintheta, maxtheta
  real(8) avtheta, fract, avweather !, convergence
  real(8), external :: sublr_amorph

  supply = 1e12
  spaceweather0 = 1e12 
  ! 1 m/Ga = 1e12 molecules/s/m^2, smooth plane
  ! 1000/(1e9*365.24*86400)/(18*1.66e-27)
  
  Tm=80.; Ta=120.

  call surface_integration(Tm,Ta,supply,spaceweather0, &
       & Tmean,avtheta,fract,Esurf,avweather,mintheta,maxtheta)
  Ebase = sublr_amorph(Tmean)

  write(*,'(a,e10.3,1x,a,e10.3)') &
       & 'Supply=',supply,'Spaceweather0=',spaceweather0
  write(*,'(2(a3,1x,f6.2,1x),a6,e10.3,a10,f10.3)') &
       & 'Tm=',Tm,'Ta=',Ta,'theta=',avtheta,'fractime=',fract

  !convergence = (rough*(avweather+Esurf)-supply)/supply  ! should be ~0
  !write(*,'(4(1x,e10.3))') Esurf,Ebase,avweather,avtheta
  !write(*,'(2(f6.2,1x),e10.3)') Tm,Ta,Esurf-Ebase
  write(*,'(f6.2,1x,g10.4,3(1x,g10.3))') &
       & Tmean,Esurf/supply,avtheta/theta0,mintheta/theta0,maxtheta/theta0
end program oscidea1

