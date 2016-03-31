module body
  ! physical properties of the body 

  ! Moon
  real(8), parameter :: solarDay=86400.*29.53  ! (s)
  real(8), parameter :: siderealDay=655.728*3600.
  real(8), parameter :: g=1.62, Rmoon=1737.e3
  real(8), parameter :: dtsec=3600.  ! thermal model time step (s)
  real(8), parameter :: semia=1.  ! (AU)
  real(8), parameter :: zmax=0.5   ! domain depth for 1D thermal model, if used
  real(8), parameter :: albedo=0.11, emiss=0.95
  real(8), parameter :: Fgeotherm=0.018  ! Langseth et al. (1976)  (W/m^2)

  ! Mercury
  !real(8), parameter :: solarDay=4222.6*3600., semia=0.3871
  !real(8), parameter :: siderealDay=1407.6*3600.
  !real(8), parameter :: g=3.7, Rmoon=2440e3
  !real(8), parameter :: dtsec=3600.*12.
  !real(8), parameter :: zmax=1.5
  !real(8), parameter :: albedo=0.07, emiss=0.9
  !real(8), parameter :: Fgeotherm=0.

  ! Ceres
  !real(8), parameter :: solarDay=9.076*3600., semia=2.77
  !real(8), parameter :: siderealDay=9.07417*3600.
  !real(8), parameter :: g=0.27, Rmoon=480e3  !(975x909km)
  !real(8), parameter :: dtsec=600. 
  !real(8), parameter :: albedo=0.09, emiss=0.9
  !real(8), parameter :: zmax=0.5, Fgeotherm=0.

  real(8), parameter :: vescape=sqrt(2*g*Rmoon)
end module body
