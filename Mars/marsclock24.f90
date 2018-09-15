subroutine marsclock24(JDUT,Deltat_J2000,Ls,dec,RM,Longitude_W,LTST)
  !***********************************************************************
  ! functionality of marsorbit.f + local time
  ! AM2000 = Allison & McEwen, PSS 48, 215 (2000)
  ! see also  www.giss.nasa.gov/tools/mars24/help/algorithm.html
  !***********************************************************************
  implicit none
  real*8, intent(IN) :: JDUT !  Julian Date
  real*8 temp1, dcor
  real*8 JDTT ! terrestrial julian date
  real*8, intent(OUT) :: Deltat_J2000  ! days since J2000
  real*8, intent(OUT) :: Ls   ! [radian]
  real*8, intent(OUT) :: dec, RM

  integer i
  real*8 A(7), tau(7), phi(7) 
  real*8 M, alpha_FMS, PBS
  real*8 numinusM   ! (degree)
  data A/ 0.0071, 0.0057, 0.0039, 0.0037, 0.0021, 0.0020, 0.0018 /
  !data A/ 0.007d0,0.006d0,0.004d0,0.004d0,0.002d0,0.002d0,0.002d0 /
  data tau/2.2353d0, 2.7543d0, 1.1177d0, 15.7866d0, 2.1354d0, 2.4694d0, 32.8493d0/
  data phi/49.409d0, 168.173d0, 191.837d0, 21.736d0, 15.704d0, 95.528d0, 49.095d0/
  real*8, parameter :: pi=3.1415926535897932d0, d2r=pi/180.

  real*8, intent(IN) :: Longitude_W ! west longitude (degree)
  real*8 EOT, MST
  real*8, intent(OUT) :: LTST  ! (hour)
      
  ! Time offset from J2000 epoch (UT)
  temp1 = (JDUT-2451545.0d0)/36525.d0
  
  ! dcor = JD_TT - JD_UTC  ! (seconds),  TT = terrestrial time
  ! (AM2000, eq. 27)
  ! dcor = (64.184d0 + 95.*temp1 + 35.*temp1**2) ! correction in sec
  ! alternate from https://www.giss.nasa.gov/tools/mars24/help/algorithm.html
  dcor = 64.184 + 59*temp1 - 51.2*temp1**2 - 67.1*temp1**3 - 16.4*temp1**4      ! (seconds)
  ! this correction is small

  JDTT = JDUT + dcor/86400.
  Deltat_J2000 = JDTT - 2451545.0   ! subtract JD(2000)

  !print *,'JDUT=',JDUT
  !print *,'JDTT=',JDTT

      
  ! Mars orbital parameters
  ! Mars mean anomaly (AM2000, eq. 16)
  M = 19.3871d0 + 0.52402073d0*Deltat_J2000 ! (degree)
  !M = 19.3870 + 0.52402075*Deltat_J2000 !  actual (degree)
      
  ! Angle of Fiction Mean Sun (AM2000, eq. 17)
  alpha_FMS = 270.3871d0 + 0.524038496d0*Deltat_J2000 ! (degree)
  !alpha_FMS = 270.3863d0 + 0.52403840d0*Deltat_J2000 ! AM2000 coefficients
  
  ! Perturbers (AM2000, eq. 18)
  PBS = 0.
  do i=1,7
     PBS = PBS+A(i)*cos(2*pi/365.25d0*Deltat_J2000/tau(i)+phi(i)*d2r)
  enddo
  
  ! Equation of Center (AM2000, eqs. 19 and 20)
  M = M*d2r 
  numinusM =(10.691d0 + 3.0d-7*Deltat_J2000)*sin(M)+0.623d0*sin(2*M) &
       &  + 0.050d0*sin(3*M) + 0.005d0*sin(4*M) + 0.0005d0*sin(5*M) + PBS ! (degree)
  
  ! Areocentric solar longitude (AM2000, eq. 19)
  Ls = alpha_FMS + numinusM
  Ls = modulo(Ls,360.d0)*d2r
  
  ! Solar declination (A1997, eq. 5)
  dec = asin(0.42565d0*sin(Ls)) + 0.25d0*d2r*sin(Ls)
  
  ! Heliocentric distance (AM2000, eq. 25, corrected)
  RM = 1.52367934d0 * (1.00436d0 - 0.09309d0*cos(M) - 0.004336d0*cos(2*M) &
       &     - 0.00031d0*cos(3*M) - 0.00003d0*cos(4*M))


  ! Mars Local Time
  ! Equation of Time (AM2000, eq. 20)
  EOT = 2.861d0*sin(2*Ls) - 0.071d0*sin(4*Ls) + 0.002d0*sin(6*Ls) - numinusM  ! (degree)
      
  ! Mean Solar Time (AM2000, eq. 22)
  MST = 24*( (JDTT-2451549.5d0)/1.02749125d0 + 44796.0d0 - 0.00072d0 )
  
  ! Local True Solar Time (AM2000, eq. 23)
  LTST = MST - Longitude_W*24.d0/360.d0 + EOT*24.d0/360.d0
  
  LTST = mod(LTST,24.d0)
  
end subroutine marsclock24
      

