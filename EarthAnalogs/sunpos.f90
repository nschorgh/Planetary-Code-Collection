module dateformat
  type cTime
     integer iYear, iMonth, iDay
     real(8) dHours, dMinutes, dSeconds
  end type cTime
end module dateformat


subroutine sunpos(udtTime, dLongitude, dLatitude, dZenithAngle, dAzimuth, R)
  ! Adapted from
  ! Blanco-Muriel et al. (2001), Computing the solar vector, Solar Energy 70, 431-441
  use dateformat
  implicit none

  type(cTime), intent(IN) :: udtTime
  real(8), intent(IN) :: dLongitude, dLatitude  ! cLocation
  real(8), intent(OUT) :: dZenithAngle, dAzimuth, R  ! cSunCoordinates

  real(8), parameter ::  pi =   3.14159265358979323846
  real(8), parameter ::  twopi = 2*pi, rad = pi/180.
  real(8), parameter ::  dEarthMeanRadius =  6371.01    ! In km
  real(8), parameter ::  dAstronomicalUnit = 149597890  ! In km

  ! Main variables
  real(8) dElapsedJulianDays, dDecimalHours
  real(8) dEclipticLongitude, dEclipticObliquity
  real(8) dRightAscension, dDeclination

  ! Auxiliary variables
  real(8) dY, dX
  real(8) dJulianDate
  integer(8) liAux1, liAux2
  real(8) dMeanLongitude, dMeanAnomaly
  real(8) dOmega, dSin_EclipticLongitude
  real(8) dGreenwichMeanSiderealTime, dLocalMeanSiderealTime
  real(8) dLatitudeInRadians, dHourAngle
  real(8) dCos_Latitude, dSin_Latitude
  real(8) dCos_HourAngle
  real(8) dParallax

  real(8), external :: distancefromsun

  ! Calculate difference in days between the current Julian Day 
  ! and JD 2451545.0, which is noon 1 January 2000 Universal Time

  ! Calculate time of the day in UT decimal hours
  dDecimalHours = udtTime%dHours + (udtTime%dMinutes + udtTime%dSeconds / 60.0 ) / 60.0
  ! Calculate current Julian Day
  liAux1 =(udtTime%iMonth-14)/12
  liAux2=(1461*(udtTime%iYear + 4800 + liAux1))/4 + (367*(udtTime%iMonth &
       & - 2-12*liAux1))/12- (3*((udtTime%iYear + 4900 & 
       & + liAux1)/100))/4+udtTime%iDay-32075
  dJulianDate=real(liAux2,8)-0.5+dDecimalHours/24.0
  ! Calculate difference between current Julian Day and JD 2451545.0 
  dElapsedJulianDays = dJulianDate-2451545.0

  ! Calculate ecliptic coordinates (ecliptic longitude and obliquity of the 
  ! ecliptic in radians but without limiting the angle to be less than 2*Pi 
  ! (i.e., the result may be greater than 2*Pi)

  dOmega=2.1429-0.0010394594*dElapsedJulianDays
  dMeanLongitude = 4.8950630+ 0.017202791698*dElapsedJulianDays ! Radians
  dMeanAnomaly = 6.2400600+ 0.0172019699*dElapsedJulianDays
  dEclipticLongitude = dMeanLongitude + 0.03341607*sin( dMeanAnomaly ) &
       & + 0.00034894*sin( 2*dMeanAnomaly )-0.0001134 &
       & -0.0000203*sin(dOmega)
  dEclipticObliquity = 0.4090928 - 6.2140e-9*dElapsedJulianDays  +0.0000396*cos(dOmega)

  ! Calculate celestial coordinates ( right ascension and declination ) in radians 
  ! but without limiting the angle to be less than 2*Pi (i.e., the result may be 
  ! greater than 2*Pi)
  dSin_EclipticLongitude= sin( dEclipticLongitude )
  dY = cos( dEclipticObliquity ) * dSin_EclipticLongitude
  dX = cos( dEclipticLongitude )
  dRightAscension = atan2( dY,dX )
  if( dRightAscension < 0.0 ) dRightAscension = dRightAscension + twopi
  dDeclination = asin( sin( dEclipticObliquity )*dSin_EclipticLongitude )

  ! Calculate local coordinates ( azimuth and zenith angle ) in degrees
  dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283*dElapsedJulianDays  + dDecimalHours
  dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime*15   + dLongitude)*rad
  dHourAngle = dLocalMeanSiderealTime - dRightAscension
  dLatitudeInRadians = dLatitude*rad
  dCos_Latitude = cos( dLatitudeInRadians )
  dSin_Latitude = sin( dLatitudeInRadians )
  dCos_HourAngle= cos( dHourAngle )
  dZenithAngle = (acos( dCos_Latitude*dCos_HourAngle*cos(dDeclination) + sin( dDeclination )*dSin_Latitude))
  dY = -sin( dHourAngle )
  dX = tan( dDeclination )*dCos_Latitude - dSin_Latitude*dCos_HourAngle
  dAzimuth = atan2( dY, dX )
  if ( dAzimuth < 0.0 ) dAzimuth = dAzimuth + twopi
  dAzimuth = dAzimuth/rad
  ! Parallax Correction
  dParallax=(dEarthMeanRadius/dAstronomicalUnit)*sin(dZenithAngle)
  dZenithAngle=(dZenithAngle + dParallax)/rad

  r = distancefromsun(dMeanAnomaly)
end subroutine sunpos


function distancefromsun(dMeanAnomaly)
  implicit none
  real(8), intent(IN) :: dMeanAnomaly  ! (radians)
  real(8) distancefromsun

  real(8), parameter :: a=1     ! semimajor axis (AU)
  real(8), parameter :: ecc=0.0167  ! orbital eccentricity
  integer j
  real*8 E,Eold
  
  ! E = eccentric anomaly 
  ! solve M = E - ecc*sin(E) by Newton Method
  E = dMeanAnomaly
  do j=1,10
     Eold = E
     E = E - (E - ecc*sin(E) - dMeanAnomaly)/(1.-ecc*cos(E))
     if (abs(E-Eold)<1.e-8) exit
  enddo
  
  ! heliocentric distance (AU)
  distancefromsun = a*(1-ecc*cos(E))
END function distancefromsun


pure subroutine addtime(udtTime)
  ! normalize overflow times and dates
  use dateformat
  implicit none
  type(cTime), intent(INOUT) :: udtTime
  integer shft

  do while (udtTime%dSeconds >= 60)
     udtTime%dSeconds = udtTime%dSeconds - 60
     udtTime%dMinutes = udtTime%dMinutes + 1
  end do
  do while (udtTime%dMinutes >= 60)
     udtTime%dMinutes = udtTime%dMinutes - 60
     udtTime%dHours = udtTime%dHours + 1
  end do
  do while (udtTime%dHours >= 24)
     udtTime%dHours = udtTime%dHours - 24
     udtTime%iDay = udtTime%iDay + 1
  end do

  select case (udtTime%iMonth)
  case(1,3,5,7,8,10,12)
     shft = 31
  case(4,6,9,11)
     shft = 30
  case(2)
     shft = 28
     if (mod(udtTime%iYear,4)==0) shft = 29   ! leap year
     if (mod(udtTime%iYear,100)==0) shft = 28 ! no leap year
     if (mod(udtTime%iYear,400)==0) shft = 29 ! leap year
  case default
     shft = -9999
     !stop 'error in addtime'
  end select

  if (udtTime%iDay>shft) then
     udtTime%iDay = udtTime%iDay - shft
     udtTime%iMonth = udtTime%iMonth + 1
  endif
  if (udtTime%iMonth > 12) then
     udtTime%iMonth = udtTime%iMonth - 12
     udtTime%iYear = udtTime%iYear + 1
  endif
end subroutine addtime


pure subroutine subtime(udtTime)
  ! normalize underflow times and dates
  use dateformat
  implicit none
  type(cTime), intent(INOUT) :: udtTime
  integer shft

  do while (udtTime%dSeconds < 0) 
     udtTime%dSeconds = udtTime%dSeconds + 60
     udtTime%dMinutes = udtTime%dMinutes - 1
  end do
  do while (udtTime%dMinutes < 0)
     udtTime%dMinutes = udtTime%dMinutes + 60
     udtTime%dHours = udtTime%dHours - 1
  end do
  do while (udtTime%dHours < 0) 
     udtTime%dHours = udtTime%dHours + 24
     udtTime%iDay = udtTime%iDay - 1
  end do

  select case (udtTime%iMonth-1)
  case(1,3,5,7,8,10,0)
     shft = 31
  case(4,6,9,11)
     shft = 30
  case(2)
     shft = 28
     if (mod(udtTime%iYear,4)==0) shft = 29   ! leap year
     if (mod(udtTime%iYear,100)==0) shft = 28 ! no leap year
     if (mod(udtTime%iYear,400)==0) shft = 29 ! leap year
  case default
     shft = -9999
     !stop 'error in subtime'
  end select

  if (udtTime%iDay <= 0) then 
     udtTime%iDay = udtTime%iDay + shft
     udtTime%iMonth = udtTime%iMonth - 1
  endif
  if (udtTime%iMonth <= 0) then
     udtTime%iMonth = udtTime%iMonth + 12
     udtTime%iYear = udtTime%iYear - 1
  endif
end subroutine subtime


subroutine timezone(time,zone)
  use dateformat
  implicit none
  type(cTime), intent(INOUT) :: time
  integer, intent(IN) :: zone

  time%dHours = time%dHours + zone
  if (zone>0) call addtime(time)
  if (zone<0) call subtime(time)
end subroutine timezone

