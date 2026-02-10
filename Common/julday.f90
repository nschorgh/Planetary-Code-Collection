function julday(mm,dd,yyyy)
  ! calculate Julian date number from calendar date
  ! Equations from D. A. Hatcher (1984). Q. J. R. Astr. Soc. 25, 53-55
  ! also see Numerical Recipes (1992)
  implicit none
  integer julday
  integer, intent(IN) :: mm,dd,yyyy
  integer, parameter :: IGREG=15+31*(10+12*1582) ! switch of calendar
  integer yp,mp,g,JDN

  yp = yyyy-int((12-mm)/10)
  mp = modulo(mm-3,12)  ! start year on March 1
  JDN = floor(365.25*(yp+4712)) + floor(30.6001*mp+0.5) + dd + 59
  if ( dd+31*(mm+12*yyyy) < IGREG ) then
     julday = JDN
  else ! after switch to Gregorian
     g = floor( ((yp/100)+49) * 0.75 ) - 38
     !g = floor( ((yp/100)+1) * 0.75) - 2  ! as long as yp>0
     julday = JDN-g
  endif

end function julday

