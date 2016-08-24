! Subroutines and functions that are called from cratersQ, fieldofviews, and shadows


elemental function horizontaldistance(i1,j1,i2,j2)
  ! distance between two points; must have same units as height
  use filemanager, only : dx,dy
  implicit none
  real(8) horizontaldistance
  integer, intent(IN) :: i1,j1,i2,j2
  horizontaldistance = sqrt(dx*dx*(i1-i2)**2+dy*dy*(j1-j2)**2)
end function horizontaldistance



elemental function azimuth(i1,j1,i2,j2)
  use filemanager, only : dx,dy
  implicit none
  real(8) azimuth
  integer, intent(IN) :: i1,j1,i2,j2
  azimuth = atan2(dx*(i2-i1),-dy*(j2-j1))   ! this is correct
end function azimuth



real(8) function viewing_angle(i0,j0,i,j,h)
!***********************************************************************
!  function that calculates angle between surface normal at (i,j)
!     and vector pointing to (ii,jj) as in findonehorizon_wsort
!***********************************************************************
  use filemanager, only : NSx,NSy
  use allinterfaces, except_this_one => viewing_angle
  implicit none
  integer, intent(IN) :: i0,j0,i,j
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), parameter :: pi=3.1415926535897931
  real(8) az, s, r, surfaceSlope, azFac, slope_along_az

  !r = sqrt(dx*dx*(i-i0)**2+dy*dy*(j-j0)**2)
  r = horizontaldistance(i,j,i0,j0)
  az = azimuth(i0,j0,i,j)
  s = (h(i,j)-h(i0,j0))/r

  call difftopo1(i0,j0,h,surfaceSlope,azFac)
  slope_along_az=atan(surfaceSlope*cos(azFac-az)) ! an angle

  viewing_angle = pi/2 - atan(s) + slope_along_az 
end function viewing_angle



pure subroutine difftopo1(i,j,h,surfaceSlope,az)
  ! calculate slopes and azimuths of surface elements
  use filemanager, only : NSx,NSy,dx,dy
  implicit none
  integer, intent(IN) :: i,j
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT) :: surfaceSlope,az
  real(8) sx,sy
  
  sx = -1e32; sy=-1e32  ! avoids compiler warning

  if (i>1 .and. i<NSx) then
     sx=(h(i+1,j)-h(i-1,j))/(2.*dx)
  else
    if (i==1) sx=(h(i+1,j)-h(i,j))/dx
    if (i==NSx) sx=(h(i,j)-h(i-1,j))/dx
  endif

  if (j>1 .and. j<NSy) then
     sy=(h(i,j+1)-h(i,j-1))/(2.*dy)
  else
     if (j==1) sy=(h(i,j+1)-h(i,j))/dy
     if (j==NSy) sy=(h(i,j)-h(i,j-1))/dy
  endif

  surfaceSlope=sqrt(sx**2+sy**2)  
  az=atan2(sx,-sy)
end subroutine difftopo1



subroutine compactoutput(unit,value,nr)
  ! output zeros without trailing zeros
  implicit none
  integer, intent(IN) :: unit,nr
  real(8), intent(IN) :: value(nr)
  integer j
  do j=1,nr
     if (value(j)==0.) then
        write(unit,'(1x,f2.0)',advance='no') value(j)
     else
        write(unit,'(1x,f6.4)',advance='no') value(j)
     endif
  enddo
  write(unit,"('')")
end subroutine compactoutput
