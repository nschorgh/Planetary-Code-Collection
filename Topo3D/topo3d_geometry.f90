! Subroutines and functions that are called from shadow* and also from fieldofview*



elemental function diffangle(a1,a2)
  implicit none
  real(8) diffangle
  real(8), intent(IN) :: a1,a2
  real(8), parameter :: pi=3.1415926535897932
  real(8) x
  x=mod(abs(a1-a2),2*pi)
  if (x<0.) x=x+2*pi
  if (2*pi-x<x) x=2*pi-x
  diffangle = x
end function diffangle



subroutine compactoutput(unit,value,nr)
  ! output zeros without trailing zeros
  implicit none
  integer, intent(IN) :: unit,nr
  real(8), intent(IN) :: value(nr)
  integer j
  do j=1,nr
     if (value(j)==0.) then
        write(unit,'(1x,f2.0)',advance='no') value(j)
     elseif (value(j)>=9.99995d0) then ! format overflow
        write(unit,'(1x,f6.4)',advance='no') 9.9999
     else
        write(unit,'(1x,f6.4)',advance='no') value(j)
     endif
  enddo
  write(unit,"('')")
end subroutine compactoutput



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



elemental function horizontaldistance1(x1,y1,x2,y2)
  ! distance between two points; must have same units as height
  ! as in topo3d_common.f90, but distance based
  implicit none
  real(8) horizontaldistance1
  real(8), intent(IN) :: x1,y1,x2,y2
  
  horizontaldistance1 = sqrt((x1-x2)**2+(y1-y2)**2)
end function horizontaldistance1



elemental function azimuth1(x1,y1,x2,y2)
  ! as in topo3d_common.f90, but distance based
  implicit none
  real(8) azimuth1
  real(8), intent(IN) :: x1,y1,x2,y2
  
  azimuth1 = atan2(x2-x1,-(y2-y1)) 
end function azimuth1

