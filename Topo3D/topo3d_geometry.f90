! Subroutines and functions that are called from shadow* and also from fieldofview*



elemental function diffangle(a1,a2)
  implicit none
  real(8) diffangle
  real(8), intent(IN) :: a1,a2
  real(8), parameter :: pi=3.1415926535897932
  real(8) x
  x = mod(abs(a1-a2),2*pi)
  if (x<0.) x=x+2*pi
  if (2*pi-x<x) x=2*pi-x
  diffangle = x
end function diffangle



subroutine compactoutput(unit,value,nr)
  ! output zeros without trailing zeros
  implicit none
  integer, intent(IN) :: unit, nr
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



elemental function horizontaldistance1(x1,y1,x2,y2)
  ! distance between two points; must have same units as height
  implicit none
  real(8) horizontaldistance1
  real(8), intent(IN) :: x1,y1,x2,y2
  
  horizontaldistance1 = sqrt((x1-x2)**2+(y1-y2)**2)
end function horizontaldistance1



elemental function azimuth1(x1,y1,x2,y2)
  implicit none
  real(8) azimuth1
  real(8), intent(IN) :: x1,y1,x2,y2
  
  azimuth1 = atan2(x2-x1, -(y2-y1))  ! this is correct
end function azimuth1



subroutine downsample(NSx,NSy,h,hhalf)
  ! downsample 2D array to half resolution
  implicit none
  integer, intent(IN) :: NSx, NSy
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT) :: hhalf(NSx/2,NSy/2) ! new dimensions
  integer :: x(NSx/2), y(NSy/2)
  integer i2, j2, ii, jj, k, r
  real(8) psum, w(NSx/2,NSy/2)
  integer, parameter :: ei(5) = (/0, 1, 0, -1, 0 /) ! counter-clockwise
  integer, parameter :: ej(5) = (/0, 0, 1, 0, -1 /)
  integer, parameter :: nan = -32000 

  do i2=1,NSx/2
     x(i2) = 2*i2
     
     do j2=1,NSy/2
        y(j2) = 2*j2

        if (x(i2)>NSx .or. y(j2)>NSy) cycle

        ! weighted average
        psum = 0.; w(i2,j2)=0.
        do k=1,5
           ii = x(i2) + ei(k)
           jj = y(j2) + ej(k)
           if (ii>NSx .or. ii<1) cycle
           if (jj>NSy .or. jj<1) cycle
           if (h(ii,jj)<=nan) cycle
           r = 1+3*(ei(k)**2+ej(k)**2)  ! should be 1 or 4
           psum = psum + h(ii,jj)/r           
           w(i2,j2) = w(i2,j2) + 1./r
        enddo
        if (w(i2,j2)>=1.5) then    ! 1/1 + 4/4 = 2
           hhalf(i2,j2) = psum/w(i2,j2)
        else ! only two or fewer valid neighbors
           hhalf(i2,j2) = nan
        endif

     enddo
  enddo
  !do j2=1,NSy/2
  !   write(6,'(9999(f3.1,1x))') w(:,j2)
  !enddo
  !write(6,*)
end subroutine downsample
