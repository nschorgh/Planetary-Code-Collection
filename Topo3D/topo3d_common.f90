! Subroutines and functions that are called from cratersQ_* and also from fieldofviews or shadows


subroutine readdem(h)
  ! read DEM
  ! (1,1) = northwest corner
  ! (NSx,1) = northeast corner
  ! (1,NSy) = southwest corner
  ! (NSx,NSy) = southeast corner
  use filemanager, only : NSx, NSy, hfn
  implicit none
  real(8), intent(OUT) :: h(NSx,NSy)
  integer j, ierr, ios

  open(unit=20,file=hfn,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'readdem: input file not found'
  do j=1,NSy
     read(20,*,iostat=ios) h(:,j)
     if (ios /= 0) stop 'readdem: unexpected number of rows'
  enddo
  close(20)
end subroutine readdem



pure function area_spherical_quadrangle(phi,theta)
! area of two triangles on sphere
  use allinterfaces, only : area_spherical_triangle
  implicit none
  real(8), intent(IN) :: phi(4), theta(4)
  real(8) E1, E2, area_spherical_quadrangle
  
  E1 = area_spherical_triangle((/ phi(1), phi(2), phi(3) /), &
       & (/ theta(1), theta(2), theta(3) /))
  E2 = area_spherical_triangle((/ phi(3), phi(4), phi(1) /), &
       & (/ theta(3), theta(4), theta(1) /))

  area_spherical_quadrangle = E1 + E2
end function area_spherical_quadrangle



pure function area_spherical_triangle(phi,theta)
  use allinterfaces, only : distanceonsphere
  implicit none
  real(8), intent(IN) :: phi(3), theta(3)
  real(8) a,b,c
  real(8) s, E, area_spherical_triangle, buf

  a=distanceonsphere(phi(1),theta(1),phi(2),theta(2))
  b=distanceonsphere(phi(2),theta(2),phi(3),theta(3))
  c=distanceonsphere(phi(3),theta(3),phi(1),theta(1))

  ! spherical excess
  s = (a+b+c)/2.
  buf=tan(s/2)*tan((s-a)/2.)*tan((s-b)/2.)*tan((s-c)/2.)
  !buf=tan((a+b+c)/4.)*tan((-a+b+c)/4.)*tan((a-b+c)/4.)*tan((a+b-c)/4.)
  if (buf<0.) buf=0.   ! roundoff
  E=4*atan(sqrt(buf))

  area_spherical_triangle = E
end function area_spherical_triangle



elemental function distanceonsphere(phi1,theta1,phi2,theta2)
! spherical distance between two points in radians
  implicit none
  real(8), intent(IN) :: phi1,phi2,theta1,theta2  ! [radians]
  real(8) distanceonsphere, buf, lat1, lat2
  real(8), parameter :: pi=3.1415926535897932
  lat1=pi/2-theta1; lat2=pi/2-theta2

  ! buf = square of half of cord length distance
  buf = sin((lat1-lat2)/2.)**2+cos(lat1)*cos(lat2)*sin((phi1-phi2)/2.)**2
  distanceonsphere = 2.*asin(sqrt(buf))
  !distanceonsphere = 2.*atan2(sqrt(buf),sqrt(1-buf))
end function distanceonsphere



subroutine slicer(NSx,ilower,iupper,extc)
  ! splits domain for parallel processing
  ! two input arguments required
  implicit none
  integer, intent(IN) :: NSx
  integer, intent(OUT) :: ilower, iupper
  character(4), intent(OUT) :: extc
  integer SLICE  ! number of slices (jobs) the domain is divided into
  integer nr, slicewidth, Mx1, Mx2  ! 1 <= nr <= SLICE
  
  call getarg(1,extc)
  read(extc,'(i4)') SLICE  ! string->integer
  call getarg(2,extc)
  read(extc,'(i4)') nr  ! string->integer
  if (SLICE>NSx-2) stop 'not that many slices available'
  if (nr<1 .or. nr>SLICE) stop 'slice id outside of range'
  slicewidth = ceiling((NSx-2)/real(SLICE))
  if (slicewidth<1) stop 'no slice width'
  print *,'Number of slices=',SLICE,'slice width=',slicewidth
  Mx1 = 2+slicewidth*(nr-1)
  Mx2 = min(2+slicewidth*nr,NSx)-1
  print *,'Working on slice Mx1=',Mx1,'Mx2=',Mx2
  if (Mx1>Mx2) stop 'FYI: out of slices - nothing left to do'
  if (Mx1<=1 .or. Mx1>=NSx .or. Mx2<=1 .or. Mx2>=NSx) stop 'argument is outside of domain'
  ilower=Mx1; iupper=Mx2

  extc = trim(extc) ! strips trailing spaces
end subroutine slicer
     
