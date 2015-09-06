! THE geographic grid
! centered on equator, poles are extra grid points

module grid
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8) dlon, dlat
  integer nlon, nlat2
  !parameter (dlon=7.5, nlon=48) 
  parameter (dlon=1, nlon=360) 

  !parameter (dlat=0.5, nlat2=180)
  parameter (dlat=2., nlat2=45)

  integer, parameter :: nlat=2*nlat2-1
  integer, parameter :: VECLEN = nlon*nlat+2
end module grid


pure subroutine lonlatgrid(longitude,latitude)
  ! set up lon-lat grid
  use grid
  implicit none
  real(8), intent(OUT) :: longitude(nlon),latitude(nlat)
  integer i,j
  do i=1,nlon
     longitude(i)=(i-1)*dlon
  enddo
  do j=1,nlat
     latitude(j)=(nlat2-j)*dlat
  enddo
end subroutine lonlatgrid


integer function inbox(r)
  ! determine what lon/lat box coordinate r is in
  ! should be shifted by 1/2
  use grid
  implicit none
  real(8), intent(IN) :: r(2)
  integer kx, ky
  kx = ceiling(r(1)/dlon+0.5)
  ky = -ceiling(r(2)/dlat-0.5)+nlat2 
  !if (kx<1) kx=kx+nlon
  if (kx>nlon) kx=kx-nlon
  inbox = 1 + kx + (ky-1)*nlon
  if (ky>nlat) inbox=veclen
  if (ky<1) inbox=1
  if (inbox<1 .or. inbox>veclen .or. kx<1 .or. kx>nlon) then
     print *,'inbox: fatal return value'
     print *,'r=',r
     print *,kx,ky,inbox
     !stop
  endif
end function inbox


pure subroutine k2lonlat(k,lon,lat)
  ! inverse of  k = 1 + i + (j-1)*nlon
  ! linear enumeration of grid cells, but first element is north pole and last element is south pole
  use grid
  implicit none
  integer, intent(IN) :: k
  real(8), intent(OUT) :: lon  ! (degree)
  real(8), intent(OUT) :: lat  ! (radian)
  integer i,j
  !if (k<1 .or. k>veclen) stop 'k2lonlat: Index k is out of bound'
  lon = 0.
  if (k>1 .and. k<veclen) then  
     i = mod(k-2,nlon)+1   ! longitude index
     lon = (i-1)*dlon  ! =longitude(i)
     j = (k-2)/nlon + 1   ! latitude index
     !if (j<1 .or. j>nlat) stop 'k2lonlat: Index j is out of bound'
     lat = (nlat2-j)*dlat*d2r   ! =latitude(j)*d2r
  endif
  if (k==1) lat=pi/2.
  if (k==veclen) lat=-pi/2.
end subroutine k2lonlat


subroutine areas(dA)
  ! areas of surface elements
  use grid, only: pi,d2r,dlat,dlon,veclen
  implicit none
  real(8), intent(OUT) :: dA(veclen)
  integer k
  real(8) lon,lat

  dA(1) = (dlat*d2r)**2*pi 
  do k=2,veclen-1
     call k2lonlat(k,lon,lat)
     dA(k) = dlat*dlon*d2r**2*cos(lat)
  enddo
  dA(veclen) = dA(1)
end subroutine areas


subroutine writeglobe(unit,Tsurf)
  use grid, only: nlon, nlat, veclen
  implicit none
  integer, intent(IN) :: unit
  real(8), intent(IN) :: Tsurf(*)
  integer i,j,k
  real(8) longitude(nlon), latitude(nlat)

  ! set up lon-lat grid
  call lonlatgrid(longitude,latitude)

  write(unit,100) 0.,90.,Tsurf(1)
  do j=1,nlat
     do i=1,nlon
        k = 1 + i + (j-1)*nlon
        write(unit,100) longitude(i),latitude(j),Tsurf(k)
     enddo
  enddo
  write(unit,100) 0.,-90.,Tsurf(veclen)
100 format (f5.1,1x,f6.2,1x,f7.3)
end subroutine writeglobe
