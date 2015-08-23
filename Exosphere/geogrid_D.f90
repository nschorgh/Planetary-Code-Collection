! THE geographic grid
! matches Diviner temperature input maps
! no grid points on poles and equator

module grid
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8) dlon, dlat
  integer nlon, nlat2
  parameter (dlon=1, nlon=360) 
  parameter (dlat=0.5, nlat2=180)
  integer, parameter :: nlat=2*nlat2
  integer, parameter :: VECLEN = nlon*nlat
end module grid


pure subroutine lonlatgrid(longitude,latitude)
  ! set up lon-lat grid
  use grid
  implicit none
  real(8), intent(OUT) :: longitude(nlon),latitude(nlat)
  integer i,j
  do i=1,nlon
     longitude(i)=(i-0.5)*dlon
  enddo
  do j=1,nlat
     latitude(j)=90.-(j-0.5)*dlat
  enddo
end subroutine lonlatgrid


integer function inbox(r) 
  ! determine what lon/lat box coordinate r is in
  use grid
  implicit none
  real(8), intent(IN) :: r(2)
  integer kx, ky
  kx = nint(r(1)/dlon+0.5)
  ky = nint((90.-r(2))/dlat+0.5)
  if (r(2)==-90.) ky=nlat  ! roundoff issue
  if (kx>nlon) kx=kx-nlon
  inbox = kx + (ky-1)*nlon
  if (ky<1 .or. ky>nlat) then
     print *,'inbox: Index ky is out of bound',ky,r(2)
     !stop
  endif
  if (kx<1 .or. kx>nlon) then
     print *,'inbox: Index kx is out of bound',kx,r(1)
     !stop
  endif
  if (inbox<1 .or. inbox>veclen) then
     print *,'inbox: fatal return value',r
     print *,kx,ky,inbox
     !stop
  endif
end function inbox


pure subroutine k2lonlat(k,lon,lat)
  ! inverse of  k = i + (j-1)*nlon
  ! linear enumeration of grid cells
  use grid
  implicit none
  integer, intent(IN) :: k
  real(8), intent(OUT) :: lon  ! (degree)
  real(8), intent(OUT) :: lat  ! (radian)
  integer i,j
  !if (k<1 .or. k>veclen) stop 'k2lonlat: Index k is out of bound'

  i = mod(k,nlon)   ! longitude index
  if (i==0) i=nlon
  lon = (i-0.5)*dlon  ! =longitude(i) in lonlatgrid

  j = (k-1)/nlon + 1   ! latitude index
  !if (j<1 .or. j>nlat) stop 'k2lonlat: Index j is out of bound'
  lat = 90.-(j-0.5)*dlat  ! =latitude(j) in lonlatgrid
  lat = lat*d2r 
end subroutine k2lonlat


subroutine areas(dA)
  ! areas of surface elements
  use grid, only : pi,d2r,dlat,dlon,veclen
  implicit none
  real(8), intent(OUT) :: dA(veclen)
  integer k
  real(8) lon,lat

  do k=1,veclen
     call k2lonlat(k,lon,lat)
     dA(k) = dlat*dlon*d2r**2*cos(lat)
  enddo
end subroutine areas
