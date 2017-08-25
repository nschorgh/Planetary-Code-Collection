pure subroutine findonehorizon(h,i0,j0,azRay,smax)
! finds horizon for one azimuth ray
  use filemanager, only : NSx,NSy,RMAX
  use allinterfaces, only : horizontaldistance, azimuth, diffangle
  implicit none
  integer, intent(IN) :: i0,j0
  real(8), intent(IN) :: h(NSx,NSy),azRay
  real(8), intent(OUT) :: smax
  integer i,j,k,in,jn
  real(8) az,az_neighbor,t,r,r_neighbor,rcut,s,hcut,d1,d2,d3,maxsep
  integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
  integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)

  smax=0.
  maxsep = horizontaldistance(1,1,2,2)
  do i=2,NSx-1
     if (horizontaldistance(i,1,i0,1)>RMAX) cycle  ! saves computations
     if (sin(azRay)*(i-i0) < 0.) cycle  ! saves computations
     do j=2,NSy-1
        if (i==i0 .and. j==j0) cycle
        r = horizontaldistance(i,j,i0,j0)
        if (r>RMAX) cycle  ! saves computations
        az = azimuth(i0,j0,i,j)
        if (cos(azRay)*(j-j0) > 0.) cycle  ! saves computations

        d1=diffangle(az,azRay)
        if (r*d1*3.1416 > maxsep*3) cycle   ! saves computations
        if (d1==0.) then  ! grid point lies on ray
           s = (h(i,j)-h(i0,j0))/r
           if (s>smax) smax=s
           cycle
        endif

        do k=1,8
           in = i+ex(k)
           jn = j+ey(k)
           if (in==i0 .and. jn==j0) cycle
           !if ((in==i0+1.and.jn==j0) .or. (in==i0-1.and.jn==j0) .or. &
           !     & (in==i0.and.jn==j0+1) .or. (in==i0.and.jn==j0-1)) cycle
           az_neighbor = azimuth(i0,j0,in,jn)

           d2=diffangle(az_neighbor,azRay)
           d3=diffangle(az,az_neighbor)

           if (d1+d2<=d3+1.d-5) then  
              if (d1>0.5*d3 .and. d3>1.d-6) cycle ! try this
              r_neighbor = horizontaldistance(in,jn,i0,j0)
              ! edge between h1,i0,j0 and h2,in,jn
              if (d3>1.d-6) then
                 t = d1/d3  ! approximation
              else
                 t = 0.5  ! dirty fix
              endif
              hcut = h(i,j)*(1-t)+h(in,jn)*t
              rcut = r*(1-t)+r_neighbor*t
              s = (hcut-h(i0,j0))/rcut
              if (s>smax) smax=s
           endif
        end do
     enddo
  enddo
end subroutine findonehorizon



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



subroutine findonehorizon_wsort(h,i0,j0,azRay,smax,visibility)
! finds horizon and determines visibility for one azimuth ray
  use filemanager, only : NSx,NSy,RMAX
  use allinterfaces, only : horizontaldistance, azimuth, diffangle, hpsort
  implicit none
  integer, intent(IN) :: i0,j0
  real(8), intent(IN) :: h(NSx,NSy),azRay
  real(8), intent(OUT) :: smax
  logical, intent(INOUT) :: visibility(NSx,NSy)
  integer i,j,k,in,jn,cc
  integer, parameter :: CCMAX = 5*(NSx+NSy) ! max # of elements on azimuth ray
  real(8) az,az_neighbor,t,r,r_neighbor,hcut,d1,d2,d3
  real(8) s, smaxlocal, surfaceSlope, azFac, slope_along_az
  real(8), dimension(CCMAX) :: rcut, slocal
  integer, dimension(CCMAX) :: celli, cellj, arr
  integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
  integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)

  cc=0
  smax=0.

  call difftopo1(i0,j0,h,surfaceSlope,azFac)
  slope_along_az=atan(surfaceSlope*cos(azFac-azRay)) ! an angle

  do i=2,NSx-1
     !if (dx*abs(i-i0)>RMAX) cycle  ! to save computational cost
     if (horizontaldistance(i,1,i0,1)>RMAX) cycle  ! to save computational cost
     do j=2,NSy-1
        if (i==i0 .and. j==j0) cycle
        !r = sqrt(dx*dx*(i-i0)**2+dy*dy*(j-j0)**2)
        r = horizontaldistance(i,j,i0,j0)
        if (r>RMAX) cycle  ! to save computational cost
        az = azimuth(i0,j0,i,j)

        d1=diffangle(az,azRay)
        if (d1==0.) then  ! grid point lies on ray
           cc = cc+1
           if (cc>CCMAX) stop 'findonehorizon_wgeom: not enough memory allocated'
           s = (h(i,j)-h(i0,j0))/r
           if (s>smax) smax=s
           
           rcut(cc) = r
           slocal(cc) = tan(atan(s)-slope_along_az)
           celli(cc)=i; cellj(cc)=j
           cycle
        endif
        
        do k=1,8
           in = i+ex(k)
           jn = j+ey(k)
           if (in==i0 .and. jn==j0) cycle
           if ((in==i0+1.and.jn==j0) .or. (in==i0-1.and.jn==j0) .or. &
                & (in==i0.and.jn==j0+1) .or. (in==i0.and.jn==j0-1)) cycle
           az_neighbor = azimuth(i0,j0,in,jn)

           d2=diffangle(az_neighbor,azRay)
           d3=diffangle(az,az_neighbor)
           ! do something about d3<=1e-6

           if (d1+d2<=d3+1e-5) then 
              if (d1>0.5*d3 .and. d3>1.e-6) cycle
              cc=cc+1
              if (cc>CCMAX) stop 'findonehorizon_wgeom: not enough memory allocated'

              !r_neighbor = sqrt(dx*dx*(in-i0)**2+dy*dy*(jn-j0)**2)
              r_neighbor = horizontaldistance(in,jn,i0,j0)
              ! edge between h1,i0,j0 and h2,in,jn
              if (d3>1.e-6) then
                 t = d1/d3  ! approximation
              else
                 t = 0.5  ! dirty fix
              endif
              hcut = h(i,j)*(1-t)+h(in,jn)*t
              rcut(cc) = r*(1-t)+r_neighbor*t  ! could be improved
              s = (hcut-h(i0,j0))/rcut(cc)
              if (s>smax) smax=s

              slocal(cc) = tan(atan(s)-slope_along_az)
              celli(cc)=i; cellj(cc)=j
           endif
        end do
     enddo
  enddo

  ! sort by distance from (i0,j0)
  call hpsort(cc,rcut,arr)
  
  smaxlocal=0.
  do i=1,cc
     j=arr(i)
     if (slocal(j)>smaxlocal) then
        smaxlocal=slocal(j)
        visibility(celli(j),cellj(j))=.true.
     endif
  enddo
end subroutine findonehorizon_wsort



subroutine findallhorizon(h,i0,j0,naz,smax)
! finds horizon for all azimuth rays
  use filemanager, only : NSx,NSy,RMAX
  use allinterfaces, only : horizontaldistance, azimuth, diffangle
  implicit none
  integer, intent(IN) :: i0,j0,naz
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT) :: smax(naz)
  integer i,j,k,in,jn,ak,ak1,ak2,buf,akak
  real(8) az,az_neighbor,t,r,r_neighbor,rcut,s,hcut,d1,d2,d3
  integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
  integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)
  real(8) f,azRay(naz)
  real(8), parameter :: pi=3.1415926535897931

  f = naz/(2*pi)
  
  ! azimuth to integer mapping
  do ak=1,naz
     azRay(ak) = (ak-1)/f
     ! inverse mapping:  ak = azRay*f+1
  enddo
  
  smax(:)=0.
  do i=2,NSx-1
     if (horizontaldistance(i,1,i0,1)>RMAX) cycle  ! saves computations
     do j=2,NSy-1
        if (i==i0 .and. j==j0) cycle
        r = horizontaldistance(i,j,i0,j0)
        if (r>RMAX) cycle  ! saves computations
        az = azimuth(i0,j0,i,j)  ! az should go from -pi ... pi

        if (floor(az*f)==ceiling(az*f)) then  ! grid point lies on ray
           ak=nint(az*f)+1
           s = (h(i,j)-h(i0,j0))/r
           if (s>smax(ak)) smax(ak)=s
           cycle
        endif

        do k=1,8
           in = i+ex(k)
           jn = j+ey(k)
           if (in==i0 .and. jn==j0) cycle
           az_neighbor = azimuth(i0,j0,in,jn)

           if (az >= az_neighbor) then 
              ak1 = floor(az*f)+1
              ak2 = ceiling(az_neighbor*f)+1
           else
              ak1 = ceiling(az*f)+1
              ak2 = floor(az_neighbor*f)+1
           endif
           if (ak1==naz/2 .and. ak2<0) ak2 = ak2+naz
           if (ak2==naz/2 .and. ak1<0) ak1 = ak1+naz
           if (ak2<ak1) then ! swap
              buf=ak1; ak1=ak2; ak2=buf;
           endif
           if (ak1>naz .or. ak2>naz) stop 'findallhorizon: index out of bound'

           d3=diffangle(az,az_neighbor)
           do akak=ak1,ak2
              ak = akak; if (ak<=0) ak = ak+naz
              d1=diffangle(az,azRay(ak))
              d2=diffangle(az_neighbor,azRay(ak))

              if (d1+d2<=d3+1.d-5) then
                 if (d1>=1.0*d3 .and. d3>1.d-6) cycle ! DIFFERENT
                 !r_neighbor = sqrt(dx*dx*(in-i0)**2+dy*dy*(jn-j0)**2)
                 r_neighbor = horizontaldistance(in,jn,i0,j0)
                 ! edge between h1,i0,j0 and h2,in,jn
                 if (d3>1.d-6) then
                    t = d1/d3  ! approximation
                 else
                    t = 0.5  ! dirty fix
                 endif
                 hcut = h(i,j)*(1-t)+h(in,jn)*t
                 rcut = r*(1-t)+r_neighbor*t
                 s = (hcut-h(i0,j0))/rcut
                 if (s>smax(ak)) smax(ak)=s
              endif

           end do  ! end ak loop
        end do ! end of k loop
        
     enddo ! end of j loop
  enddo ! end of i loop
end subroutine findallhorizon



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

