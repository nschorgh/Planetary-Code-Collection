MODULE findvisibletopo
!***********************************************************************
! identify hidden versus visible topography
!***********************************************************************
  use filemanager, only : NSx, NSy, dx, dy

  integer, parameter, public :: naz = 360

  integer, parameter, private :: CCMAX = 6*(NSx+NSy) 
  ! CCMAX is the max # of elements in one azimuth ray
  integer, private :: cc(naz)
  real(8), dimension(naz,CCMAX), private :: rcut, slocal 
  integer, dimension(naz,CCMAX), private :: celli, cellj 

  integer, private :: ak
  real(8), parameter, private :: pi = 3.1415926535897932
  real(8), parameter, private :: f = naz/(2*pi)
  real(8), parameter, public :: azRay(naz) = (/ ( (ak-1)/f, ak=1,naz) /)
  
contains

  subroutine findallhorizon_wsort_v3(h,i0,j0,smax,visibility)
    ! finds horizons and determines visibility along all azimuth rays
    use filemanager, only : NSx, NSy, RMAX, dx, dy
    use allinterfaces, only : horizontaldistance1
    implicit none
    real(8), intent(IN) :: h(NSx,NSy)
    integer, intent(IN) :: i0, j0
    real(8), intent(OUT) :: smax(naz)
    logical, intent(OUT) :: visibility(NSx,NSy)
    integer i, j, ak
    real(8) surfaceSlope, azFac, resolved
    real(8) smaxlocal, x0, y0, h00
    integer, dimension(CCMAX) :: arr

    x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)
    resolved = naz*min(dx,dy)/(2*pi)
    !if ( resolved**2 < (NSx*dx)**2 + (NSy*dy)**2 ) then
    !   print *,'findvisibletopo: visibility may be underresolved. increase naz.'
    !endif
    
    cc(:)=0
    smax(:)=0.

    call difftopo1(NSx,NSy,i0,j0,h,dx,dy,surfaceSlope,azFac)
    visibility(:,:) = .false.
    
    do i=2,NSx-1
       if (horizontaldistance1(i*dx,1*dy,x0,1*dy)>RMAX) cycle 
       do j=2,NSy-1
          if (i==i0 .and. j==j0) cycle
          if (horizontaldistance1(i*dx,j*dy,x0,y0)>RMAX) cycle
          
          call horizon_core_wsort(x0,y0,h00,smax,surfaceSlope,azFac,i,j,h)
       end do  ! end of j loop
    end do  ! end of i loop
    
    do ak = 1,naz
       ! sort by distance from (i0,j0)
       call hpsort(cc(ak),rcut(ak,:),arr)
       
       smaxlocal=0.
       do i=1,cc(ak)
          j=arr(i)
          ! avoid obstruction by nearby interpolated points
          if (rcut(ak,j) < sqrt(dx**2+dy**2)) cycle 
          if (slocal(ak,j)>smaxlocal) then
             smaxlocal = slocal(ak,j)
             visibility(celli(ak,j),cellj(ak,j)) = .true.
          endif
       end do
    end do
    
  end subroutine findallhorizon_wsort_v3



  subroutine horizon_core_wsort(x0,y0,h00,smax,surfaceSlope,azFac,i,j,h)
    ! very similar to horizon_core
    use filemanager, only : NSx, NSy, dx, dy
    use allinterfaces, only : horizontaldistance1, azimuth1, diffangle
    implicit none
    integer, intent(IN) :: i, j
    real(8), intent(IN) :: x0, y0, h00, surfaceslope, azFac, h(NSx,NSy)
    real(8), intent(INOUT) :: smax(naz)
    
    integer k, in, jn, ak, ak1, ak2, buf, akak
    real(8) az, az_neighbor, t, r, r_neighbor, hcut, d1, d2, d3
    real(8) s, slope_along_az
    integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
    integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)
    
    r = horizontaldistance1(i*dx,j*dy,x0,y0)
    if (r==0.) return
    az = azimuth1(x0,y0,i*dx,j*dy)  ! az should go from -pi ... pi
    
    if (floor(az*f)==ceiling(az*f)) then  ! grid point lies on ray
       ak = nint(az*f)+1
       
       cc(ak) = cc(ak)+1
       if (cc(ak)>CCMAX) error stop &
            & 'horizon_core_wsort: not enough memory allocated'
       s = (h(i,j)-h00)/r
       if (s>smax(ak)) smax(ak)=s
       
       slope_along_az = surfaceSlope*cos(azFac-azRay(ak))
       !angle_along_az = atan(slope_along_az)
       
       rcut(ak,cc(ak)) = r
       !slocal(ak,cc(ak)) = tan(atan(s)-atan(slope_along_az))
       slocal(ak,cc(ak)) = (s-slope_along_az)/(1+s*slope_along_az)

       ! atan(slocal) is not the same as asin(cosv) because it considers
       ! the slope in a *vertical* plane whereas cos(v) is 3D
       
       celli(ak,cc(ak))=i; cellj(ak,cc(ak))=j
       return
    endif
    
    do k=1,8
       in = i+ex(k)
       jn = j+ey(k)
       !if (in<1 .or. in>NSx .or. jn<1 .or. jn>NSy) cycle
       if (horizontaldistance1(x0,y0,in*dx,jn*dy) == 0.) cycle
       !if ((in-i0)**2+(jn-j0)**2 <= 1) cycle  ! stricter than in findallhorizon
       az_neighbor = azimuth1(x0,y0,in*dx,jn*dy)
       
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
       if (ak1>naz .or. ak2>naz) error stop &
            & 'horizon_core_wsort: index out of bound'
       
       d3 = diffangle(az,az_neighbor)
       do akak=ak1,ak2
          ak = akak; if (ak<=0) ak = ak+naz
          
          d1 = diffangle(az,azRay(ak))
          d2 = diffangle(az_neighbor,azRay(ak))
        
          if (d1+d2 <= d3+1.d-5) then
             if (d1>0.5*d3 .and. d3>1.d-6) cycle
             ! in findallhorizon 0.5 is 1.0 instead,
             ! but this leads to missing visibilities along zero azimth (ak=1)
             cc(ak) = cc(ak)+1
             if (cc(ak)>CCMAX) error stop &
                  & 'horizon_core_wsort: not enough memory allocated'
             
             r_neighbor = horizontaldistance1(in*dx,jn*dy,x0,y0)
             ! edge between h1,i0,j0 and h2,in,jn
             if (d3>1.d-6) then
                t = d1/d3  ! approximation
             else
                t = 0.5  ! dirty fix
             endif
             hcut = h(i,j)*(1-t) + h(in,jn)*t
             rcut(ak,cc(ak)) = r*(1-t) + r_neighbor*t  ! could be improved
             s = (hcut-h00) / rcut(ak,cc(ak))
             if (s>smax(ak)) smax(ak)=s
             
             slope_along_az = surfaceSlope*cos(azFac-azRay(ak))
             !slocal(ak,cc(ak)) = tan(atan(s)-atan(slope_along_az))
             slocal(ak,cc(ak)) = (s-slope_along_az)/(1+s*slope_along_az)
             
             celli(ak,cc(ak))=i; cellj(ak,cc(ak))=j
          endif
        
       end do  ! end of ak loop
    end do  ! end of k loop
    
  end subroutine horizon_core_wsort

END MODULE findvisibletopo



subroutine findviewfactors(h,i0,j0,unit,visibility)
  ! calculate subtended spherical angles and view factors of
  !     all quadrangles visible from (i0*dx,j0*dy,h(i0,j0))
  ! write view factors to file
  use filemanager
  use allinterfaces, except_this_one => findviewfactors
  implicit none
  integer, intent(IN) :: i0, j0, unit
  real(8), intent(IN) :: h(NSx,NSy)
  logical, intent(IN) :: visibility(NSx,NSy)
  real(8), parameter :: pi=3.1415926535897932
  integer i, j, cc
  real(8) r, dOh, landsize, cosv, VF, viewsize
  real(8) surfaceSlope, azFac
  integer, parameter :: CCMAX = NSx*NSy 
  integer, dimension(CCMAX) :: cellx, celly
  real(8), dimension(CCMAX) :: dOstack, VFstack
  logical, parameter :: full = .false.  ! output visible and invisible cells
  
  cc=0

  call difftopo1(NSx,NSy,i0,j0,h,dx,dy,surfaceSlope,azFac)
  
  do i=2,NSx-1
     if (.not.full) then
        !r = dx*abs(i-i0)
        r = horizontaldistance1(i*dx,1*dy,i0*dx,1*dy)
        if (r>RMAX) cycle  ! to save computational cost
     end if
     do j=2,NSy-1
        dOh = 0.
        VF = 0.
        if (i==i0 .and. j==j0 .and. .not.full) cycle
        if (.not.visibility(i,j) .and. .not.full) cycle
        r = horizontaldistance1(i*dx,j*dy,i0*dx,j0*dy)
        if (r>RMAX .and. .not.full) cycle  ! to save computational cost

        dOh = spherical_area(h,i,j,i0,j0)
        
        if (i==i0 .and. j==j0 .and. full) dOh = 0.
        
        ! cos(v)
        cosv = cos_viewing_angle(i0*dx,j0*dy,h(i0,j0), &
             & surfaceSlope,azFac,i*dx,j*dy,h(i,j))
        VF = dOh*cosv/pi  ! view factor

        !write(33,'(4(i5,1x),f7.2,1x,2(g10.4,1x),l)') &
        !     & i0,j0,i,j,h(i,j),dOh,VF,visibility(i,j)

        if (VF>0. .and. visibility(i,j)) then
           cc = cc+1
           cellx(cc)=i; celly(cc)=j
           dOstack(cc) = dOh
           VFstack(cc) = VF
        elseif (full) then
           cc = cc+1
           dOstack(cc) = 0.
           VFstack(cc) = 0.
        endif
        
     end do
  end do
  if (full .and. cc /= (NSx-2)*(NSy-2) ) error stop 'entries missing'
  
  landsize = sum(dOstack(1:cc))   ! 2*pi - (size of sky)
  viewsize = sum(VFstack(1:cc)) 

  write(unit,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') &
       & i0, j0, cc, viewsize
  do i=1,cc
     if (full) then
        if (VFstack(i)/=0.) then
           write(unit,'(g10.4,1x)',advance='no') VFstack(i)
        else ! compact zeros
           write(unit,'(f2.0,1x)',advance='no') 0.
        end if
     else
        write(unit,'(2(i5,1x),g10.4,1x)',advance='no') &
             & cellx(i),celly(i),VFstack(i)
     end if
  end do
  write(unit,"('')")
  
end subroutine findviewfactors



function spherical_area(h,i,j,i0,j0)
  use filemanager, only : NSx, NSy, dx, dy
  use allinterfaces, except_this_one => spherical_area
  implicit none
  real(8) spherical_area
  real(8), intent(IN) :: h(NSx,NSy)
  integer, intent(IN) :: i, j, i0, j0
  real(8), dimension(4) :: hq, xq, yq, theta, phi  ! quadrangle corners
  
  !call xyz2thetaphi(dx*(i-i0),dy*(j-j0),h(i,j)-h(i0,j0),thetac,phic)

!-get quadrangle corners
  ! upper right
  hq(1) = (h(i+1,j)+h(i+1,j+1)+h(i,j+1)+h(i,j))/4.
  xq(1) = dx*(i-i0+1./2.)
  yq(1) = dy*(j-j0+1./2.)
  
  ! upper left
  hq(2) = (h(i-1,j)+h(i-1,j+1)+h(i,j+1)+h(i,j))/4.
  xq(2) = dx*(i-i0-1./2.)
  yq(2) = dy*(j-j0+1./2.)
  
  ! lower left
  hq(3) = (h(i-1,j)+h(i-1,j-1)+h(i,j-1)+h(i,j))/4.
  xq(3) = dx*(i-i0-1./2.)
  yq(3) = dy*(j-j0-1./2.)
  
  ! lower right
  hq(4) = (h(i+1,j)+h(i+1,j-1)+h(i,j-1)+h(i,j))/4.
  xq(4) = dx*(i-i0+1./2.)
  yq(4) = dy*(j-j0-1./2.)
  
  hq(:) = hq(:) - h(i0,j0)
        
!-calculate spherical angle
  call xyz2thetaphi(xq,yq,hq,theta,phi) ! elemental

  spherical_area = area_spherical_quadrangle(phi,theta)
end function spherical_area



pure subroutine refinevisibility(i0,j0,h,visibility)
!***********************************************************************
! refinevisibility: This correction is necessary because azimuth rays 
!    for horizon calculation are not the same as azimuth rays connecting 
!    surface elements
!***********************************************************************
  use filemanager, only : NSx, NSy, dx, dy
  use allinterfaces, only : difftopo1, cos_viewing_angle
  implicit none
  integer, intent(IN) :: i0, j0
  real(8), intent(IN) :: h(NSx,NSy)
  logical, intent(INOUT) :: visibility(NSx,NSy)
  integer ii, jj
  real(8) cosv  ! cos(v)
  real(8) surfaceSlope00, azFac00, x0, y0, h00
  real(8) surfaceSlope, azFac, xB, yB, hB

  x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)
  call difftopo1(NSx,NSy,i0,j0,h,dx,dy,surfaceSlope00,azFac00)

  do ii=1,NSx
     do jj=1,NSy
        if (.not. visibility(ii,jj)) cycle

        xB = ii*dx; yB = jj*dy; hB = h(ii,jj)
        
        cosv = cos_viewing_angle(x0,y0,h00,surfaceSlope00,azFac00,xB,yB,hB)
        if (cosv<=0.) then
           visibility(ii,jj) = .false.
           cycle
        end if
        
        call difftopo1(NSx,NSy,ii,jj,h,dx,dy,surfaceSlope,azFac)
        cosv = cos_viewing_angle(xB,yB,hB,surfaceSlope,azFac,x0,y0,h00)
        ! sometimes happens for first surface element beyond cusp at horizon
        if (cosv<0.) visibility(ii,jj) = .false.
     end do
  end do
end subroutine refinevisibility
