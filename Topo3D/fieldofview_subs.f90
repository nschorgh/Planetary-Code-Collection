MODULE findvisibletopo
!***********************************************************************
! identify hidden versus visible topography
!***********************************************************************
  use filemanager, only : NSx, NSy, dx, dy

  integer, parameter, public :: naz = 360

  integer, parameter, private :: CCMAX = 6*(NSx+NSy) ! max # of elements in one azimuth ray
  integer, private :: cc(CCMAX)
  real(8), dimension(naz,CCMAX), private :: rcut, slocal 
  integer, dimension(naz,CCMAX), private :: celli, cellj 

  integer, private :: ak
  real(8), parameter, private :: pi = 3.1415926535897931
  real(8), parameter, private :: f = naz/(2*pi)
  real(8), parameter, public :: azRay(naz) = (/ ( (ak-1)/f, ak=1,naz) /)
  
contains


  subroutine findallhorizon_wsort_v3(h,i0,j0,smax,visibility)
    ! finds horizon and determines visibility for all azimuth rays
    use filemanager, only : NSx, NSy, RMAX, dx, dy
    use allinterfaces, only : horizontaldistance1, hpsort
    implicit none
    real(8), intent(IN) :: h(NSx,NSy)
    integer, intent(IN) :: i0, j0
    real(8), intent(OUT) :: smax(naz)
    logical, intent(OUT) :: visibility(NSx,NSy)
    integer i, j, ak
    real(8) surfaceSlope, azFac
    real(8) smaxlocal, x0, y0, h00
    integer, dimension(CCMAX) :: arr
    
    cc(:)=0
    smax(:)=0.

    call difftopo1(NSx,NSy,i0,j0,h,dx,dy,surfaceSlope,azFac)
    visibility(:,:) = .false.
    
    x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)
    
    do i=2,NSx-1
       if (horizontaldistance1(i*dx,1*dy,x0,1*dy)>RMAX) cycle  ! saves computations
       do j=2,NSy-1
          if (i==i0 .and. j==j0) cycle
          if (horizontaldistance1(i*dx,j*dy,x0,y0)>RMAX) cycle  ! saves computations
          
          call horizon_core_wsort(x0,y0,h00,smax,surfaceSlope,azFac,i,j,h)
          
       end do  ! end of j loop
    end do  ! end of i loop
    
    do ak = 1,naz
       ! sort by distance from (i0,j0)
       call hpsort(cc(ak),rcut(ak,:),arr)
       
       smaxlocal=0.
       do i=1,cc(ak)
          j=arr(i)
          if (rcut(ak,j) < sqrt(dx**2+dy**2)) cycle ! avoid obstruction by nearby interpolated points
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
    
    integer k,in,jn,ak,ak1,ak2,buf,akak
    real(8) az,az_neighbor,t,r,r_neighbor,hcut,d1,d2,d3
    real(8) s, slope_along_az
    integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
    integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)
    
    r = horizontaldistance1(i*dx,j*dy,x0,y0)
    if (r==0.) return
    az = azimuth1(x0,y0,i*dx,j*dy)  ! az should go from -pi ... pi
    
    if (floor(az*f)==ceiling(az*f)) then  ! grid point lies on ray
       ak = nint(az*f)+1
       
       cc(ak) = cc(ak)+1
       if (cc(ak)>CCMAX) error stop 'horizon_core_wsort: not enough memory allocated'
       s = (h(i,j)-h00)/r
       if (s>smax(ak)) smax(ak)=s
       
       slope_along_az = surfaceSlope*cos(azFac-azRay(ak))
       !angle_along_az = atan(slope_along_az)
       
       rcut(ak,cc(ak)) = r
       !slocal(ak,cc(ak)) = tan(atan(s)-angle_along_az))
       slocal(ak,cc(ak)) = (s-slope_along_az)/(1+s*slope_along_az)
       
       celli(ak,cc(ak))=i; cellj(ak,cc(ak))=j
       return
    endif
    
    do k=1,8
       in = i+ex(k)
       jn = j+ey(k)
       !if (in<1 .or. in>nx .or. jn<1 .or. jn>ny) cycle
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
       if (ak1>naz .or. ak2>naz) error stop 'horizon_core_wsort: index out of bound'
       
       d3=diffangle(az,az_neighbor)
       do akak=ak1,ak2
          ak = akak; if (ak<=0) ak = ak+naz
          
          d1=diffangle(az,azRay(ak))
          d2=diffangle(az_neighbor,azRay(ak))
        
          if (d1+d2<=d3+1.d-5) then
             if (d1>0.5*d3 .and. d3>1.d-6) cycle
             ! in findallhorizon 0.5 is 1.0 instead,
             ! but this leads to missing visibilities along zero azimth (ak=1)
             cc(ak) = cc(ak)+1
             if (cc(ak)>CCMAX) error stop 'horizon_core_wsort: not enough memory allocated'
             
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
             !slocal(ak,cc(ak)) = tan(atan(s)-angle_along_az)
             slocal(ak,cc(ak)) = (s-slope_along_az)/(1+s*slope_along_az)
             
             celli(ak,cc(ak))=i; cellj(ak,cc(ak))=j
          endif
        
       end do  ! end of ak loop
    end do  ! end of k loop
    
  end subroutine horizon_core_wsort


END MODULE findvisibletopo



subroutine find3dangle(h,i0,j0,unit,visibility)
  ! calculate subtended spherical angles and view factors of
  !     all quadrangles visible from (i0*dx,j0*dy,h(i0,j0))
  ! write view factors to file
  use filemanager
  use allinterfaces, except_this_one => find3dangle
  implicit none
  integer, intent(IN) :: i0, j0, unit
  real(8), intent(IN) :: h(NSx,NSy)
  logical, intent(IN) :: visibility(NSx,NSy)
  real(8), parameter :: pi=3.1415926535897931
  integer i, j, cc
  real(8) r, dOh, landsize, cosv, VF, viewsize
  real(8) surfaceSlope, azFac
  integer, parameter :: CCMAX = NSx*NSy 
  integer, dimension(CCMAX) :: cellx, celly
  real(8), dimension(CCMAX) :: dOstack, VFstack
  real(8), dimension(4) :: hq, xq, yq, theta, phi  ! quadrangle corners
  logical, parameter :: verbose = .false.

  cc=0

  call difftopo1(NSx,NSy,i0,j0,h,dx,dy,surfaceSlope,azFac)
  
  do i=2,NSx-1
     !r = dx*abs(i-i0)
     r = horizontaldistance1(i*dx,1*dy,i0*dx,1*dy)
     if (r>RMAX) cycle  ! to save computational cost
     do j=2,NSy-1
        dOh = 0.
        VF = 0.
        if (i==i0 .and. j==j0 .and. .not.verbose) cycle
        if (.not.visibility(i,j) .and. .not.verbose) cycle
        r = horizontaldistance1(i*dx,j*dy,i0*dx,j0*dy)
        if (r>RMAX .and. .not.verbose) cycle  ! to save computational cost

        !call xyz2thetaphi(dx*(i-i0),dy*(j-j0),h(i,j)-h(i0,j0),thetac,phic)

!-------get quadrangle corners
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
        
!-------calculate spherical angle
        call xyz2thetaphi(xq,yq,hq,theta,phi) ! elemental

        dOh = area_spherical_quadrangle(phi,theta)
        if (i==i0 .and. j==j0 .and. verbose) dOh = 0.
        
        if (verbose) then
           write(23,'(4(i5,1x),f7.2,1x,g10.4,1x,l)') &
                & i0,j0,i,j,h(i,j),dOh,visibility(i,j)
        endif

        ! cos(v)
        cosv = cos_viewing_angle(i0*dx,j0*dy,h(i0,j0),surfaceSlope,azFac,i*dx,j*dy,h(i,j))
        VF = dOh*cosv/pi  ! view factor

        !if (dOh<0.) stop 'Does this ever happen?'
        if (dOh>0.) then
           cc = cc+1   
           if (cc>CCMAX) stop 'find3dangle: not enough memory allocated'
           cellx(cc)=i; celly(cc)=j
           dOstack(cc) = dOh
           VFstack(cc) = VF
        endif
        
     end do
  end do

  landsize = sum(dOstack(1:cc))   ! 2*pi - (size of sky)
  viewsize = sum(VFstack(1:cc)) 

  !write(unit-1,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') i0, j0, cc, landsize
  !do i=1,cc
  !   write(unit-1,'(2(i5,1x),g10.4,1x)',advance='no') cellx(i),celly(i),dOstack(i)
  !end do
  !write(unit-1,"('')")
  
  write(unit,'(2(i5,1x),i6,1x,2(f7.5,1x))',advance='no') i0, j0, cc, landsize, viewsize
  do i=1,cc
     write(unit,'(2(i5,1x),g10.4,1x)',advance='no') cellx(i),celly(i),VFstack(i)
  end do
  write(unit,"('')")
  
end subroutine find3dangle



pure subroutine difftopo1(NSx,NSy,i,j,h,dx,dy,surfaceSlope,az)
! calculate slope and azimuth of surface element
  implicit none
  integer, intent(IN) :: NSx, NSy, i, j
  real(8), intent(IN) :: h(NSx,NSy), dx, dy
  real(8), intent(OUT) :: surfaceSlope, az
  real(8) sx, sy
  
  sx = -1e32; sy = -1e32  ! avoids compiler warning

  if (i>1 .and. i<NSx) then
     sx = (h(i+1,j)-h(i-1,j))/(2.*dx)
  else
    if (i==1)   sx = (h(i+1,j)-h(i,j))/dx
    if (i==NSx) sx = (h(i,j)-h(i-1,j))/dx
  endif

  if (j>1 .and. j<NSy) then
     sy = (h(i,j+1)-h(i,j-1))/(2.*dy)
  else
     if (j==1)   sy = (h(i,j+1)-h(i,j))/dy
     if (j==NSy) sy = (h(i,j)-h(i,j-1))/dy
  endif

  surfaceSlope = sqrt(sx**2+sy**2)  
  az = atan2(sx,-sy)  ! north is up, clockwise
end subroutine difftopo1



elemental function cos_viewing_angle(x0,y0,h00,surfaceSlope,azFac,xB,yB,hB)
!***********************************************************************
!  function that calculates angle between surface normal at (x0,y0,h00)
!     and vector pointing to (xB,yB,h(xB,yB)) as in horizon_core_wsort
!***********************************************************************
  use allinterfaces, only : horizontaldistance1, azimuth1
  implicit none
  real(8) cos_viewing_angle
  real(8), intent(IN) :: x0, y0, h00, surfaceSlope, azFac
  real(8), intent(IN) :: xB, yB, hB
  real(8) az, s, r, slope_along_az

  r = horizontaldistance1(xB,yB,x0,y0)
  az = azimuth1(x0,y0,xB,yB)
  s = (hB-h00)/r

  slope_along_az = surfaceSlope*cos(azFac-az)

  !viewing_angle = pi/2 - atan(s) + atan(slope_along_az)
  cos_viewing_angle = (s-slope_along_az)/sqrt(1+s**2)/sqrt(1+slope_along_az**2)
end function cos_viewing_angle



elemental subroutine xyz2thetaphi(x,y,z,theta,phi)
  ! cartesian -> polar coordinates
  implicit none
  real(8), intent(IN) :: x,y,z
  real(8), intent(OUT) :: theta,phi
  theta = acos(z/sqrt(x**2+y**2+z**2))
  phi = atan2(y,x)
end subroutine xyz2thetaphi



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
