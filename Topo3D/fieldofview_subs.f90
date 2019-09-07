subroutine findallhorizon_wsort(h,i0,j0,naz,smax,visibility)
! finds horizon and determines visibility for all azimuth rays
  use filemanager, only : NSx,NSy,RMAX
  use allinterfaces, only : horizontaldistance, azimuth, diffangle, hpsort
  implicit none
  integer, intent(IN) :: i0,j0,naz
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT) :: smax(naz)
  logical, intent(OUT) :: visibility(NSx,NSy)
  integer, parameter :: CCMAX = 6*(NSx+NSy) ! max # of elements on azimuth ray
  integer i,j,k,in,jn,cc(CCMAX), ak,ak1,ak2,buf,akak
  real(8) az,az_neighbor,t,r,r_neighbor,hcut,d1,d2,d3, f,azRay(naz)
  real(8) s, smaxlocal, surfaceSlope, azFac, slope_along_az
  real(8), dimension(:,:), allocatable :: rcut, slocal
  integer, dimension(:,:), allocatable :: celli, cellj
  integer, dimension(CCMAX) :: arr
  integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
  integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)
  real(8), parameter :: pi=3.1415926535897931

  allocate(rcut(naz,CCMAX), slocal(naz,CCMAX), celli(naz,CCMAX), cellj(naz, CCMAX))
  
  f = naz/(2*pi)
  
  ! azimuth to integer mapping
  do ak=1,naz
     azRay(ak) = (ak-1)/f   ! inverse mapping:  ak = azRay*f+1
  enddo
  
  cc(:)=0
  smax(:)=0.

  call difftopo1(i0,j0,h,surfaceSlope,azFac)
  visibility(:,:) = .false.
  
  do i=2,NSx-1
     if (horizontaldistance(i,1,i0,1)>RMAX) cycle  ! saves computations
     do j=2,NSy-1
        if (i==i0 .and. j==j0) cycle
        r = horizontaldistance(i,j,i0,j0)
        if (r>RMAX) cycle  ! saves computations
        az = azimuth(i0,j0,i,j)  ! az should go from -pi ... pi

        if (floor(az*f)==ceiling(az*f)) then  ! grid point lies on ray
           ak=nint(az*f)+1

           cc(ak) = cc(ak)+1
           if (cc(ak)>CCMAX) stop 'findallhorizons_wsort: not enough memory allocated'
           s = (h(i,j)-h(i0,j0))/r
           if (s>smax(ak)) smax(ak)=s

           slope_along_az = surfaceSlope*cos(azFac-azRay(ak))
           !angle_along_az = atan(slope_along_az)
           
           rcut(ak,cc(ak)) = r
           !slocal(ak,cc(ak)) = tan(atan(s)-angle_along_az))
           slocal(ak,cc(ak)) = (s-slope_along_az)/(1+s*slope_along_az)

           celli(ak,cc(ak))=i; cellj(ak,cc(ak))=j
           cycle
        endif
        
        do k=1,8
           in = i+ex(k)
           jn = j+ey(k)
           if (in==i0 .and. jn==j0) cycle
           !if ((in-i0)**2+(jn-j0)**2 <= 1) cycle  ! stricter than in findallhorizon
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
           if (ak1>naz .or. ak2>naz) stop 'findallhorizons_wsort: index out of bound'

           d3=diffangle(az,az_neighbor)
           do akak=ak1,ak2
              ak = akak; if (ak<=0) ak = ak+naz

              d1=diffangle(az,azRay(ak))
              d2=diffangle(az_neighbor,azRay(ak))
           
              if (d1+d2<=d3+1.d-5) then
                 if (d1>0.5*d3 .and. d3>1.d-6) cycle
                 ! in findallhorizon 0.5 is 1.0 instead,
                 ! but this leads to missing visibilities along zero azimth (ak=1)
                 cc(ak)=cc(ak)+1
                 if (cc(ak)>CCMAX) stop 'findonehorizons_wsort: not enough memory allocated'
                 
                 r_neighbor = horizontaldistance(in,jn,i0,j0)
                 ! edge between h1,i0,j0 and h2,in,jn
                 if (d3>1.d-6) then
                    t = d1/d3  ! approximation
                 else
                    t = 0.5  ! dirty fix
                 endif
                 hcut = h(i,j)*(1-t)+h(in,jn)*t
                 rcut(ak,cc(ak)) = r*(1-t)+r_neighbor*t  ! could be improved
                 s = (hcut-h(i0,j0))/rcut(ak,cc(ak))
                 if (s>smax(ak)) smax(ak)=s

                 slope_along_az = surfaceSlope*cos(azFac-azRay(ak))
                 !slocal(ak,cc(ak)) = tan(atan(s)-angle_along_az)
                 slocal(ak,cc(ak)) = (s-slope_along_az)/(1+s*slope_along_az)
                 
                 celli(ak,cc(ak))=i; cellj(ak,cc(ak))=j
              endif
              
           end do  ! end of ak loop
        end do  ! end of k loop
        
     enddo  ! end of j loop
  enddo  ! end of i loop

  do ak = 1,naz
     ! sort by distance from (i0,j0)
     call hpsort(cc(ak),rcut(ak,:),arr)
     
     smaxlocal=0.
     do i=1,cc(ak)
        j=arr(i)
        if (slocal(ak,j)>smaxlocal) then
           smaxlocal=slocal(ak,j)
           visibility(celli(ak,j),cellj(ak,j))=.true.
        endif
     enddo
  end do
  
end subroutine findallhorizon_wsort



subroutine find3dangle(h,i0,j0,unit,visibility)
  ! calculate subtended spherical angle and view factor
  use filemanager
  use allinterfaces, except_this_one => find3dangle
  implicit none
  integer, intent(IN) :: i0,j0,unit
  real(8), intent(IN) :: h(NSx,NSy)
  logical, intent(IN) :: visibility(NSx,NSy)
  real(8), parameter :: pi=3.1415926535897931
  integer i, j, k, cc
  real(8) r, thetac, phic, dOh, landsize, cosv, viewsize
  integer, parameter :: CCMAX = NSx*NSy 
  integer, dimension(CCMAX) :: cellx, celly
  !real(8), dimension(CCMAX) :: thetastack, phistack
  real(8), dimension(CCMAX) :: dOstack, VFstack
  real(8), dimension(4) :: hq,xq,yq,theta,phi  ! quadrangle centers
  logical, parameter :: verbose = .false.

  cc=0
  do i=2,NSx-1
     if (dx*abs(i-i0)>RMAX) cycle  ! to save computational cost
     do j=2,NSy-1
        dOh=0.
        if (i==i0 .and. j==j0 .and. .not.verbose) cycle
        if (.not.visibility(i,j) .and. .not.verbose) cycle
        r = sqrt(dx*dx*(i-i0)**2+dy*dy*(j-j0)**2)
        if (r>RMAX .and. .not.verbose) cycle  ! to save computational cost

        call xyz2thetaphi(dx*(i-i0),dy*(j-j0),h(i,j)-h(i0,j0),thetac,phic)

!-------get quadrangle centers
        ! upper right
        hq(1)=(h(i+1,j)+h(i+1,j+1)+h(i,j+1)+h(i,j))/4.
        xq(1)=dx*(i-i0+1./2.)
        yq(1)=dy*(j-j0+1./2.)

        ! upper left
        hq(2)=(h(i-1,j)+h(i-1,j+1)+h(i,j+1)+h(i,j))/4.
        xq(2)=dx*(i-i0-1./2.)
        yq(2)=dy*(j-j0+1./2.)

        ! lower left
        hq(3)=(h(i-1,j)+h(i-1,j-1)+h(i,j-1)+h(i,j))/4.
        xq(3)=dx*(i-i0-1./2.)
        yq(3)=dy*(j-j0-1./2.)

        ! lower right
        hq(4)=(h(i+1,j)+h(i+1,j-1)+h(i,j-1)+h(i,j))/4.
        xq(4)=dx*(i-i0+1./2.)
        yq(4)=dy*(j-j0-1./2.)

        hq = hq - h(i0,j0)

!-------get spherical angle
        do k=1,4
           call xyz2thetaphi(xq(k),yq(k),hq(k),theta(k),phi(k))
        enddo

        if (visibility(i,j)) then
           dOh = area_spherical_quadrangle(phi,theta)
        else
           dOh = 0.  ! for verbose mode
        end if
        if (i==i0 .and. j==j0 .and. verbose) then
           dOh=0.
           thetac=3.1416/2.; phic=0.
        endif
        
        if (verbose) then
           !write(23,'(4(i5,1x),f7.2,1x,g10.4,1x,l)') &
           !     & i0,j0,i,j,h(i,j),dOh,visibility(i,j)
           if (visibility(i,j)) then
              write(23,'(4(i5,1x),f7.2,1x,g10.4,1x,i1)') i0,j0,i,j,h(i,j),dOh,1
           else
              write(23,'(4(i5,1x),f7.2,1x,g10.4,1x,i1)') i0,j0,i,j,h(i,j),dOh,0
           endif
        endif
        
        if (dOh>0.) then
           cc = cc+1   
           if (cc>CCMAX) stop 'find3dangle: not enough memory allocated'
           cellx(cc)=i; celly(cc)=j
           !thetastack(cc)=thetac; phistack(cc)=phic  ! for optional output
           dOstack(cc)=dOh

           cosv = cos_viewing_angle(i0,j0,i,j,h)  ! cos(v)
           VFstack(cc) = dOh*cosv/pi  ! view factor
        endif
     enddo
  enddo

  landsize = sum(dOstack(1:cc))   ! 2*pi- size of sky
  viewsize = sum(VFstack(1:cc)) 

  !write(unit-1,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') i0, j0, cc, landsize
  !do i=1,cc
  !   write(unit-1,'(2(i5,1x),g10.4,1x)',advance='no') cellx(i),celly(i),dOstack(i)
  !enddo
  !write(unit-1,"('')")

  write(unit,'(2(i5,1x),i6,1x,2(f7.5,1x))',advance='no') i0, j0, cc, landsize, viewsize
  do i=1,cc
     write(unit,'(2(i5,1x),g10.4,1x)',advance='no') cellx(i),celly(i),VFstack(i)
  enddo
  write(unit,"('')")
  
  if (minval(cellx(1:cc))<1 .or. minval(celly(1:cc))<1) stop 'find3dangle: index out of boud'
  if (maxval(cellx(1:cc))>NSx .or. maxval(celly(1:cc))>NSy) stop 'find3dangle: index out of boud'
end subroutine find3dangle



pure subroutine difftopo1(i,j,h,surfaceSlope,az)
  ! calculate slope and azimuth of surface element
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



pure function cos_viewing_angle(i0,j0,i,j,h)
!***********************************************************************
!  function that calculates angle between surface normal at (i,j)
!     and vector pointing to (ii,jj) as in findallhorizon_wsort
!***********************************************************************
  use filemanager, only : NSx,NSy
  use allinterfaces, except_this_one => cos_viewing_angle
  implicit none
  real(8) cos_viewing_angle
  integer, intent(IN) :: i0,j0,i,j
  real(8), intent(IN) :: h(NSx,NSy)
  !real(8), parameter :: pi=3.1415926535897931
  real(8) az, s, r, surfaceSlope, azFac, slope_along_az

  !r = sqrt(dx*dx*(i-i0)**2+dy*dy*(j-j0)**2)
  r = horizontaldistance(i,j,i0,j0)
  az = azimuth(i0,j0,i,j)
  s = (h(i,j)-h(i0,j0))/r

  call difftopo1(i0,j0,h,surfaceSlope,azFac)
  slope_along_az = surfaceSlope*cos(azFac-az)

  !viewing_angle = pi/2 - atan(s) + atan(slope_along_az)
  cos_viewing_angle = (s-slope_along_az)/sqrt(1+s**2)/sqrt(1+slope_along_az**2)
end function cos_viewing_angle



elemental subroutine xyz2thetaphi(x,y,z,theta,phi)
  implicit none
  real(8), intent(IN) :: x,y,z
  real(8), intent(OUT) :: theta,phi
  theta=acos(z/sqrt(x**2+y**2+z**2))
  phi=atan2(y,x)
end subroutine xyz2thetaphi



subroutine refinevisibility(i0,j0,h,visibility)
!***********************************************************************
! refinevisibility: This correction is necessary because azimuth rays 
!    for horizon calculation are not the same as azimuth rays connecting 
!    surface elements
!***********************************************************************
  use filemanager, only : NSx,NSy
  use allinterfaces, only : cos_viewing_angle
  implicit none
  integer, intent(IN) :: i0,j0
  real(8), intent(IN) :: h(NSx,NSy)
  logical, intent(INOUT) :: visibility(NSx,NSy)
  integer ii,jj
  real(8) cosv  ! cos(v)

  do ii=1,NSx
     do jj=1,NSy
        cosv = cos_viewing_angle(i0,j0,ii,jj,h)
        if (cosv<=0.) visibility(ii,jj)=.false.
        cosv = cos_viewing_angle(ii,jj,i0,j0,h)
        ! sometimes happens for first surface element beyond cusp/horizon
        if (cosv<0.) visibility(ii,jj)=.false.
     enddo
  enddo
end subroutine refinevisibility
