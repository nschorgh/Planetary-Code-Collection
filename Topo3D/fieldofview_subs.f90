subroutine find3dangle(h,i0,j0,unit,visibility)
  use filemanager
  use allinterfaces, only : xyz2thetaphi, area_spherical_quadrangle
  implicit none
  integer, intent(IN) :: i0,j0,unit
  real(8), intent(IN) :: h(NSx,NSy)
  logical, intent(IN) :: visibility(NSx,NSy)
  integer i, j, k, cc
  real(8) r, thetac, phic, dOh, skysize
  integer, parameter :: CCMAX = NSx*NSy 
  integer, dimension(CCMAX) :: cellx, celly
  real(8), dimension(CCMAX) :: thetastack, phistack, dOstack
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
           thetastack(cc)=thetac; phistack(cc)=phic
           dOstack(cc)=dOh
        endif
     enddo
  enddo

  skysize = sum(dOstack(1:cc))   ! 2*pi- size of sky

  write(unit,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') i0, j0, cc, skysize
  do i=1,cc
     write(unit,'(2(i5,1x),g10.4,1x)',advance='no') cellx(i),celly(i),dOstack(i)
  enddo
  write(unit,"('')")
  if (minval(cellx(1:cc))<1 .or. minval(celly(1:cc))<1) stop 'find3dangle: index out of boud'
  if (maxval(cellx(1:cc))>NSx .or. maxval(celly(1:cc))>NSy) stop 'find3dangle: index out of boud'
end subroutine find3dangle



elemental subroutine xyz2thetaphi(x,y,z,theta,phi)
  implicit none
  real(8), intent(IN) :: x,y,z
  real(8), intent(OUT) :: theta,phi
  theta=acos(z/sqrt(x**2+y**2+z**2))
  phi=atan2(y,x)
end subroutine xyz2thetaphi



pure function area_spherical_quadrangle(phi,theta)
! area of two triangles on sphere
  implicit none
  real(8), intent(IN) :: phi(4), theta(4)
  real(8) E1, E2, area_spherical_quadrangle
  interface
     pure function area_spherical_triangle(phi,theta)
       real(8), intent(IN) :: phi(3), theta(3)
       real(8) area_spherical_triangle
     end function area_spherical_triangle
  end interface
  
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
  real(8), intent(IN) :: phi1,phi2,theta1,theta2  ! in radians
  real(8) distanceonsphere, buf, lat1, lat2
  real(8), parameter :: pi=3.1415926535897932
  lat1=pi/2-theta1; lat2=pi/2-theta2

  ! buf = square of half of cord length distance
  buf = sin((lat1-lat2)/2.)**2+cos(lat1)*cos(lat2)*sin((phi1-phi2)/2.)**2
  distanceonsphere = 2.*asin(sqrt(buf))
  !distanceonsphere = 2.*atan2(sqrt(buf),sqrt(1-buf))
end function distanceonsphere



subroutine refinevisibility(i0,j0,h,visibility)
!***********************************************************************
! refinevisibility: This correction is necessary, because azimuth rays 
!    for horizon calculation are not the same as azimuth rays connecting 
!    surface elements
!***********************************************************************
  use filemanager, only : NSx,NSy
  use allinterfaces, only : viewing_angle
  implicit none
  integer, intent(IN) :: i0,j0
  real(8), intent(IN) :: h(NSx,NSy)
  logical, intent(INOUT) :: visibility(NSx,NSy)
  integer ii,jj
  real(8) v

  do ii=1,NSx
     do jj=1,NSy
        v = viewing_angle(i0,j0,ii,jj,h)
        if (cos(v)<=0.) visibility(ii,jj)=.false.
        v = viewing_angle(ii,jj,i0,j0,h)
        ! sometimes happens for first surface element beyond cusp/horizon
        if (cos(v)<0.) visibility(ii,jj)=.false.
     enddo
  enddo
end subroutine refinevisibility
