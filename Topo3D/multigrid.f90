! subroutines for multigrid acceleration


subroutine downsample(NSx,NSy,h,hhalf)
  implicit none
  integer, intent(IN) :: NSx,NSy
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT) :: hhalf(NSx/2,NSy/2) ! new dimensions
  integer :: x(NSx/2),y(NSy/2)
  integer i2,j2,ii,jj,k,r
  real(8) psum, w(NSx/2,NSy/2)
  integer, parameter :: ei(5)= (/0, 1, 0, -1, 0 /) ! counter-clockwise
  integer, parameter :: ej(5)= (/0, 0, 1, 0, -1 /)

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
           r = 1+3*(ei(k)**2+ej(k)**2)  ! should be 1 or 4
           psum = psum + h(ii,jj)/r           
           w(i2,j2) = w(i2,j2) + 1./r
        enddo
        hhalf(i2,j2) = psum/w(i2,j2)

     enddo
  enddo
  !do j2=1,NSy/2
  !   write(6,'(9999(f3.1,1x))') w(:,j2)
  !enddo
  !write(6,*)
end subroutine downsample



real(8) elemental function horizontaldistance1(x1,y1,x2,y2)
  ! distance between two points; must have same units as height
  ! as in crater_common.f90, but distance based
  implicit none
  real(8), intent(IN) :: x1,y1,x2,y2
  horizontaldistance1 = sqrt((x1-x2)**2+(y1-y2)**2)
  !if (horizontaldistance1>0. .and. horizontaldistance1<1d-6) stop 'nogood'
end function horizontaldistance1



real(8) elemental function azimuth1(x1,y1,x2,y2)
  ! as in crater_common.f90, but distance based
  implicit none
  real(8), intent(IN) :: x1,y1,x2,y2
  azimuth1 = atan2(x2-x1,-(y2-y1)) 
end function azimuth1



subroutine findallhorizon_MG1(h,i0,j0,naz,smax)
  ! find all horizon heights, without use of multigrid
  use filemanager, only : NSx,NSy,dx,dy
  use allinterfaces, only : horizon_MG_core
  implicit none
  integer, intent(IN) :: i0,j0,naz
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT) :: smax(naz)
  integer i,j
  real(8) x0,y0,h00

  smax=0.
  x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)
  do i=2,NSx-1
     do j=2,NSy-1
        call horizon_MG_core(x0,y0,h00,naz,smax,i,j,h,1)
     end do
  end do

end subroutine findallhorizon_MG1



subroutine findallhorizon_MG3(h,h2,h3,i0,j0,naz,smax,RMG)
  ! based on shadow_subs.f90, but for multigrid method
  use filemanager, only : NSx,NSy,dx,dy
  use allinterfaces, only : horizontaldistance1, horizon_MG_core
  implicit none
  integer, intent(IN) :: i0,j0,naz
  real(8), intent(IN) :: h(NSx,NSy),h2(NSx/2,NSy/2),h3(NSx/4,NSy/4)
  real(8), intent(IN) :: RMG
  real(8), intent(OUT) :: smax(naz)
  integer i1,j1,i2,j2,i3,j3
  real(8) r,x0,y0,h00
  logical :: VERBOSE = .false.

  smax=0.

  x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)

  ! start with the coarsest grid
  ! if distance is too close, switch to finer grid
  do i3=1,NSx/4 + 1  ! +1 so a last odd one gets included too
     !if (sin(azRay)*(4*i3-i0) < 0.) cycle  ! TO BE TESTED
     do j3=1,NSy/4 + 1  ! +1 so a last odd one gets included too
        r = horizontaldistance1(4*i3*dx,4*j3*dy,x0,y0)
        !r = horizontaldistance(4*i3,4*j3,i0,j0)
        if (r>=2*RMG) then  ! do 1 coarse grid cell
           call horizon_MG_core(x0,y0,h00,naz,smax,i3,j3,h3,4)
           if (VERBOSE) write(6,*) 'Level3',4*i3,4*j3,smax(90)
        else ! do 4 finer cells
           do i2=2*i3-1,2*i3
              do j2=2*j3-1,2*j3
                 r = horizontaldistance1(2*i2*dx,2*j2*dy,x0,y0)
                 !r = horizontaldistance(2*i2,2*j2,i0,j0)
                 if (r>=RMG) then
                    call horizon_MG_core(x0,y0,h00,naz,smax,i2,j2,h2,2)
                    if (VERBOSE) write(6,*) 'Level2',2*i2,2*j2,smax(90)
                 else ! do 4 finest cells
                    do i1=2*i2-1,2*i2
                       do j1=2*j2-1,2*j2
                          if (i1<=1 .or. j1<=1 .or. i1>=NSx .or. j1>=NSy) cycle
                          call horizon_MG_core(x0,y0,h00,naz,smax,i1,j1,h,1)
                          if (VERBOSE) write(6,*) 'Level1',i1,j1,smax(90)
                       enddo
                    enddo
                 endif
              enddo
           enddo
        endif
     enddo
  enddo
end subroutine findallhorizon_MG3



subroutine horizon_MG_core(x0,y0,h00,naz,smax,i,j,h,P)
  use filemanager, only : NSx,NSy,dx,dy
  use allinterfaces, only : horizontaldistance1, azimuth1, diffangle
  implicit none
  real(8), intent(IN) :: x0,y0,h00  ! on fine grid
  integer, intent(IN) :: naz
  real(8), intent(INOUT) :: smax(naz)
  integer, intent(IN) :: i,j,P  ! fine or coarse grid
  real(8), intent(IN) :: h(NSx/P,NSy/P) ! fine or coarse grid

  integer k,in,jn,ak,ak1,ak2,buf,akak
  integer nx,ny    ! grid-level specific grid size
  real(8) az,az_neighbor,t,r,r_neighbor,rcut,s,hcut,d1,d2,d3
  real(8) dxl,dyl  ! grid-level specific resolution
  integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
  integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)
  real(8) f,azRay(naz)
  real(8), parameter :: pi=3.1415926535897931

  nx = NSx/P; ny = NSy/P
  dxl = P*dx; dyl = P*dy
  ! no grid-level specific variables are used below

  if (i>nx .or. j>ny) return

  r = horizontaldistance1(i*dxl,j*dyl,x0,y0)
  if (r==0.) return
  az = azimuth1(x0,y0,i*dxl,j*dyl)

  f = naz/(2*pi)
  azRay = (/ ( (ak-1)/f, ak=1,naz) /)
  
  if (floor(az*f)==ceiling(az*f)) then  ! grid point lies on ray
     ak=nint(az*f)+1
     s = (h(i,j)-h00)/r
     if (s>smax(ak)) smax(ak)=s
     return
  endif
  
  do k=1,8
     in = i+ex(k)
     jn = j+ey(k)
     if (in>nx .or. jn>ny) cycle
     if (horizontaldistance1(x0,y0,in*dxl,jn*dyl) == 0.) cycle
     az_neighbor = azimuth1(x0,y0,in*dxl,jn*dyl)

     if (az >= az_neighbor) then 
        ak1 = floor(az*f)+1
        ak2 = ceiling(az_neighbor*f)+1
     else
        ak1 = ceiling(az*f)+1
        ak2 = floor(az_neighbor*f)+1
     endif
     if (abs(ak1-naz/2)<=P/2 .and. ak2<0) ak2 = ak2+naz
     if (abs(ak2-naz/2)<=P/2 .and. ak1<0) ak1 = ak1+naz
     if (ak2<ak1) then ! swap
        buf=ak1; ak1=ak2; ak2=buf;
     endif
     if (ak1>naz .or. ak2>naz) stop 'horizon_MG_core: index out of bound'

     d3=diffangle(az,az_neighbor)
     do akak=ak1,ak2
     !do akak=1,naz
        ak = akak; if (ak<=0) ak = ak+naz
        if (azRay(ak)==180) print *,'bad ray'
        d1=diffangle(az,azRay(ak))
        d2=diffangle(az_neighbor,azRay(ak))
        if (d1+d2<=d3+1.d-5) then  
           if (d1>1.0*d3 .and. d3>1.d-6) cycle 
           r_neighbor = horizontaldistance1(in*dxl,jn*dyl,x0,y0)
           ! edge between h1,i0,j0 and h2,in,jn
           if (d3>1.d-6) then
              t = d1/d3  ! approximation
           else
              t = 0.5  ! dirty fix
           endif
           hcut = h(i,j)*(1-t)+h(in,jn)*t
           rcut = r*(1-t)+r_neighbor*t
           s = (hcut-h00)/rcut
           if (s>smax(ak)) smax(ak)=s
        endif
     end do
  end do

end subroutine horizon_MG_core
