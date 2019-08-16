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



real(8) elemental function horizontaldistance1(x1,y1,x2,y2)
  ! distance between two points; must have same units as height
  ! as in topo3d_common.f90, but distance based
  implicit none
  real(8), intent(IN) :: x1,y1,x2,y2
  horizontaldistance1 = sqrt((x1-x2)**2+(y1-y2)**2)
end function horizontaldistance1



real(8) elemental function azimuth1(x1,y1,x2,y2)
  ! as in topo3d_common.f90, but distance based
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



pure subroutine horizon_MG_core(x0,y0,h00,naz,smax,i,j,h,P)
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
     if (in<1 .or. jn<1) cycle
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
     if (ak1>naz .or. ak2>naz) error stop 'horizon_MG_core: index out of bound'

     d3=diffangle(az,az_neighbor)
     do akak=ak1,ak2
     !do akak=1,naz
        ak = akak; if (ak<=0) ak = ak+naz
        !if (azRay(ak)==180) print *,'bad ray'
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



module newmultigrid
  use filemanager, only : NSx,NSy,dx,dy
  real(8), allocatable :: h2(:,:), h3(:,:), h4(:,:), h5(:,:), h6(:,:)
  real(8), allocatable :: h7(:,:), h8(:,:), h9(:,:), h10(:,:)
  
contains
  subroutine downsample_all(h,LMAX,LACT)
    use allinterfaces, only: downsample
    implicit none
    integer, intent(IN) :: LMAX
    real(8), intent(IN) :: h(NSx,NSy)
    integer, intent(OUT) :: LACT  ! levels actually allocated

    LACT = 1
    if (LMAX<2) return
    allocate(h2(NSx/2,NSy/2))
    call downsample(NSx,NSy,h,h2)
    LACT = 2
    if (min(NSx,NSy)/4>10 .and. LMAX>=3) then
       allocate(h3(NSx/4,NSy/4))
       call downsample(NSx/2,NSy/2,h2,h3)
       LACT = 3
       if (min(NSx,NSy)/8>10 .and. LMAX>=4) then
          allocate(h4(NSx/8,NSy/8))
          call downsample(NSx/4,NSy/4,h3,h4)
          LACT = 4
          if (min(NSx,NSy)/16>10 .and. LMAX>=5) then
             allocate(h5(NSx/16,NSy/16));
             call downsample(NSx/8,NSy/8,h4,h5)
             LACT = 5
             if (min(NSx,NSy)/32>10 .and. LMAX>=6) then
                allocate(h6(NSx/32,NSy/32))
                call downsample(NSx/16,NSy/16,h5,h6)
                LACT = 6
                if (min(NSx,NSy)/64>10 .and. LMAX>=7) then
                   allocate(h7(NSx/64,NSy/64))
                   call downsample(NSx/32,NSy/32,h6,h7)
                   LACT = 7
                   if (min(NSx,NSy)/128>10 .and. LMAX>=8) then
                      allocate(h8(NSx/128,NSy/128))
                      call downsample(NSx/64,NSy/64,h7,h8)
                      LACT = 8
                      if (min(NSx,NSy)/256>10 .and. LMAX>=9) then
                         allocate(h9(NSx/256,NSy/256))
                         call downsample(NSx/128,NSy/128,h8,h9)
                         LACT = 9
                         if (min(NSx,NSy)/512>10 .and. LMAX>=10) then
                            allocate(h10(NSx/512,NSy/512))
                            call downsample(NSx/256,NSy/256,h9,h10)
                            LACT = 10
                         endif
                      endif
                   endif
                endif
             endif
          endif
       endif
    endif
  end subroutine downsample_all

  
  subroutine findallhorizon_MGR(h,i0,j0,naz,smax,RMG,L)
    ! based on shadow_subs.f90, but for multigrid method
    ! recursive implementation
    use allinterfaces, only: horizontaldistance1, horizon_MG_core
    implicit none
    integer, intent(IN) :: i0,j0,naz,L
    real(8), intent(IN) :: h(NSx,NSy)
    real(8), intent(IN) :: RMG
    real(8), intent(OUT) :: smax(naz)
    integer P,ii,jj
    real(8) r,x0,y0,h00
    logical, parameter :: VERBOSE = .false.

    smax=0.
    
    x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)
    
    P=2**(L-1)
    if (L>10 .or. L<2) error stop 'findallhorizon_MGR: invalid grid level'
    ! The top loop for the coarsest grid is different from all others
    do ii=1,NSx/P + 1; do jj=1,NSy/P+1  ! +1 so a last odd one gets included too
       r = horizontaldistance1(P*ii*dx,P*jj*dy,x0,y0)
       if (r>=P/2*RMG) then  ! do 1 coarse grid cell
          select case (L)
          case (10)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h10,P)
          case (9)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h9,P)
          case (8)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h8,P)
          case (7)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h7,P)
          case (6)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h6,P)
          case (5)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h5,P)
          case (4)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h4,P)
          case (3)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h3,P)
          case (2)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h2,P)
          end select
          if (VERBOSE) write(6,*) 'Level',L,P*ii,P*jj,smax(90)
       else
          call findallhorizon_recursive(ii,jj,h,x0,y0,h00,naz,smax,RMG,L-1)
       endif
    enddo; enddo
  end subroutine findallhorizon_MGR


  pure recursive subroutine findallhorizon_recursive(i,j,h,x0,y0,h00,naz,smax,RMG,L)
    use allinterfaces, only: horizontaldistance1, horizon_MG_core
    implicit none
    integer, intent(IN) :: naz,L,i,j
    real(8), intent(IN) :: h(NSx,NSy)
    real(8), intent(IN) :: x0,y0,h00,RMG
    real(8), intent(INOUT) :: smax(naz)
    integer ii,jj,P
    real(8) r
    
    ! start with the coarsest grid
    ! if distance is too close, switch to finer grid
    
    P = 2**(L-1)
    if (L>9) error stop 'findallhorizon_recursive: exceeds maximum number of levels'
    if (L<1) error stop 'findallhorizon_recursive: level is zero or negative'
    do ii=2*i-1,2*i; do jj=2*j-1,2*j
       r = horizontaldistance1(P*ii*dx,P*jj*dy,x0,y0)
       if (P==1) r = 0.  ! should be redundant
       if (r>=P/2*RMG) then  ! do 1 coarse grid cell
          select case (L)
          case (9)  ! one below the highest
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h9,256)
          case (8) 
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h8,128)
          case (7)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h7,64)
          case (6)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h6,32)
          case (5)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h5,16)
          case (4)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h4,8)
          case (3)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h3,4)
          case (2)
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h2,2)
          case (1)
             if (ii<=1 .or. jj<=1 .or. ii>=NSx .or. jj>=NSy) cycle
             call horizon_MG_core(x0,y0,h00,naz,smax,ii,jj,h,1)
          end select
       else ! do 4 finer cells
          call findallhorizon_recursive(ii,jj,h,x0,y0,h00,naz,smax,RMG,L-1)
       endif
    enddo; enddo
  end subroutine findallhorizon_recursive

end module newmultigrid
