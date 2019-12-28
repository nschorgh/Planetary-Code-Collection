module azRays
  implicit none
  integer, parameter :: naz = 180  ! # of azimuths
  real(8), parameter :: pi = 3.1415926535897931
  real(8), parameter :: f = naz/(2*pi)
  
  integer, private :: ak
  real(8), parameter :: azRay(naz) = (/ ( (ak-1)/f, ak=1,naz) /)
  ! inverse mapping:  ak = azRay*f+1
  ! pre-defining azRay leads to performance improvement in horizon_core
end module azRays



subroutine findallhorizon1(h,i0,j0,naz,smax)
  ! find all horizon heights, without use of multigrid
  use filemanager, only : NSx,NSy,dx,dy,RMAX
  use allinterfaces, only : horizontaldistance1
  implicit none
  integer, intent(IN) :: i0,j0,naz
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT) :: smax(naz)
  integer i,j
  real(8) x0,y0,h00
  
  smax(:)=0.
  x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)
  do i=2,NSx-1
     !if (horizontaldistance(i,1,i0,1)>RMAX) cycle  ! saves computations
     if (horizontaldistance1(i*dx,1*dy,x0,1*dy)>RMAX) cycle  ! saves computations
     do j=2,NSy-1
        !if (horizontaldistance(i,j,i0,j0)>RMAX) cycle  ! saves computations
        if (horizontaldistance1(i*dx,j*dy,x0,y0)>RMAX) cycle  ! saves computations
        call horizon_core(x0,y0,h00,smax,i,j,h,1)
     end do
  end do

end subroutine findallhorizon1



pure subroutine horizon_core(x0,y0,h00,smax,i,j,h,P)
  ! can be used for single- or multi-grid
  use filemanager, only : NSx, NSy, dx, dy
  use allinterfaces, only : horizontaldistance1, azimuth1, diffangle
  use azRays
  implicit none
  real(8), intent(IN) :: x0, y0, h00  ! on fine grid
  real(8), intent(INOUT) :: smax(naz)
  integer, intent(IN) :: i, j, P  ! fine or coarse grid
  real(8), intent(IN) :: h(NSx/P,NSy/P) ! fine or coarse grid

  integer k,in,jn,ak,ak1,ak2,buf,akak
  integer nx, ny    ! grid-level specific grid size
  real(8) az,az_neighbor,t,r,r_neighbor,rcut,s,hcut,d1,d2,d3
  real(8) dxl, dyl  ! grid-level specific resolution
  integer, parameter :: ex(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)
  integer, parameter :: ey(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)
  !real(8) f,azRay(naz)
  !real(8), parameter :: pi=3.1415926535897931

  nx = NSx/P; ny = NSy/P
  dxl = P*dx; dyl = P*dy
  ! no grid-level specific variables are used below

  if (i>nx .or. j>ny) return

  r = horizontaldistance1(i*dxl,j*dyl,x0,y0)
  if (r==0.) return
  az = azimuth1(x0,y0,i*dxl,j*dyl)

  if (floor(az*f)==ceiling(az*f)) then  ! grid point lies on ray
     ak = nint(az*f)+1
     s = (h(i,j)-h00) / r
     if (s>smax(ak)) smax(ak)=s
     return
  endif
  
  do k=1,8
     in = i+ex(k)
     jn = j+ey(k)
     if (in<1 .or. in>nx .or. jn<1 .or. jn>ny) cycle
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
     if (ak1>naz .or. ak2>naz) error stop 'horizon_core: index out of bound'

     d3=diffangle(az,az_neighbor)
     do akak=ak1,ak2
        ak = akak; if (ak<=0) ak = ak+naz
        d1=diffangle(az,azRay(ak))
        d2=diffangle(az_neighbor,azRay(ak))
        
        if (d1+d2 <= d3+1.d-5) then  
           if (d1>1.0*d3 .and. d3>1.d-6) cycle 
           r_neighbor = horizontaldistance1(in*dxl,jn*dyl,x0,y0)
           ! edge between h1,i0,j0 and h2,in,jn
           if (d3>1.d-6) then
              t = d1/d3  ! approximation
           else
              t = 0.5  ! dirty fix
           endif
           hcut = h(i,j)*(1-t) + h(in,jn)*t
           rcut = r*(1-t) + r_neighbor*t
           s = (hcut-h00) / rcut
           if (s>smax(ak)) smax(ak)=s
        endif
     end do  ! end of ak loop
  end do  ! end of k loop

end subroutine horizon_core

