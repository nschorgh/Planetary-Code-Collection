program maketopo
  ! create artifical topography 
  implicit none
  !integer, parameter :: NSx=30, NSy=20
  integer, parameter :: NSx=81, NSy=81
  integer i,j,i0,j0
  real(8) h(NSx,NSy), dx, x, y
  real(8) depthtodiam, R, z0
  
  dx = 5.  ! horizontal resolution

  depthtodiam = 1./5.
  R= 50*dx
  i0=ceiling(NSx/2.); j0=ceiling(NSy/2.)

  !depthtodiam = 0.5
  !R=15*dx
  !i0=NSx/2; j0=NSx*4./5.

  ! crater depth  
  z0 = R*(1-4*depthtodiam**2)/(1+4*depthtodiam**2)
  ! crater diameter
  !D = z0/depthtodiam;  
    
  open(unit=20,file='topo.xyz',status='unknown',action='write')
  do j=1,NSy
     do i=1,NSx
        x = dx*(i-i0)
        y = dx*(j-j0)

        ! spherical crater
        if (x**2+y**2<R**2) then
           h(i,j) = z0-sqrt(R**2-x**2-y**2)
        else
           h(i,j) = 0.
        endif
        if (h(i,j)>0) h(i,j)=0.

     enddo
     write(20,'(999(f7.3,1x))') h(1:NSx,j)
  enddo
end program maketopo

