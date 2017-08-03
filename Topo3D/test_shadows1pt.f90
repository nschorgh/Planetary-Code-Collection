program test_shadows1pt
  ! test shadows for one point
  use filemanager, only : NSx,NSy,fileext,dx,dy,RMAX
  use allinterfaces
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer i, j, k
  integer, parameter :: nres=180   ! # of azimuths
  real(8) h(NSx,NSy), azSun, smax(nres)
  !real(8) RMG,h2(NSx/2,NSy/2),h3(NSx/4,NSy/4)  ! multigrid
  !logical visibility(NSx,NSy) 

  ! azimuth in degrees east of north, 0=north facing, 0...2*pi
  
  call readdem(h)
  print *,'...finished reading topography... ',fileext
  print *,'# domain size =',NSx*dx,'x',NSy*dy
  print *,'# azimuth rays = ',nres
  print *,'# fully sampled radius =',min(dx,dy)*nres/(2*pi)
  
  ! downsample grid
  !call downsample(NSx,NSy,h,h2)
  !call downsample(NSx/2,NSy/2,h2,h3)
  
  print *,'...creating file horizons.dat...'
  open(unit=31,file='horizons.dat',status='unknown',action='write')

  !RMG = 50.   ! RMG*2*pi = nres*max(dx,dy)
  !RMG = nres*max(dx,dy)/(2*pi)
  !print *,'domain size=',NSx*dx,'x',NSy*dy,'RMG=',RMG

  i=17; j=31
  print *,i,j
  call findallhorizon(h,i,j,nres,smax(:))
  !call findallhorizon_MG1(h,i,j,nres,smax(:))
  !call findallhorizon_MG3(h,h2,h3,i,j,nres,smax,RMG)
  do k=1,nres
     azSun = 360./real(nres)*(k-1)*d2r
     !call findonehorizon(h,i,j,azSun,smax(k))
     !call findonehorizon_wsort(h,i,j,azSun,smax(k),visibility)
     write(31,'(f4.0,2x,f6.4)') azSun/d2r,smax(k)
  end do

  close(31)
end program test_shadows1pt
