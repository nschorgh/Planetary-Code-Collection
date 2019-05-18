program toposhadows
!***********************************************************************
! calculates horizon for every location and every azimuth
!
! written by Norbert Schorghofer 2010-2018
!***********************************************************************
  use filemanager, only : NSx,NSy,fileext,dx,dy,RMAX
  use allinterfaces
  use newmultigrid
  implicit none
  real(8), parameter :: pi=3.1415926535897932
  integer i, j, narg, ilower, iupper, LACT, LMAX
  integer, parameter :: naz=180      ! # of azimuths
  real(8) h(NSx,NSy), smax(naz)
  character(5) extc
  logical :: MULTIGRID = .true.
  real(8) RMG   ! distance of finest multigrid transition
  
  narg = COMMAND_ARGUMENT_COUNT()
  print *,'narg=',narg
  
  ! azimuth in degrees east of north, 0=north facing, 0...2*pi

  call readdem(h)
  print *,'...finished reading topography... ',fileext
  print *,'# domain size (pixels)',NSx,'x',NSy
  print *,'# domain size (m)',NSx*dx,'x',NSy*dy
  print *,'# azimuth rays = ',naz
  print *,'# fully sampled radius =',min(dx,dy)*naz/(2*pi)

  if (MULTIGRID) then
     RMG = naz*min(dx,dy)/(2*pi)
     print *,'# multgrid transition radius RMG =',RMG
     print *,'# grid levels =',ceiling(log(max(NSx*dx,NSy*dy)/RMG)/log(2.))
     LMAX = floor(log(sqrt((NSx*dx)**2+(NSy*dy)**2)/RMG)/log(2.))
     print *,'# log2(domain size/RMG) =',LMAX
     LMAX = min(10,LMAX)
     call downsample_all(h,LMAX,LACT)
     LMAX = min(LACT,LMAX)
     print *,'# levels allocated = ',LMAX
     
  else
     print *,'# cutoff radius RMAX =',RMAX
  endif

  select case(narg)
  case (0)  ! serial implementation
     print *,'...creating file horizons.dat'
     open(unit=21,file='horizons.dat',status='unknown',action='write')
     ilower = 2; iupper = NSx-1

  case (2)  ! parallel implementation
     call slicer(NSx,ilower,iupper,extc)
     print *,'...creating file horizon....'
     open(unit=21,file='horizon.'//extc,status='unknown',action='write')

  case default
     stop 'Number of command line arguments must be zero or two' 
  end select

  do i=ilower,iupper
     print *,i
     do j=2,NSy-1
        if (.not.MULTIGRID) then
           call findallhorizon(h,i,j,naz,smax)
        else
           call findallhorizon_MGR(h,i,j,naz,smax,RMG,LMAX)
        endif
        write(21,'(2(i0,1x))',advance='no') i,j
        call compactoutput(21,smax,naz)
     enddo
  enddo
  
  close(21)
end program toposhadows
