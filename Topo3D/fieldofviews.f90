PROGRAM fieldofviews
!*****************************************************************************
! 1. calculates horizon for every location and every azimuth
! 2. determines visibility of surface elements
! 3. calculates view factors for every surface element
!
! Without command line arguments, it executes serially
! With two integers as input arguments, it executes parallel implementation
!     e.g., 'a.out 6 2'  processes second out of six slices of domain
!
! written by Norbert Schorghofer 2010-2019
!*****************************************************************************
  use filemanager, only : NSx,NSy
  use allinterfaces
  use findviewfactors
  implicit none
  integer i, j, narg, ilower, iupper
  ! # of azimuths (naz) defined in module findviewfactors
  real(8) h(NSx,NSy), smax(naz)
  character(5) extc
  logical visibility(NSx,NSy)

  ! azimuth in degrees east of north, 0=north facing, 0...2*pi

  narg = COMMAND_ARGUMENT_COUNT()
  print *,'narg=',narg

  call readdem(h)
  print *,'...finished reading topography...'
  print *,'Nsx=',NSx,'Nsy=',NSy,'# azimuths=',naz

  select case(narg)
  case (0)  ! serial implementation
     print *,'...creating files horizons.dat and viewfactors.dat'
     !open(unit=21,file='horizons.dat',status='unknown',action='write')
     !open(unit=22,file='fieldofviews.dat',status='unknown',action='write')
     open(unit=23,file='viewfactors.dat',status='unknown',action='write')
     ilower = 2; iupper = NSx-1

  case (2)  ! parallel implementation
     call slicer(NSx,ilower,iupper,extc)
     print *,'...creating files horizon and viewfactor...',extc
     !open(unit=21,file='horizon.'//extc,status='unknown',action='write')
     !open(unit=22,file='fieldofview.'//extc,status='unknown',action='write')
     open(unit=23,file='viewfactor.'//extc,status='unknown',action='write')
     ! merge and sort output files with script
     
  case default
     stop 'Number of command line arguments must be zero or two' 
  end select
  
  if (ilower<=1 .or. iupper>=NSx) stop

  do i=ilower,iupper
     print *,i
     do j=2,NSy-1
        visibility(:,:) = .false.
        call findallhorizon_wsort_v3(h,i,j,smax,visibility)
        !write(21,'(2(i0,1x))',advance='no') i,j
        !call compactoutput(21,smax,naz)
        call refinevisibility(i,j,h,visibility)
        call find3dangle(h,i,j,23,visibility)
     enddo
  enddo

  !close(21)
  !close(22)
  close(23)
  print *,'Finished'

END PROGRAM fieldofviews


