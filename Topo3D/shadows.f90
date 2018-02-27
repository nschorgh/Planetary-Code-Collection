program toposhadows
!***********************************************************************
! calculates horizon for every location and every azimuth
! written by Norbert Schorghofer 2010-2018
!***********************************************************************
  use filemanager, only : NSx,NSy,fileext,dx,dy,RMAX
  use allinterfaces
  use newmultigrid
  implicit none
  real(8), parameter :: pi=3.1415926535897932
  integer i, j, ext, narg, LACT
  integer, parameter :: naz=180      ! # of azimuths
  real(8) h(NSx,NSy), smax(naz)
  character(5) extc
  logical :: MULTIGRID = .false.
  real(8) RMG   ! distance of finest multigrid transition
  
  narg = COMMAND_ARGUMENT_COUNT()
  print *,'narg=',narg
  
  ! azimuth in degrees east of north, 0=north facing, 0...2*pi

  call readdem(h)
  print *,'...finished reading topography... ',fileext
  print *,'# domain size =',NSx*dx,'x',NSy*dy
  print *,'# azimuth rays = ',naz
  print *,'# fully sampled radius =',min(dx,dy)*naz/(2*pi)

  if (MULTIGRID) then
     !RMG = 50.
     RMG = naz*min(dx,dy)/(2*pi)
     print *,'# multgrid transition radius RMG =',RMG
     print *,'# grid levels =',ceiling(log(max(NSx*dx,NSy*dy)/RMG)/log(2.))

     call downsample_all(h,5,LACT)
     print *,'# levels allocated = ',LACT
  else
     print *,'# cutoff radius RMAX =',RMAX
  endif
  
  if (narg==0) then  ! serial implementation
     print *,'...creating file horizons.dat'
     open(unit=21,file='horizons.dat',status='unknown',action='write')

     do i=2,NSx-1
        print *,i
        do j=2,NSy-1
           if (.not.MULTIGRID) then
              call findallhorizon(h,i,j,naz,smax)
           else
              call findallhorizon_MGR(h,i,j,naz,smax,RMG,5)
           endif
           !write(21,'(2(i4,1x),9999(1x,f6.4))') i,j,smax(:)
           write(21,'(2(i4,1x))',advance='no') i,j
           call compactoutput(21,smax,naz)

        enddo
     enddo

  else  ! parallel implementation
     call getarg(1,extc)
     read(extc,'(i4)') ext  ! string->integer
     if (ext<=1 .or. ext>=NSx) stop 'argument is out of bounds'

     print *,'...creating file horizon....'
     open(unit=21,file='horizon.'//extc,status='unknown',action='write')
     
     i = ext  ! replaces loop over i=2,...,NSx-1
     print *,i
     do j=2,NSy-1
        call findallhorizon(h,i,j,naz,smax)
        !write(21,'(2(i4,1x),9999(1x,f6.4))') i,j,smax(:)
        write(21,'(2(i4,1x))',advance='no') i,j
        call compactoutput(21,smax,naz)
     enddo

  end if
  close(21)
end program toposhadows
