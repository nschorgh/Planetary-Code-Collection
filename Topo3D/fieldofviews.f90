program fieldofviews
!*****************************************************************************
! 1. calculates horizon for every location and every azimuth
! 2. determines visibility of surface elements
! 3. calculates subtended angle for every surface element
!
! Without input parameter, it executes serially
! With integer as input parameter, it executes parallel implementation
! For example, a.out 2    executes second row
!
! written by Norbert Schorghofer 2010-2015  
!*****************************************************************************
  use filemanager, only : NSx,NSy
  use allinterfaces
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer i, j, k, ext, narg
  integer, parameter :: nres=180   ! # of azimuths
  real(8) h(NSx,NSy), azSun, smax(nres)
  character(5) extc
  logical visibility(NSx,NSy)

  ! azimuth in degrees east of north, 0=north facing, 0...2*pi

  narg = COMMAND_ARGUMENT_COUNT()
  print *,'narg=',narg
  
  call readdem(h)
  print *,'...finished reading topography...'

  if (narg==0) then  ! serial implementation
     print *,'...creating files horizons.dat and fieldofviews.dat'
     open(unit=21,file='horizons.dat',status='unknown',action='write')
     open(unit=22,file='fieldofviews.dat',status='unknown',action='write')
     print *,'Nsx=',nsx,'Nsy=',nsy,'# azimuths=',nres
     do i=2,NSx-1
        print *,i
        do j=2,NSy-1
           visibility(:,:) = .false.
           do k=1,nres
              azSun = 360./real(nres)*(k-1)*d2r
              call findonehorizon_wsort(h,i,j,azSun,smax(k),visibility)
           enddo
           !write(21,'(2(i4,1x),9999(1x,f6.4))') i,j,smax(:)
           write(21,'(2(i4,1x))',advance='no') i,j
           call compactoutput(21,smax,nres)
           call refinevisibility(i,j,h,visibility)
           call find3dangle(h,i,j,22,visibility)
        enddo
     enddo

  else  ! parallel implementation
     call getarg(1,extc)
     read(extc,'(i4)') ext  ! string->integer
     if (ext<=1 .or. ext>=NSx) stop 'argument is out of bounds'  
     i = ext  ! replaces loop over i=2,...,NSx-1
     print *,'...creating file horizon and fieldofview...',i
     open(unit=21,file='horizon.'//extc,status='unknown',action='write')
     open(unit=22,file='fieldofview.'//extc,status='unknown',action='write')
     do j=2,NSy-1
        visibility(:,:) = .false.
        do k=1,nres
           azSun = (360./real(nres))*(k-1)*d2r
           call findonehorizon_wsort(h,i,j,azSun,smax(k),visibility)
        enddo
        write(21,'(2(i4,1x))',advance='no') i,j
        call compactoutput(21,smax,nres)
        call refinevisibility(i,j,h,visibility)
        call find3dangle(h,i,j,22,visibility)
     enddo
     ! merge and sort output files with script
  endif
  
  close(21)
  close(22)
  print *,'Finished'
end program fieldofviews


