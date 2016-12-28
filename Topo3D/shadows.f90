program toposhadows
!***********************************************************************
! calculates horizon for every location and every azimuth
! written by Norbert Schorghofer 2010-2015
!***********************************************************************
  use filemanager, only : NSx,NSy,fileext,dx,dy,RMAX
  use allinterfaces
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer i, j, k, ext, narg
  integer, parameter :: nres=360   ! # of azimuths
  real(8) h(NSx,NSy), azSun, smax(nres)
  character(5) extc

  narg = COMMAND_ARGUMENT_COUNT()
  print *,'narg=',narg
  
  ! azimuth in degrees east of north, 0=north facing, 0...2*pi

  call readdem(h)
  print *,'...finished reading topography... ',fileext
  print *,'# domain size =',NSx*dx,'x',NSy*dy
  print *,'# azimuth rays = ',nres
  print *,'# fully sampled radius =',min(dx,dy)*nres/(2*pi)
  print *,'# cutoff radius RMAX =',RMAX  

  if (narg==0) then  ! serial implementation
     print *,'...creating file horizons.dat...'
     open(unit=21,file='horizons.dat',status='unknown',action='write')
     
     do i=2,NSx-1
        print *,i
        do j=2,NSy-1
           do k=1,nres
              azSun = 360./real(nres)*(k-1)*d2r
              call findonehorizon(h,i,j,azSun,smax(k))
              !write(31,'(2(i4,1x),1x,f4.0,2x,f6.4)') i,j,azSun/d2r,smax(k)
           end do
           !write(21,'(2(i4,1x),9999(1x,f6.4))') i,j,smax(:)
           write(21,'(2(i4,1x))',advance='no') i,j
           call compactoutput(21,smax,nres)
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
        do k=1,nres
           azSun = 360./real(nres)*(k-1)*d2r
           call findonehorizon(h,i,j,azSun,smax(k))
        enddo
        !write(21,'(2(i4,1x),9999(1x,f6.4))') i,j,smax(:)
        write(21,'(2(i4,1x))',advance='no') i,j
        call compactoutput(21,smax,nres)
     enddo

  end if
  close(21)
end program toposhadows
