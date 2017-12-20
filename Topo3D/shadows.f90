program toposhadows
!***********************************************************************
! calculates horizon for every location and every azimuth
! written by Norbert Schorghofer 2010-2017
!***********************************************************************
  use filemanager, only : NSx,NSy,fileext,dx,dy,RMAX
  use allinterfaces
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer i, j, k, ext, narg
  integer, parameter :: naz=180      ! # of azimuths
  real(8) h(NSx,NSy), azSun, smax(naz)
  character(5) extc
  logical :: MULTIGRID = .false.
  real(8) RMG   ! distance of finest multigrid transition
  real(8), allocatable :: h2(:,:), h3(:,:), h4(:,:), h5(:,:)  ! multigrid

  
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

     allocate(h2(NSx/2,NSy/2)); call downsample(NSx,NSy,h,h2)  
     if (min(NSx,NSy)/4>10) then
        allocate(h3(NSx/4,NSy/4)); call downsample(NSx/2,NSy/2,h2,h3)
        if (min(NSx,NSy)/8>10) then
           allocate(h4(NSx/8,NSy/8)); call downsample(NSx/4,NSy/4,h3,h4)
           if (min(NSx,NSy)/16>10) then
              allocate(h5(NSx/16,NSy/16)); call downsample(NSx/8,NSy/8,h4,h5)
           endif
        endif
     endif
  else
     print *,'# cutoff radius RMAX =',RMAX
  endif
  
  if (narg==0) then  ! serial implementation
     print *,'...creating file horizons.dat...'
     open(unit=21,file='horizons.dat',status='unknown',action='write')

     do i=2,NSx-1
        print *,i
        do j=2,NSy-1
           if (.not.MULTIGRID) then
              call findallhorizon(h,i,j,naz,smax)
           else
              call findallhorizon_MG3(h,h2,h3,i,j,naz,smax,RMG)
           endif
           do k=1,naz
              azSun = 360./real(naz)*(k-1)*d2r
              !write(32,'(2(i4,1x),1x,f4.0,2x,f6.4)') i,j,azSun/d2r,smax(k)
              !call findonehorizon(h,i,j,azSun,smax(k))
              !write(31,'(2(i4,1x),1x,f4.0,2x,f6.4)') i,j,azSun/d2r,smax(k)
              !call findonehorizon_wsort(h,i,j,azSun,smax(k),visibility)
           end do
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
