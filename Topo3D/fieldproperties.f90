program fieldproperties
!***********************************************************************
!   read horizons file and, optionally, field of views file
!   and calculate skysizes from them
!***********************************************************************
  use filemanager
  use allinterfaces
  use newhorizons

  implicit none
  integer i, j, CCMAX
  integer, dimension(NSx,NSy) :: cc
  real(8), dimension(NSx,NSy) :: h, skysize1, landsize, skysize2 ! skysize+landsize=2*pi
  !real(8), dimension(NSx,NSy) :: gterm, surfaceSlope, azFac
  integer(2), dimension(:,:,:), allocatable :: ii,jj
  real(4), dimension(:,:,:), allocatable :: dO12
  real(8), parameter :: pi=3.1415926535897932
  logical, parameter :: fieldofview=.false.
  character(len = *), parameter :: ffn='fieldofviews.'//fileext
  
  call readdem(h)
  !call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)
  
  print *,'...reading horizons file...'
  call readhorizons
  do i=2,NSx-1
     do j=2,NSy-1
        skysize1(i,j) = getoneskysize(i,j)
        skysize2(i,j) = getoneskysize_v2(i,j)
        !gterm(i,j ) = getoneGterm(i,j,surfaceSlope(i,j),azFac(i,j))
     enddo
  enddo
  
  if (fieldofview) then
     print *,'...reading huge fieldofviews file...'
     CCMAX = getmaxfieldsize(NSx,NSy,ffn)
     print *,'... max field of view size=',CCMAX
     allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), dO12(NSx,NSy,CCMAX))
     call getfieldofview(NSx,NSy,ffn,cc,ii,jj,dO12,landsize,CCMAX)
     !call getviewfactors(NSx,NSy,vfn,cc,ii,jj,dO12,landsize,CCMAX)
  else
     landsize = -9.
  end if

  print *,'...writing sky view factors...'
  open(unit=21,file='tmp.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,1x,f5.3,2x,f6.3)') &
             & i,j,h(i,j),2*pi-skysize1(i,j),landsize(i,j)
        !write(21,'(2(i4,1x),f9.2,1x,f6.4,1x,f6.3,3(2x,f6.3))') &
        !     & i,j,h(i,j),surfaceSlope(i,j),azFac(i,j),2*pi-skysize2(i,j),gterm(i,j)
     enddo
  enddo
  close(21)

end program fieldproperties

