program fieldproperties
!***********************************************************************
!   read horizons file and, optionally, field of views file
!   and calculate skysizes from them
!***********************************************************************
  use filemanager
  use allinterfaces
  !use newhorizons

  implicit none
  integer i, j, CCMAX
  real(8), dimension(NSx,NSy) :: h, skysize, landsize ! skysize+landsize=2*pi
  integer, dimension(NSx,NSy) :: cc
  integer(2), dimension(:,:,:), allocatable :: ii,jj
  real(4), dimension(:,:,:), allocatable :: dO12
  real(8), parameter :: pi=3.1415926535897932
  logical, parameter :: fieldofview=.false.

  call readdem(h)

  print *,'...reading horizons file...'
  !call readhorizons(skysize,'horizons.'//fileext)
  !call getskysize(skysize)
  call getskysize(skysize,'horizons.'//fileext)
  
  if (fieldofview) then
     print *,'...reading huge fieldofviews file...'
     call getmaxfieldsize(NSx,NSy,fileext,CCMAX)
     print *,'... max field of view size=',CCMAX
     allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), dO12(NSx,NSy,CCMAX))
     call getfieldofview(NSx,NSy,fileext,cc,ii,jj,dO12,landsize,CCMAX)
  else
     landsize = -9.
  end if

  print *,'...writing sky view factors...'
  open(unit=21,file='tmp.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,1x,f5.3,2x,f6.3)') &
             & i,j,h(i,j),2*pi-skysize(i,j),landsize(i,j)
     enddo
  enddo
  close(21)

end program fieldproperties


subroutine getskysize(skysize,fn)
!***********************************************************************
!   reads horizons file and calculates sky size (steradian)
!   there is an abbreviationed version in module newhorizons
!***********************************************************************
  use filemanager, only : NSx,NSy,fileext
  use allinterfaces
  implicit none
  character(len=*), intent(IN) :: fn
  real(8), intent(OUT) :: skysize(NSx,NSy) 
  real(8), parameter :: pi=3.1415926535897932
  integer, parameter :: naz=180*2   ! # of azimuths
  real(8) smax(naz)
  integer i, j, ii, jj, ierr, k, kp1
  real(8) phi(3), theta(3), dOmega, landsize0

  ! azimuth in degrees east of north, 0=north facing, 0...2*pi

  print *,'# azimuth rays = ',naz
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  
  print *,'...reading horizons file ...'
  open(unit=21,file=fn,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getskysize: Input file not found'
  
  do i=2,NSx-1
     do j=2,NSy-1
        read(21,*) ii,jj,smax(:)
        if (ii/=i .or. jj/=j) stop 'getskysize: index mismatch'

        ! approximate
        landsize0 = sum(atan(smax))*2*pi/naz

        ! exact
        skysize(i,j)=0
        phi = (/ 0, 0, 1 /) *2*pi/naz
        do k=1,naz
           !kp1 = k+1; if (kp1>naz) kp1=kp1-naz
           kp1 = mod(k,naz)+1
           theta = (/ 0.d0, atan(1/smax(k)), atan(1/smax(kp1)) /) ! from zenith
           dOmega = area_spherical_triangle(phi,theta)
           skysize(i,j) = skysize(i,j) + dOmega
        enddo

        !print *,i,j,2*pi-landsize0,skysize(i,j)
     enddo
  enddo

  close(21)
end subroutine getskysize
