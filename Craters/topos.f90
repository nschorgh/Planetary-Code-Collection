
! (1,1) = northwest corner
! (NSx,1) = northeast corner
! (1,NSy) = southwest corner
! (NSx,NSy) = southeast corner

module filemanager 
  implicit none
  ! NSx = # pixels in x-direction
  ! NSy = # pixels in y-direction 
  ! dx, dy = horizontal resolution (meter)
  ! RMAX = max radius of consideration (meter)

  !integer, parameter :: NSx=30, NSy=20
  !character(len=*), parameter :: fileext = 'topo20'
  !real(8), parameter :: dx=5., dy=5., RMAX=1e35

  integer, parameter :: NSx=40, NSy=40
  character(len=*), parameter :: fileext = 'topo40'
  !integer, parameter :: NSx=81, NSy=81
  !character(len=*), parameter :: fileext = 'topo81'
  !integer, parameter :: NSx=80, NSy=80
  !character(len=*), parameter :: fileext = 'boulder'
  real(8), parameter :: dx=5., dy=5., RMAX=1e35

  !integer, parameter :: NSx=131, NSy=191 
  !character(len = *), parameter :: fileext = 'summit_large'
  !real(8), parameter :: dx=10., dy=10., RMAX=2e3
  
  !integer, parameter :: NSx=71, NSy=81  ! x=lon, y=lat
  !character(len = *), parameter :: fileext = 'summit_small'
  !real(8), parameter :: dx=10., dy=10., RMAX=2e3

  !integer, parameter :: NSx=901, NSy=801  
  !character(len = *), parameter :: fileext = 'summit_very'
  !real(8), parameter :: dx=10., dy=10., RMAX=2e3

  !integer, parameter :: NSx=800, NSy=700  ! x=lon, y=lat
  !character(len = *), parameter :: fileext = 'summitregion'
  !real(8), parameter :: dx=4.4971, dy=4.4971, RMAX=2e3

  !integer, parameter :: NSx=400, NSy=350  ! x=lon, y=lat
  !character(len = *), parameter :: fileext = 'summitregion_9m'
  !real(8), parameter :: dx=8.9941, dy=8.9941, RMAX=2e3

  !integer, parameter :: NSx=220, NSy=230  ! x=lon, y=lat
  !character(len = *), parameter :: fileext = 'summitregion_9m_smaller'
  !real(8), parameter :: dx=8.9941, dy=8.9941, RMAX=1e3

  logical, parameter :: verbose=.false.

  ! For subsurface model
  !real(8), parameter :: solarDay = 86400., Fgeotherm = 0.065  ! Earth
  !integer, parameter :: nz=40
  real(8), parameter :: solarDay = 29.53*86400., Fgeotherm = 0.029  ! The Moon
  integer, parameter :: nz=25
end module filemanager



module filemanager_ll
  implicit none
  ! NSx = # pixels in lon-direction
  ! Nsy = # pixels in lat-direction 
  ! RMAX = max radius of consideration (meter)
  ! PER = longitude coordinate periodic or not

  !integer, parameter :: NSx=1440, Nsy=60    ! whole
  integer, parameter :: NSx=720, Nsy=60    ! half
  !character(*), parameter :: fileext = '/Users/norbert/Moon/Herschel/topo4s_combined.txt'
  !character(*), parameter :: fileext = '/Users/norbert/Moon/Herschel/Data/topo4s_combined_2009-10-09.txt'
  character(*), parameter :: fileext = '/Users/norbert/Moon/Herschel/Data/topo4s_combined_2013-07-29.txt'

  !integer, parameter :: NSx=5760, Nsy=240  ! whole
  !integer, parameter :: NSx=2880, Nsy=240  ! half
  !character(len=*), parameter :: fileext = '/Users/norbert/Moon/Herschel/tmp16s_combined.txt'

  logical, parameter :: verbose=.FALSE.
  logical, parameter :: PER=.FALSE.
  real(8), parameter :: Rbody = 1737.1e3, RMAX=250e3  ! Moon
end module filemanager_ll



subroutine readdem(NSx,NSy,h,fileext)
  ! read DEM
  implicit none
  integer, intent(IN) :: NSx,NSy
  real(8), intent(OUT) :: h(NSx,NSy)
  character(len=*), intent(IN) :: fileext
  integer j, ierr
  open(unit=20,file=fileext//'.xyz',status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'readdem: input file not found'
  do j=1,NSy
     read(20,*) h(:,j)
  enddo
  close(20)
end subroutine readdem



!subroutine readdem_ll(NSx,NSy,fileext,lon,lat,h,az,el,NrAz)
subroutine readdem_ll(NSx,NSy,fileext,lon,lat,h)
  ! read elevations on lon-lat grid
  implicit none
  integer, intent(IN) :: NSx,NSy !,NrAz
  character*(*), intent(IN) :: fileext
  real(8), dimension(NSx,NSy), intent(OUT) :: h,lon,lat
  !real(8), dimension(NSx,NSy,NrAz), intent(OUT), optional :: az,el
  integer i, j, ierr
  open(unit=20,file=fileext,status='old',action='read',iostat=ierr)
  if (ierr>0) then
     print *,fileext
     stop 'readdem_ll: input file not found'
  endif
  do j=1,NSy
     do i=1,NSx
        read(20,*) lon(i,j),lat(i,j),h(i,j)
        !read(20,*) lon(i,j),lat(i,j),h(i,j),az(i,j,:),el(i,j,:)
     enddo
  enddo
  close(20)
end subroutine readdem_ll

