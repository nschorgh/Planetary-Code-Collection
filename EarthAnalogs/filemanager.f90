
! (1,1) = northwest corner
! (NSx,1) = northeast corner
! (1,NSy) = southwest corner
! (NSx,NSy) = southeast corner

module filemanager 
  implicit none
  ! NSx = # pixels in x-direction (longitude)
  ! NSy = # pixels in y-direction (latitude)
  ! dx, dy = horizontal resolution [meter]
  ! RMAX = max radius of consideration [meter]

  !integer, parameter :: NSx=40, NSy=40
  !character(len=*), parameter :: fileext = 'topo40'
  !integer, parameter :: NSx=81, NSy=81
  !character(len=*), parameter :: fileext = 'topo81'
  !real(8), parameter :: dx=5., dy=5., RMAX=1e35

  
  !integer, parameter :: NSx=131, NSy=191 
  !character(len = *), parameter :: fileext = 'summit_large'
  !real(8), parameter :: dx=10., dy=10., RMAX=5e3
  
  !integer, parameter :: NSx=71, NSy=81 
  !character(len = *), parameter :: fileext = 'summit_small'
  !real(8), parameter :: dx=10., dy=10., RMAX=2e3

  !integer, parameter :: NSx=800, NSy=700  
  !character(len = *), parameter :: fileext = 'summitregion'
  !real(8), parameter :: dx=4.4971, dy=4.4971, RMAX=5e3

  !integer, parameter :: NSx=475, NSy=350
  !character(len = *), parameter :: fileext = 'summitregion3_9m'
  !real(8), parameter :: dx=8.9941, dy=8.9941, RMAX=5e3

  integer, parameter :: NSx=950, NSy=700  
  character(len = *), parameter :: fileext = 'summitregion3'
  real(8), parameter :: dx=4.4971, dy=4.4971, RMAX=5e3
  
  
  
  ! paths and default filenames for inputs
  character(len = *), parameter :: &
       hfn='./'//fileext//'.xyz' ,& ! path and filename for topography (*.xyz)
       sfn='./horizons.'//fileext  ! path and filename for horizons.* input 
  ! outputs are written to working directory
end module filemanager

