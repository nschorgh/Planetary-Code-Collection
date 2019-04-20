module filemanager 
  implicit none
  ! NSx = # pixels in x-direction (longitude)
  ! NSy = # pixels in y-direction (latitude)
  ! dx, dy = horizontal resolution [meter]
  ! RMAX = max radius of consideration [meter]

  !integer, parameter :: NSx=30, NSy=20
  !character(len=*), parameter :: fileext = 'topo20'
  !real(8), parameter :: dx=5., dy=5., RMAX=1e35

  integer, parameter :: NSx=40, NSy=40
  character(len=*), parameter :: fileext = 'topo40'
  !integer, parameter :: NSx=81, NSy=81
  !character(len=*), parameter :: fileext = 'topo81'
  !character(len=*), parameter :: fileext = 'boulder'
  real(8), parameter :: dx=5., dy=5., RMAX=1e35

  !integer, parameter :: NSx=131, NSy=191 
  !character(len = *), parameter :: fileext = 'summit_large'
  !real(8), parameter :: dx=10., dy=10., RMAX=2e3
  
  !integer, parameter :: NSx=71, NSy=81 
  !character(len = *), parameter :: fileext = 'summit_small'
  !real(8), parameter :: dx=10., dy=10., RMAX=2e3

  !integer, parameter :: NSx=800, NSy=700  
  !character(len = *), parameter :: fileext = 'summitregion'
  !real(8), parameter :: dx=4.4971, dy=4.4971, RMAX=2e3

  !integer, parameter :: NSx=400, NSy=350  
  !character(len = *), parameter :: fileext = 'summitregion_9m'
  !real(8), parameter :: dx=8.9941, dy=8.9941, RMAX=2e3
  
  !integer, parameter :: NSx=950, NSy=700  
  !character(len = *), parameter :: fileext = 'summitregion3'
  !real(8), parameter :: dx=4.4971, dy=4.4971, RMAX=4e3

  !integer, parameter :: NSx=150, NSy=180 
  !character(len = *), parameter :: fileext = 'rsl4'
  !real(8), parameter :: dx=4., dy=4., RMAX=1e3
  !integer, parameter :: MARGIN = 20
  !character(len = *), parameter :: fileext = 'rsl1'
  !real(8), parameter :: dx=1., dy=1., RMAX=1e3
  !integer, parameter :: MARGIN = 20*4
  
  !integer, parameter :: NSx=128, NSy=128
  !character(len = *), parameter :: fileext = 'moongauss3.rms=0.5'
  !real(8), parameter :: dx=1., dy=1., RMAX=1e5

  !integer, parameter :: NSx=1792, NSy=4680
  !character(len = *), parameter :: fileext = 'PSP005943_4'
  !real(8), parameter :: dx=4., dy=4., RMAX=1e30

  !integer, parameter :: NSx=896, NSy=2340
  !character(len = *), parameter :: fileext = 'PSP005943_8'
  !real(8), parameter :: dx=8., dy=8., RMAX=1e30

  !integer, parameter :: NSx=448, NSy=1170
  !character(len = *), parameter :: fileext = 'PSP005943_16'
  !real(8), parameter :: dx=16., dy=16., RMAX=1e30

  !integer, parameter :: NSx=224, NSy=585
  !character(len = *), parameter :: fileext = 'PSP005943_32'
  !real(8), parameter :: dx=32., dy=32., RMAX=1e30

  
  ! paths and default filenames for inputs
  character(len = *), parameter :: &
       hfn='./'//fileext//'.xyz' ,& ! path and filename for topography (*.xyz)
       sfn='./horizons.'//fileext ,& ! path and filename for horizons.* input 
       ffn='./fieldofviews.'//fileext  ! path and filename for fieldofviews.* input
  ! outputs are written to working directory
end module filemanager


