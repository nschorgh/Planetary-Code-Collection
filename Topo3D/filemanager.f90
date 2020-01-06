MODULE filemanager 
  implicit none
  ! NSx = # pixels in x-direction (longitude)
  ! NSy = # pixels in y-direction (latitude)
  ! dx, dy = horizontal resolution [meter]
  ! RMAX = max radius of consideration [meter]

!-begin idealized topographies
  integer, parameter :: NSx=40, NSy=40
  character(len=*), parameter :: fileext = 'topo40'
  !integer, parameter :: NSx=81, NSy=81
  !character(len=*), parameter :: fileext = 'topo81'
  !integer, parameter :: NSx=81, NSy=71
  !character(len=*), parameter :: fileext = 'boulder2'
  real(8), parameter :: dx=5., dy=5., RMAX=1e35

  
!-begin Mauna Kea topographies
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
  
  !integer, parameter :: NSx=950, NSy=700  
  !character(len = *), parameter :: fileext = 'summitregion3'
  !real(8), parameter :: dx=4.4971, dy=4.4971, RMAX=5e3
  
  
!-begin Moon topographies
  !integer, parameter :: NSx=128, NSy=128
  !character(len = *), parameter :: fileext = 'moongauss3.rms=0.5'
  !real(8), parameter :: dx=1., dy=1., RMAX=1e5

  
!-begin Mars topographies
  !integer, parameter :: NSx=112, NSy=292
  !character(len = *), parameter :: fileext = 'PSP005943_64'
  !real(8), parameter :: dx=64., dy=64., RMAX=1e30

  !integer, parameter :: NSx=224, NSy=585
  !character(len = *), parameter :: fileext = 'PSP005943_32'
  !real(8), parameter :: dx=32., dy=32., RMAX=1e30

  !integer, parameter :: NSx=448, NSy=1170
  !character(len = *), parameter :: fileext = 'PSP005943_16'
  !real(8), parameter :: dx=16., dy=16., RMAX=1e30

  !integer, parameter :: NSx=896, NSy=2340
  !character(len = *), parameter :: fileext = 'PSP005943_8'
  !real(8), parameter :: dx=8., dy=8., RMAX=1e30
  
  !integer, parameter :: NSx=1792, NSy=4680
  !character(len = *), parameter :: fileext = 'PSP005943_4'
  !real(8), parameter :: dx=4., dy=4., RMAX=1e30

  !integer, parameter :: NSx=3584, NSy=9360
  !character(len = *), parameter :: fileext = 'PSP005943_2'
  !real(8), parameter :: dx=2., dy=2., RMAX=1e30

  !integer, parameter :: NSx=7169, NSy=18720
  !character(len = *), parameter :: fileext = 'PSP005943_1'
  !real(8), parameter :: dx=1., dy=1., RMAX=1e30

  !integer, parameter :: NSx=1050, NSy=1200
  !character(len = *), parameter :: fileext = 'ctxwider'
  !real(8), parameter :: dx=18., dy=18., RMAX=1e30

  !integer, parameter :: NSx=1050/2, NSy=1200/2
  !character(len = *), parameter :: fileext = 'ctxwider_36m'
  !real(8), parameter :: dx=36., dy=36., RMAX=1e30

  
  
  ! paths and default filenames for inputs
  character(len = *), parameter :: &
       hfn='./'//fileext//'.xyz' ,& ! path and filename for topography (*.xyz)
       sfn='./horizons.'//fileext ,& ! path and filename for horizons.* input 
       vfn='./viewfactors.'//fileext  ! path and filename for viewfactors.* input
  ! outputs are written to working directory
END MODULE filemanager

