module azRays
  implicit none
  integer, parameter :: naz=180  ! # of azimuths
  real(8), parameter :: pi=3.1415926535897931
  real(8), parameter :: f=naz/(2*pi)
  
  integer, private :: ak
  real(8), parameter :: azRay(naz) = (/ ( (ak-1)/f, ak=1,naz) /)
  ! inverse mapping:  ak = azRay*f+1
  ! pre-defining azRay leads to performance improvement in horizon_core
end module azRays
