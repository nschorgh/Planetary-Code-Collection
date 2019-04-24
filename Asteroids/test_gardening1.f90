program testgardening
  ! test impactstirring with constant time step and nothing else
  implicit none
  integer nz, n
  real(8) zfac, zmax, bigstep
  !parameter(nz=160, zfac=1.d0, zmax=0.5) 
  parameter(nz=160, zfac=1.05d0, zmax=10.) 
  real(8) z(nz), sigma(nz)
  real(8) m0, m1, m2
  real(8), external :: colint

  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z.dat',action='write',status='unknown')
  write(30,'(999(f8.5,1x))') z(1:nz)
  close(30)

  bigstep = 1e7  ! years

  sigma = 0.
  sigma(1) = 400./0.00021*1e-2

  !sigma = 400.
  !where(z<0.10) sigma=0.

  do n=1,400
     call impactstirring(nz,z,bigstep,sigma,1000)

     ! optional diagnostics
     m0 = colint(sigma,z,nz,1,nz)
     m1 = colint(sigma*z,z,nz,1,nz)/m0
     m2 = colint(sigma*z*z,z,nz,1,nz)/m0

     write(30,'(999(f8.3,1x))') sigma(:)
     write(31,*) n*bigstep,m0,m1,m2
  enddo
  !write(30,'(999(f8.3,1x))') sigma(:)
end program testgardening



