program exper_fast
!***********************************************************************
! Retreat and growth of subsurface ice 
! prescribed surface temperatures
!***********************************************************************
  use miscparameters
  use allinterfaces
  implicit none
  integer nz, i
  real(8) zfac, zmax, z(NMAX), icetime, porosity
  real(8) pfrost   ! partial pressure of H2O
  real(8) thIn, rhoc, timestep, porefill(NMAX)
  real(8) zdepthT, zdepthF, zdepthE, zdepthG
  real(8) Tmean, Tampl, Diff, tmax, icefrac
  namelist /stuff/ nz,zmax,zfac,thIn,rhoc,pfrost,Tmean,Tampl,Diff, &
       & porosity,tmax,timestep
  !real(8), external :: colint  ! defined in allinterfaces.mod

  ! parameters that do not change
  open(unit=20,file='input_fast.par')
  read(20, nml=stuff)
  write(6,nml=stuff)
  call outputmoduleparameters

  ! set eternal grid
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z',action='write',status='unknown')
  write(30,'(999(f8.5,1x))') z(1:nz)
  close(30)

  icetime=0.
  porefill(1:nz) = 0.
  !zdepthT=-9999.
  zdepthT = 2.
  icefrac= 0.9
  zdepthF = -9999.
  if (zdepthT>=0.) then
     where(z(1:nz)>zdepthT)
        porefill(1:nz) = icefrac/porosity  ! for output
     end where
  end if
  print *,icetime,zdepthT,zdepthF
  write(36,'(f10.0,2x)',advance='no') icetime
  call compactoutput(36,porefill,nz)

  do i=1,10000000
     icetime = i*timestep
     if (icetime>tmax) exit
     !if (icetime>10000.) pfrost=0.
     call icelayer_exper(timestep,nz,thIn,rhoc,z,porosity,pfrost, &
          & zdepthT, icefrac, zdepthF, zdepthE, porefill(1:nz),  &
          & Tmean,Tampl,Diff,zdepthG)
     if (abs(mod(icetime/2000.,1.d0))<1.e-3) then ! output every 2000 years
        !write(36,*) icetime,porefill(1:nz)
        write(36,'(f10.0,2x)',advance='no') icetime
        call compactoutput(36,porefill,nz)
     endif
     write(35,*) icetime,zdepthT,zdepthF,zdepthE, &
          & colint(porefill(1:nz),z,nz,1,nz),zdepthG
     print *,icetime,zdepthT,zdepthF
  enddo
end program exper_fast










