program mars_fast
!***********************************************************************
! Retreat and growth of subsurface ice on Mars
! orbital elements remain constant
!***********************************************************************
  use miscparameters, only : pi, d2r, NMAX, marsDay, solsperyear 
  use allinterfaces
  implicit none
  integer, parameter :: NP=1   ! # of sites
  integer nz, i, k, ierr
  real(8) zfac, zmax, delta, z(NMAX), icetime, porosity, icefrac
  real(8), dimension(NP) :: latitude, albedo, thIn, rhoc
  real(8), dimension(NP) :: pfrost, p0, htopo
  real(8) newti, stretch, newrhoc, ecc, omega, eps, timestep
  real(8), dimension(NP) :: zdepthF, zdepthE, zdepthT, zdepthG
  real(8), dimension(NMAX,NP) :: porefill
  real(8), dimension(NP) ::  Tb, Tmean1, Tmean3, avrho1, mCO2mean
  real(8) tmax, tlast
  character(10) ext
  real(8), external :: smartzfac

  if (iargc() /= 1) then
     stop 'USAGE: icages ext'
  endif
  call getarg( 1, ext )

  if (NP>100) stop 'subroutine icelayer_mars cannot handle this many sites'

  ! parameters that never ever change
  nz=60
  porosity = 0.4d0  ! porosity of till
  !rhoc(:) = 1500.*800.  ! will be overwritten
  zdepthT(:) = -9999.
  icefrac = 0.9
  tmax = 10000.
  tlast = 0.

  open(unit=21,file='lats.'//ext,action='read',status='old',iostat=ierr)
  if (ierr /= 0) then
     print *,'File lats.'//ext,'not found'
     stop
  endif
  do k=1,NP
     read(21,*) latitude(k),albedo(k),thIn(k),htopo(k)
     ! empirical relation from Mellon & Jakosky
     rhoc(k) = 800.*(150.+100.*sqrt(34.2+0.714*thIn(k))) 
  enddo
  close(21)

  ! set eternal grid
  zmax = 8.
  zfac = smartzfac(nz,zmax,6,0.032d0)
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z.'//ext,action='write',status='unknown')
  write(30,'(999(f8.5,1x))') z(1:nz)
  close(30)

  ecc = 0.0934;  eps = 25.19*d2r;  omega = 250.87*d2r   ! today
  pfrost(:) = 0.16d0
  ! total atmospheric pressure
  !p0(:) = 600.
  ! presently 520 Pa at zero elevation (Smith & Zuber, 1998)
  do k=1,NP
     p0(k)=520*exp(-htopo(k)/10800.)
  enddo
  timestep = 50.  ! must be integer fraction of 1 ka
  icetime = -tmax-timestep  ! earth years
  
  ! initializations 
  Tb(:) = -9999.
  zdepthF(:) = -9999.

  porefill(1:nz,1:NP) = 0.
  zdepthT(1:NP) = -9999.
  !zdepthT(1:NP) = 1.

  print *,'RUNNING MARS_FAST'
  print *,'Global model parameters:'
  print *,'nz=',nz,' zfac=',zfac,'zmax=',zmax
  print *,'porosity=',porosity
  print *,'starting at time',icetime,'years'
  print *,'time step=',timestep,'years'
  print *,'eps=',eps/d2r,'ecc=',ecc,'omega=',omega/d2r
  print *,'number of sites=',NP
  print *,'Site specific parameters:'
  do k=1,NP
     if (NP>1) print *,'  Site ',k
     print *,'  latitude (deg)',latitude(k),' rho*c (J/m^3/K)',rhoc(k),' thIn=',thIn(k)
     print *,'  total pressure=',p0(k),'partial pressure=',pfrost(k)
     delta = thIn(k)/rhoc(k)*sqrt(marsDay/pi)
     print *,'  skin depths (m)',delta,delta*sqrt(solsperyear)
     call soilthprop(porosity,1.d0,rhoc(k),thIn(k),1,newrhoc,newti)
     stretch = (newti/thIn(k))*(rhoc(k)/newrhoc)
     do i=1,nz
        if (z(i)<delta) cycle
        print *,'  ',i-1,' grid points within diurnal skin depth'
        exit
     enddo
     print *,'  ',zmax/(sqrt(solsperyear)*delta),'times seasonal dry skin depth'
     print *,'  ',zmax/(sqrt(solsperyear)*delta*stretch),'times seasonal filled skin depth'
     print *,'  Initial ice depth=',zdepthT(k)
     print *
  enddo
  call outputmoduleparameters
  print *

  ! open and name all output files
  open(unit=34,file='subout.'//ext,action='write',status='unknown')
  open(unit=36,file='depthF.'//ext,action='write',status='unknown')
  open(unit=37,file='depths.'//ext,action='write',status='unknown')

  print *,'Equilibrating initial temperature'
  do i=1,4
     call icelayer_mars(0d0,nz,NP,thIn,rhoc,z,porosity,pfrost,Tb,zdepthF, &
       &  zdepthE,porefill(1:nz,:),Tmean1,Tmean3,zdepthG, &
       &  latitude,albedo,p0,ecc,omega,eps,icefrac,zdepthT,avrho1,mCO2mean)
  enddo

  print *,'History begins here'
  porefill(1:nz,1:NP) = 0.
  zdepthT(1:NP) = -9999.
  !zdepthT(1:NP) = 1.

  do
     call icelayer_mars(timestep,nz,NP,thIn,rhoc,z,porosity,pfrost,Tb,zdepthF, &
          & zdepthE,porefill(1:nz,:),Tmean1,Tmean3,zdepthG, & 
          & latitude,albedo,p0,ecc,omega,eps,icefrac,zdepthT,avrho1,mCO2mean)
     icetime = icetime+timestep
     if (abs(mod(icetime/100.,1.d0))<1.e-3) then ! output every 1000 years
        do k=1,NP
           !write(36,*) icetime,latitude(k),zdepthF(k),porefill(1:nz,k)
           ! compact output format
           write(36,'(f10.0,2x,f7.3,1x,f11.5,1x)',advance='no') & 
                & icetime,latitude(k),zdepthF(k)
           call compactoutput(36,porefill(:,k),nz)
           write(37,501) icetime,latitude(k),zdepthT(k), &
                & Tmean1(k),Tmean3(k),zdepthG(k),avrho1(k)
        enddo
     endif
     print *,icetime
     if (icetime>=tlast) exit
  enddo

  close(34)
  close(36); close(37)

501 format (f10.0,2x,f7.3,2x,f10.4,2(2x,f6.2),2x,f9.3,2x,g11.4)
  
end program mars_fast


