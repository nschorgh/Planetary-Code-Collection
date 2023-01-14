PROGRAM asteroid_fast1
!*************************************************************************
! Asynchronously coupled model of temperature and near-surface ice retreat
!    on airless bodies
!
!  - diurnally resolved 1D thermal model
!  - ice loss by vapor diffusion through porous layer
!  - mixture of silicates and ice (no void spaces)
!  - increasing solar luminosity
!  - allows for deflation (ice fraction larger than porosity)  
!  x no redistribution of ice within ice-rich layer
!  x no impact stirring
!
! mostly a simplified version of asteroid_fast2    2016-2017
!*************************************************************************
  use constants, only : pi, d2r
  use body, only : nz, zfac, zmax, ecc, eps
  use allinterfaces
  implicit none
  integer, parameter :: NP=1    ! # of sites
  integer SPINUPN   ! # number of spin-up steps
  real(8) spinupfac 
  parameter(SPINUPN=20, spinupfac=2.)
  integer i, k, earliest, ierr
  real(8) tstart  ! (earth) years
  real(8) z(nz), icetime, timestep
  real(8) bigstep, bssum, omega, porosity, icefrac
  real(8), dimension(NP) :: latitude, albedo, zdepthT
  real(8), dimension(NP) :: Tmean1, Tmean3, Tmin, Tmax
  character(4) ext

  ! latitudes
  if (command_argument_count() /= 1) stop 'USAGE: asteroid_fast ext'
  call get_command_argument( 1, ext )
  open(unit=21,file='lats.'//ext,action='read',status='old',iostat=ierr)
  if (ierr /= 0) then
     print *,'File lats.'//ext,'not found'
     stop
  endif
  do k=1,NP
     read(21,*,iostat=ierr) latitude(k),albedo(k)
  end do
  close(21)

  tstart = 4.5e9  ! Earth years
  timestep = 2e5  ! Earth years
  zdepthT(:) = 0.  ! initial ice depth

  omega = 0.*d2r

  ! set eternal grid
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z.'//ext,action='write',status='unknown')
  write(30,'(*(f8.5,1x))') z(1:nz)
  close(30)

  ! porosity can decrease with depth, but should be constant within stirring depth
  porosity = 0.4d0   ! ice-free porosity
  icefrac  = porosity

  print *,'RUNNING FAST ASTEROID MODEL'
  print *,'Starting at time',tstart,'years'
  print *,'Time step=',timestep,'years'
  print *,'Spinup',SPINUPN,spinupfac
  print *,'eps=',eps/d2r,'ecc=',ecc,'omega=',omega/d2r
  print *,'Number of sites=',NP
  print *,'Site specific parameters:'
  do k=1,NP
     if (NP>1) print *,'  Site ',k
     print *,'  Latitude (deg)',latitude(k)
     call outputskindepths(nz,z,zmax,porosity,icefrac)
     print *,'  Initial ice depth=',zdepthT(k)
     print *,'  Porosity=',porosity,' Ice fraction=',icefrac
     print *
  end do
  call outputmoduleparameters
  print *

  open(unit=34,file='subout.'//ext,action='write',status='unknown')
  open(unit=37,file='depths.'//ext,action='write',status='unknown')

  earliest = nint(tstart/timestep)
  icetime = -earliest*timestep

  print *,'Equilibrating initial temperature'
  !icetime = -tstart
  call icelayer_asteroid(0d0,NP,z,porosity,icefrac,.true., &
       &        zdepthT,Tmean1,Tmean3,Tmin,Tmax,latitude,albedo, &
       &        ecc,omega,eps,faintsun(icetime))
  do k=1,NP
     write(37,501) icetime,latitude(k),zdepthT(k), &
          &  Tmean1(k),Tmean3(k),Tmin(k),Tmax(k)
  end do

  icetime = - earliest*timestep
  print *,icetime
  print *,'Spin-up begins here'
  ! sum_{j=1,n} a^j = a (a^n-1)/(a-1)
  bssum = spinupfac*(spinupfac**SPINUPN-1)/(spinupfac-1.)
  print *,'Spin-up', SPINUPN,'steps over',timestep,'years'
  do i=1,SPINUPN
     bigstep = spinupfac**i/bssum*timestep
     icetime = icetime + bigstep
     call icelayer_asteroid(bigstep,NP,z,porosity,icefrac,.false., &
          & zdepthT,Tmean1,Tmean3,Tmin,Tmax,latitude,albedo, &
          & ecc,omega,eps,faintsun(icetime))
     print *,i,'of',SPINUPN,'  ',bigstep,zdepthT,omega/d2r
     do k=1,NP
        ! variables were evaluated at previous time step
        write(37,501) icetime,latitude(k),zdepthT(k), &
             & Tmean1(k),Tmean3(k),Tmin(k),Tmax(k)
     end do
     omega = mod(omega + 36.*d2r,2*pi)  ! sweep
  end do

  
  icetime = -(earliest-1)*timestep
  print *,icetime
  do 
     icetime = icetime + timestep
     call icelayer_asteroid(timestep,NP,z,porosity,icefrac,.false., &
          & zdepthT,Tmean1,Tmean3,Tmin,Tmax,latitude,albedo, &
          & ecc,omega,eps,faintsun(icetime))
     do k=1,NP
        ! variables were evaluated at previous time step
        write(37,501) icetime,latitude(k),zdepthT(k), &
             & Tmean1(k),Tmean3(k),Tmin(k),Tmax(k)
     end do
     print *,icetime
     !if (any(-icetime == (/ 4.498d9, 4.450d9, 4d9 /))) then  ! with 1e5
     if (any(-icetime == (/ 4.498d9, 4.460d9, 4d9 /))) then   ! with 2e5
     !if (any(-icetime == 2.0d9 - (/ 0.002d9, 0.05d9, 0.5d9 /))) then
        timestep = 10.*timestep
     endif
     if (icetime>=0.) exit
     omega = mod(omega + 36.*d2r,2*pi)  ! sweep
  end do

  close(34)
  close(37)

501 format (f12.0,2x,f7.3,4x,f12.6,4(2x,f6.2)) 

END PROGRAM asteroid_fast1



subroutine outputskindepths(nz,z,zmax,porosity,icefrac)
  ! diagnostics only
  use constants, only : pi
  use body, only : solarDay, semia, Tnominal
  use allinterfaces
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), zmax, porosity, icefrac
  integer i
  real(8) delta, stretch, newrhoc, newti, rhoc
  real(8) rhocv(nz), ti(nz), thIn, solsperyear

  call assignthermalproperties1(nz,z,Tnominal,porosity,ti,rhocv)
  thIn = ti(1); rhoc=rhocv(1)
  print *,'Thermal inertia=',thIn
  call assignthermalproperties1(nz,z,Tnominal,porosity,ti,rhocv,icefrac,0.d0)
  newti = ti(1); newrhoc=rhocv(1)

  solsperyear = sols_per_year(semia,solarDay)
  delta = thIn/rhoc*sqrt(solarDay/pi)
  stretch = (newti/thIn)*(rhoc/newrhoc)
  print *,'  ice-free skin depth - diurnal/seasonal',delta,delta*sqrt(solsperyear)
  do i=1,nz
     if (z(i)<delta) cycle
     print *,'  ',i-1,' grid points within ice-free diurnal skin depth'
     exit
  end do
  print *,'  zmax=',zmax/(sqrt(solsperyear)*delta), &
       & 'times seasonal ice-free skin depth'
  print *,'  zmax=',zmax/(sqrt(solsperyear)*delta*stretch), &
       & 'times seasonal filled skin depth'
  write(*,'(3x,a,3(1x,f6.1))') 'Nominal thermal inertia extremes',thIn,newti
  if (i<=5) stop 'Not enough grid points'
end subroutine outputskindepths

