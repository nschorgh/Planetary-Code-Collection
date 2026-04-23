PROGRAM trojans_fast3
!*************************************************************************
! Asynchronously coupled model of temperature and near-surface ice retreat
!    on airless bodies
!
!  - diurnally resolved 1D thermal model
!  - ice loss by vapor diffusion through porous layer
!  - mixture of silicates, ice, and void spaces
!  - increasing solar luminosity
!  x no deflation (ice fraction larger than porosity)  
!  x no redistribution of ice within ice-rich layer
!  x no impact stirring
!  x only a single location, NP=1
!
! builds on asteroid_fast1
! written by Norbert Schorghofer 2025-2026
!*************************************************************************
  use body, only : pi, d2r, nz, zfac, zmax, Orbit
  use allinterfaces
  implicit none
  integer SPINUPN   ! # number of spin-up steps
  real(8) spinupfac 
  parameter(SPINUPN=15, spinupfac=2.)
  integer i, earliest, ierr
  real(8) tstart  ! [Earth years]
  real(8) z(nz), icetime
  real(8) timestep ! time step at the end of spin-up [Earth years]
  real(8) bigstep, bssum, porosity(nz)
  real(8) icefrac(nz) ! volumetric ice concentration relative to total volume
  real(8) latitude, zdepthT
  real(8), dimension(0:nz) :: Tmean, Tmini, Tmaxi
  character(4) ext

  ! latitudes
  if (command_argument_count() /= 1) stop 'USAGE: asteroid_fast ext'
  call get_command_argument( 1, ext )
  open(unit=21,file='lats.'//ext,action='read',status='old',iostat=ierr)
  if (ierr /= 0) then
     print *,'File lats.'//ext,'not found'
     stop
  endif
  read(21,*,iostat=ierr) latitude
  if (latitude>90. .or. latitude<-90.) then
     print *,latitude
     stop 'Latitude out of range'
  end if
  close(21)

  tstart = 4.4e9  ! Earth years
  timestep = 1e8  ! Earth years
  zdepthT = 0. ! initial ice depth

  orbit%omega = 0.*d2r

  ! set eternal grid
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z.'//ext,action='write',status='unknown')
  write(30,'(*(f8.5,1x))') z(1:nz)
  close(30)

  porosity(:) = 0.6d0   ! ice-free porosity
  icefrac(:)  = 0.3     ! initial ice fraction relative to total volume

  if (maxval(icefrac + porosity)>1.) stop 'icefrac>porosity>1'
  
  print *,'RUNNING FAST ASTEROID MODEL'
  print *,'Starting at time',tstart,'years'
  print *,'Time step=',timestep,'years'
  print *,'Spinup',SPINUPN,spinupfac
  print *,'eps=',Orbit%eps/d2r,'ecc=',Orbit%ecc,'omega=',orbit%omega/d2r
  print *,'Site specific parameters:'
  print *,'  Latitude (deg)',latitude
  call outputskindepths(z,porosity,icefrac)
  print *,'  Initial ice depth=',zdepthT
  print *,'  Porosity=',minval(porosity),maxval(porosity)
  print *,'  Ice fraction=',minval(icefrac),maxval(icefrac)
  print *
  call outputmoduleparameters
  print *

  open(unit=34,file='subout.'//ext,action='write',status='unknown')
  open(unit=37,file='depths.'//ext,action='write',status='unknown')

  earliest = nint(tstart/timestep)
  icetime = -earliest*timestep

  print *,'Equilibrating initial temperature'
  !icetime = -tstart
  call icelayer_asteroid(0d0,z,porosity,icefrac,.true., &
       & zdepthT,Tmean,Tmini,Tmaxi,latitude,Orbit,faintsun(icetime))
  write(37,501) icetime,latitude,zdepthT,Tmean(0),Tmean(nz),Tmini(0),Tmaxi(0)

  icetime = - earliest*timestep
  print *,icetime
  print *,'Spin-up begins here'
  ! sum_{j=1,n} a^j = a (a^n-1)/(a-1)
  bssum = spinupfac*(spinupfac**SPINUPN-1)/(spinupfac-1.)
  print *,'Spin-up', SPINUPN,'steps over',timestep,'years'
  do i=1,SPINUPN
     bigstep = spinupfac**i/bssum*timestep
     icetime = icetime + bigstep
     call icelayer_asteroid(bigstep,z,porosity,icefrac,.false., &
          & zdepthT,Tmean,Tmini,Tmaxi,latitude,Orbit,faintsun(icetime))
     print *,i,'of',SPINUPN,'  ',bigstep,zdepthT,orbit%omega/d2r
     ! variables were evaluated at previous time step
     write(37,501) icetime,latitude,zdepthT, &
          & Tmean(0),Tmean(nz),Tmini(0),Tmaxi(0)
     orbit%omega = mod(orbit%omega + 36.*d2r, 2*pi)  ! sweep
  end do

  icetime = -(earliest-1)*timestep
  print *,icetime
  do 
     icetime = icetime + timestep
     call icelayer_asteroid(timestep,z,porosity,icefrac,.false., &
          & zdepthT,Tmean,Tmini,Tmaxi,latitude,Orbit,faintsun(icetime))
     ! variables were evaluated at previous time step
     write(37,501) icetime,latitude,zdepthT,Tmean(0),Tmean(nz),Tmini(0),Tmaxi(0)
     print *,icetime
     if (icetime>=0.) exit
     orbit%omega = mod(orbit%omega + 36.*d2r, 2*pi)  ! sweep
  end do
  close(34)
  close(37)
501 format (f12.0,2x,f7.3,4x,f12.6,4(2x,f6.2))
  
  block
    integer j
    open(unit=38,file='profiles.'//ext,action='write',status='unknown')
    write(38,502) latitude,0.,Tmini(0),Tmean(0),Tmaxi(0)
    do j=1,nz
       write(38,502) latitude,z(j),Tmini(j),Tmean(j),Tmaxi(j)
    enddo
    close(38)
502 format (f7.3,1x,f12.6,3(1x,f7.3)) 
  end block
  
END PROGRAM trojans_fast3



subroutine outputskindepths(z,porosity,icefrac)
  ! diagnostics only
  use body, only : pi, nz, zmax, Tnominal, Orbit
  use allinterfaces
  implicit none
  real(8), intent(IN) :: z(nz), porosity(nz), icefrac(nz)
  integer i
  real(8) delta, stretch, newrhoc, newti, rhoc
  real(8) rhocv(nz), ti(nz), thIn, T(nz), solsperorbit

  T = spread(Tnominal,1,nz)
  call assignthermalproperties3(nz,z,T,porosity,ti,rhocv)
  thIn = ti(1); rhoc=rhocv(1)
  print *,'Thermal inertia=',thIn
  call assignthermalproperties3(nz,z,T,porosity,ti,rhocv,icefrac,0.d0)
  newti = ti(1); newrhoc=rhocv(1)

  solsperorbit = sols_per_orbit( Orbit%semia, Orbit%solarDay )
  delta = thIn/rhoc*sqrt( Orbit%solarDay / pi )
  stretch = (newti/thIn)*(rhoc/newrhoc)
  print *,'  ice-free skin depths - diurnal & seasonal', &
       & delta,delta*sqrt(solsperorbit)
  print *,'  ice-filled skin depths - diurnal & seasonal', &
       & delta*stretch,delta*sqrt(solsperorbit)*stretch
  do i=1,nz
     if (z(i)<delta) cycle
     print *,'  ',i-1,' grid points within ice-free diurnal skin depth'
     exit
  end do
  print *,'  zmax=',zmax/(sqrt(solsperorbit)*delta), &
       & 'times seasonal ice-free skin depth'
  print *,'  zmax=',zmax/(delta*stretch), &
       & 'times diurnal ice-filled skin depth'
  print *,'  zmax=',zmax/(sqrt(solsperorbit)*delta*stretch), &
       & 'times seasonal ice-filled skin depth'
  write(*,'(3x,a,3(1x,f6.1))') 'Nominal thermal inertia extremes',thIn,newti
  if (i<=5) stop 'Not enough grid points'
end subroutine outputskindepths

