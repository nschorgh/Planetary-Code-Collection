module globalparams
  implicit none
  real(8), parameter :: sigSB = 5.6704d-8
  real(8), parameter :: rockdensity = 2500.
  real(8), parameter :: Fgeotherm = 0.018
  integer, parameter :: STEPSPERSOL = 96
end module globalparams


program ensemble
!***********************************************************************
! run ensemble of thermal models for the Moon
!***********************************************************************
  use globalparams
  implicit none
  integer :: k, idum=-123
  real(8) lat, H, ell, albedo, emiss, x
  real(8) porosity_s, porosity_d, k_s, k_d
  real(8), external :: ran2
  character(6) extc

  if (command_argument_count() > 0) then
     ! take latitude from command line argument instead
     ! a.out +5
     call get_command_argument( 1, extc )
     read(extc,'(f6.2)') lat
  else
     lat = 0.5
     extc = ''
  endif
  print *,'Latitude=',lat

  open(unit=40,file='params'//trim(extc)//'.model',action='write')
  open(unit=22,file='Tsurface'//trim(extc)//'.model',action='write')

  ! default values which may or may not be overwritten
  H = 0.07
  ell = 200e-6
  albedo = 0.12
  emiss = 0.95
  porosity_s = 0.56; porosity_d = 0.28
  k_s = 7.4e-4; k_d = 3.4e-3
  k_d = 6e-3
  
  do k=1,10000  ! size of ensemble
     !x = ran2(idum)
     !H = 0.01*x + 0.15*(1-x)
     x = ran2(idum)
     ell = 0.01e-3*x + 1e-3*(1-x)
     x = ran2(idum)
     albedo = 0.05*x + 0.25*(1-x)
     !x = ran2(idum)
     !porosity_s = 0.90*x + 0.30*(1-x)
     !x = ran2(idum)
     !porosity_d = 0.50*x + 0.20*(1-x)
     x = ran2(idum)
     k_s = 2e-4*x + 10e-4*(1-x)
     x = ran2(idum)
     k_d = 1e-3*x + 10e-3*(1-x)
     !x = ran2(idum)
     !emiss = 0.9*x + 1.0*(1-x)
     
     call one(lat,H,ell,albedo,porosity_s,porosity_d,k_s,k_d,emiss)
     write(40,'(f6.4,1x,g9.3,1x,f6.4,2(1x,f6.4),2(1x,g9.3),1x,f5.3)') &
          & H,ell,albedo,porosity_s,porosity_d,k_s,k_d,emiss
  end do

  close(22)
  close(40)
end program ensemble



subroutine one(lat,H,ell,albedo,porosity_s,porosity_d,k_s,k_d,emiss)
!***********************************************************************
! run one lunar thermal model
!***********************************************************************
  use globalparams
  implicit none
  real(8), intent(IN) :: lat, H, ell, albedo, porosity_s, porosity_d, k_s, k_d, emiss
  real(8) Tsurf_out(STEPSPERSOL)
  real(8), parameter :: pi = 3.1415926535897932, d2r = pi/180.
  real(8), parameter :: lunarDay = 86400.*29.53
  integer nz
  real(8) zmax, zfac
  parameter(nz=30, zmax=0.5, zfac=1.07)
  !parameter(nz=40, zmax=0.5, zfac=1.05)
  integer nsteps, n, i, j
  real(8) T(nz), tmax, time, Tsurf, Qn, ltime
  real(8) Qnp1, dtsec, sunR, decl, HA
  real(8), dimension(nz) :: ti, rhocv, z, k
  real(8) Fsurf, porosity, rho, kc
  real(8), external :: flux_moon, heatcapacity, radconductivity, twolayers
  logical :: firsttime = .TRUE.
  
  dtsec = lunarDay/STEPSPERSOL
  tmax = 50.  ! lunations
  
  decl = 0.
  !decl = -1.5
  sunR = 1.
  
  print *,'Latitude=',lat,'ell=',ell,'albedo=',albedo
  if (lat<-90. .or. lat>+90.) stop 'latitude out of range'
  if (ell<1e-7 .or. ell>1) stop 'unreasonable particle size ell'
  
  nsteps = nint(tmax*lunarDay/dtsec) + STEPSPERSOL/2  ! calculate total number of timesteps

  call setgrid(nz,z,zmax,zfac)
  if (firsttime) then
     open(unit=20,file='z',action='write')
     write(20,'(999(f8.6,1x))') z(1:nz)
     close(20)

     open(unit=21,file='ltimes',action='write')
  end if
  
  block  ! initialize temperatures
    real(8) geof, Tinit
    real(8) Tnom, delta, knz, thIn, rhoc, ztell

    geof = cos(lat*d2r)/pi
    Tinit = (1370*(1.-albedo)*geof/sigSB)**0.25 - 25. 
    if (geof<=0.) Tinit=50.  ! the poles
    T(1:nz) = Tinit
    Tsurf = Tinit

    ! extra checks
    Tnom = Tinit
    ztell = z(10)
    kc = twolayers(k_s, k_d, ztell, H)
    porosity = twolayers(porosity_s, porosity_d, ztell, H)
    knz = kc + radconductivity(Tnom,ell,emiss,porosity)
    rhoc = rockdensity * (1-porosity) * heatcapacity(Tnom)
    thIn = sqrt( knz*rhoc )
    delta = (thIn/rhoc)*sqrt(lunarDay/pi)  ! skin depth
    
    do i=1,nz
       if (z(i)<delta) cycle
       !print *,i-1,' grid points within diurnal skin depth',delta
       exit
    enddo
  end block
  
  ! time=0, n=0
  HA = 0.
  Qn = flux_moon(sunR,decl*d2r,lat*d2r,HA,albedo)
  j = 0

  do n=0,nsteps-1  ! loop over time steps
     time = (n+1)*dtsec   ! time at n+1
     
     HA = time/lunarDay*2*pi  ! hour angle
     Qnp1 = flux_moon(sunR,decl*d2r,lat*d2r,HA,albedo)

     do i=1,nz ! update thermal properties
        !kc = solidconductivity(z(i),H)
        kc = twolayers(k_s, k_d, z(i), H)
        porosity = twolayers(porosity_s, porosity_d, z(i), H)
        k(i) = kc + radconductivity(T(i), ell, emiss, porosity)
        rho = rockdensity*(1-porosity)  ! 1100 ... 1800
        rhocv(i) = rho * heatcapacity(T(i))
        ti(i) = sqrt( k(i)*rhocv(i) )
     enddo
     
     call conductionQ(nz,z,dtsec,Qn,Qnp1,T(:),ti(:),rhocv(:),emiss,Tsurf,Fgeotherm,Fsurf)
     Qn = Qnp1
     
     if (n >= nsteps-STEPSPERSOL) then  ! last lunation
        ltime = mod(time/lunarDay*360-1e-10 + 180, 360.d0)
        !write(22,'(f6.3,1x,f8.3,1(1x,f7.2))') ltime*24./360.,Qn,Tsurf
        j = j+1
        Tsurf_out(j) = Tsurf
        if (firsttime) write(21,'(f6.3)') ltime*24./360.
     end if
  end do

  if (firsttime) close(21)
  write(22,'(*(1x,f7.2))') Tsurf_out

  firsttime = .FALSE.
end subroutine one

