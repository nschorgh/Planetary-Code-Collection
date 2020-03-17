!************************************************************************
! Subroutines for fast asteroid method
!************************************************************************


subroutine icelayer_asteroid(bigstep,NP,z,porosity,icefrac,Tinit, &
     & zdepthT,Tmean1,Tmean3,Tmin,Tmax,latitude,albedo,ecc,omega,eps,S0)
!************************************************************************
! bigstep = time step [Earth years]
! latitude  [degree]
! eps = axis tilt [radians]  
! S0 = solar constant relative to present
!************************************************************************
  use constants, only : d2r, NMAX
  use body, only : Tnominal, nz
  use allinterfaces
  implicit none
  integer, intent(IN) :: NP
  real(8), intent(IN) :: bigstep
  real(8), intent(IN) :: z(NMAX), porosity, icefrac
  logical, intent(IN) :: Tinit
  real(8), intent(INOUT) :: zdepthT(NP), Tmean1(NP), Tmean3(NP)
  real(8), intent(OUT) :: Tmin(NP), Tmax(NP)
  real(8), intent(IN) :: latitude(NP), albedo(NP), ecc, omega, eps, S0

  integer k, typeT, j, jump
  real(8) ti(NMAX), rhocv(NMAX), diam
  real(8) Deff, deltaz, Diff0, avrho
  real(8), dimension(nz) :: Diff
  real(8), SAVE :: zdepth_old(100)  ! NP<=100

  do k=1,NP   ! big loop over sites

     typeT = gettype(zdepthT(k),nz,z)

     ! assign/update property profiles
     call assignthermalproperties1(nz,z,Tnominal,porosity,ti,rhocv,icefrac,zdepthT(k))
     diam = 100e-6  ! assumed grain diameter in mantle, used to calculate mean-free path
     Diff0 = vapordiffusivity(diam,porosity,Tnominal) ! surface
     do j=1,nz
        !if (z(j)>0.5) diam=1e-3  ! coarser below 0.5m
        Diff(j) = vapordiffusivity(diam,porosity,Tnominal) 
        if (z(j)>zdepthT(k)) then
           Diff(j) = 0.
        endif
     end do
     
     ! run thermal model
     call ajsub_asteroid(latitude(k)*d2r, albedo(k), z, ti, rhocv, & 
          &     ecc, omega, eps, S0, typeT, avrho, &
          &     Tinit, Tmean1(k), Tmean3(k), Tmin(k), Tmax(k))

     ! run ice evolution model
     if (typeT<=1) then
        Deff = Diff0
     else
        deltaz = colint(spread(1d0,1,nz),z,nz,1,typeT-1)  ! for normalization
        Deff = deltaz/colint(1./Diff,z,nz,1,typeT-1) 
        if (minval(Diff(1:typeT-1))<=0.) then
           stop 'D_EFF PROBLEM'
        endif
     endif
     call icechanges(nz,z(:),avrho,Deff,bigstep,zdepthT(k),porosity,icefrac)

     ! diagnose
     if (zdepthT(k)>=0.) then
        jump = 0
        do j=1,nz
           if (zdepth_old(k)<z(j).and.zdepthT(k)>z(j)) jump=jump+1
        end do
     else
        jump=-9
     endif
     write(34,'(f12.2,1x,f6.2,1x,f11.5,1x,g11.4,1x,i3,1x,g10.4)') &
          &        bigstep,latitude(k),zdepthT(k),avrho,jump,Deff
     zdepth_old(k) = zdepthT(k)

  end do  ! end of big loop
end subroutine icelayer_asteroid



subroutine ajsub_asteroid(latitude, albedo, z, ti, rhocv, ecc, omega, &
     &     eps, S0, typeT, rhosatav, Tinit, Tmean1, Tmean3, Tmin, Tmaxi)
!***********************************************************************
!  A 1D thermal model that also returns various time-averaged quantities
!
!  Tinit = initalize if .true., otherwise use Tmean1 and Tmean3
!***********************************************************************
  use constants
  use body, only : EQUILTIME, dt, solsperyear, Fgeotherm, semia, nz, emiss, solarDay
  use allinterfaces, only : flux_noatm
  implicit none
  real(8), intent(IN) :: latitude  ! [radians]
  real(8), intent(IN) :: albedo, z(NMAX)
  real(8), intent(IN) :: ti(NMAX), rhocv(NMAX)
  real(8), intent(IN) :: ecc, omega, eps, S0
  integer, intent(IN) :: typeT
  real(8), intent(OUT) :: rhosatav   ! annual mean vapor density
  logical, intent(IN) :: Tinit
  real(8), intent(INOUT) :: Tmean1, Tmean3
  real(8), intent(OUT) :: Tmin, Tmaxi
  integer nsteps, n, j, nm
  real(8), parameter :: zero = 0.
  real(8) tmax, time, Qn, Qnp1, tdays
  real(8) orbitR, orbitLs, orbitDec, HA
  real(8) Tsurf, Fsurf, T(NMAX)
  real(8) Tmean0, S1, coslat, Evap
  real(8), external :: psv, sublrate
  
  ! initialize
  if (Tinit) then 
     S1=S0*1365./semia**2  ! must match solar constant defined in flux_noatm
     coslat = max(cos(latitude),cos(latitude+eps),cos(latitude-eps))
     Tmean0 = (S1*(1.-albedo)*coslat/(pi*emiss*sigSB))**0.25 ! estimate
     Tmean0 = Tmean0-5.
     if (Tmean0<50.) Tmean0=50.
     print *,Tmean0,S1,latitude,cos(latitude)
     write(*,*) '# initialized with temperature estimate of',Tmean0,'K'
     write(34,*) '# initialized with temperature estimate of',Tmean0,'K'
     T(1:nz) = Tmean0 
     Tsurf = Tmean0
     tmax = 3*EQUILTIME*solsperyear
  else
     do concurrent (j=1:nz)
        T(j) = (Tmean1*(z(nz)-z(j))+Tmean3*z(j))/z(nz)
     end do
     Tsurf = Tmean1
     tmax = EQUILTIME*solsperyear
  endif
  Fsurf=0.

  nsteps=int(tmax/dt)       ! calculate total number of timesteps

  nm=0
  Tmean1=0.; Tmean3=0.
  rhosatav = 0.
  Tmin=+1e32; Tmaxi=-9.
  Evap = 0.

  time=0.
  call generalorbit(0.d0,semia,ecc,omega,eps,orbitLs,orbitDec,orbitR)
  HA = 2.*pi*time            ! hour angle
  Qn = S0*(1-albedo)*flux_noatm(orbitR,orbitDec,latitude,HA,zero,zero)
  !----loop over time steps 
  do n=0,nsteps-1
     time = (n+1)*dt         !   time at n+1 
     tdays = time*(solarDay/86400.) ! parenthesis may improve roundoff
     call generalorbit(tdays,semia,ecc,omega,eps,orbitLs,orbitDec,orbitR)
     HA = 2.*pi*mod(time,1.d0)  ! hour angle
     Qnp1 = S0*(1-albedo)*flux_noatm(orbitR,orbitDec,latitude,HA,zero,zero)
     
     call conductionQ(nz,z,dt*solarDay,Qn,Qnp1,T,ti,rhocv,emiss, &
          &           Tsurf,Fgeotherm,Fsurf)
     Qn = Qnp1
     
     if (time>=tmax-solsperyear) then
        Tmean1 = Tmean1+Tsurf
        Tmean3 = Tmean3+T(nz)
        if (typeT>0 .and. typeT<=nz) then
           rhosatav = rhosatav+psv(T(typeT))/T(typeT)
        end if
        Evap = Evap + sublrate(Tsurf)
        nm = nm+1

        if (Tsurf<Tmin) Tmin=Tsurf
        if (Tsurf>Tmaxi) Tmaxi=Tsurf
     endif

  end do  ! end of time loop
  
  Tmean1 = Tmean1/nm; Tmean3 = Tmean3/nm
  rhosatav = rhosatav/nm
  rhosatav = rhosatav*18./8314.

  if (typeT<=0 .or. typeT>nz) rhosatav = -9999.
end subroutine ajsub_asteroid



subroutine outputmoduleparameters
  use body
  implicit none
  print *,'Global parameters stored in modules'
  print *,'  Ice bulk density',icedensity,'kg/m^3'
  print *,'  dt=',dt,'solar days'
  print *,'  Fgeotherm=',Fgeotherm,'W/m^2'
  print *,'  Emissivity of surface=',emiss
  print *,'  Thermal model equilibration time',EQUILTIME,'orbits'
  print *,'  Semimajor axis',semia
  print *,'  Solar day',solarDay,'Sols per year',solsperyear
  print *,'  Vertical grid: nz=',nz,' zfac=',zfac,'zmax=',zmax
end subroutine outputmoduleparameters



subroutine icechanges(nz,z,avrho,Deff,bigstep,zdepthT,porosity,icefrac)
!***********************************************************************
! advances ice interface, allows for deflation
!***********************************************************************
  use body, only : icedensity
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), avrho, Deff, bigstep, porosity, icefrac
  real(8), intent(INOUT) :: zdepthT
  real(8) zdepthTnew, buf, bigdtsec, beta

  if (zdepthT<0.) return   ! no ice anywhere

  bigdtsec = bigstep*86400*365.24

  ! advance ice table
  if (icefrac>porosity) then  ! deflation
     beta = (1-icefrac)/(1-porosity)
  else  ! no deflation
     beta = 1.
  endif
  buf = Deff*avrho*beta/(icefrac*icedensity)
  zdepthTnew = sqrt(2*buf*bigdtsec + zdepthT**2)
  print *,'# advance of ice table',zdepthT,zdepthTnew

  zdepthT = zdepthTnew
  if (zdepthT>z(nz)) zdepthT=-9999.
  
end subroutine icechanges



subroutine assignthermalproperties1(nz,z,Tnom,porosity,ti,rhocv,icefrac,zdepthT)
!***********************************************************************
! assign thermal properties of soil
!***********************************************************************
  use body, only : icedensity
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), Tnom, porosity
  real(8), intent(OUT) :: ti(nz), rhocv(nz)
  real(8), intent(IN), optional :: icefrac, zdepthT
  real(8), parameter :: rhodry = 2500  ! bulk density
  real(8), parameter :: kice=4.6, cice=1145   ! 140K
  integer j
  real(8) cdry  ! heat capacity of dry regolith
  real(8) k(nz)  ! thermal conductivity
  real(8) thIn
  real(8), external :: heatcapacity

  cdry = heatcapacity(Tnom)
  thIn = 15.
  do j=1,nz
     rhocv(j) = (1.-porosity)*rhodry*cdry
     k(j) = thIn**2/rhocv(j) 
  end do
  if (present(icefrac) .and. present(zdepthT)) then
     do j=1,nz
        if (z(j)>zdepthT) then
           !k(j) = 1./((1.-icefrac)/k(j) + icefrac/kice) 
           k(j) = k(j) + icefrac*kice  ! icefrac <= porosity
           rhocv(j) = rhocv(j) + icedensity*cice*icefrac
        endif
     end do
  end if

  ti(1:nz) = sqrt(k(1:nz)*rhocv(1:nz))
end subroutine assignthermalproperties1

