!************************************************************************
! Subroutines for retreat of subsurface ice on airless bodies,
!      based on time-averaged vapor transport equations
!
! written by Norbert Schorghofer 2025-2026
!************************************************************************


subroutine icelayer_asteroid(bigstep,z,porosity,icefrac,Tinit, &
     & zdepthT,Tmean,Tmini,Tmaxi,latitude,Orbit,S0)
!************************************************************************
! bigstep = time step [Earth years]
! z = depths [m]
! zdetphT = depth of ice table or shallowest depth with perennial ice [m]  
! latitude  [degree]
! eps = axis tilt [radians]  
! S0 = solar constant relative to present
!************************************************************************
  use body, only : d2r, Tnominal, nz, diam, orbitp
  use allinterfaces
  implicit none
  real(8), intent(IN) :: bigstep
  real(8), intent(IN) :: z(nz), porosity(nz), icefrac(nz)
  logical, intent(IN) :: Tinit
  real(8), intent(INOUT) :: zdepthT, Tmean(0:nz)
  real(8), intent(OUT) :: Tmini(0:nz), Tmaxi(0:nz)
  real(8), intent(IN) :: latitude, S0
  type(orbitp), intent(IN) :: Orbit
  
  integer typeT, j, jump
  real(8) ti(nz), rhocv(nz), T(nz)
  real(8) ell0, elleff, deltaz, avSice
  real(8), dimension(nz) :: ell
  real(8), SAVE :: zdepth_old

  typeT = getfirst(zdepthT,nz,z)

  ! assign/update property profiles
  if (Tinit) then
     T(:) = spread(Tnominal,1,nz)
  else
     !T(:) = ( Tmean1*(z(nz)-z(:)) + Tmean3*z(:) ) / z(nz)
     T(:) = Tmean(1:nz)
  end if
  call assignthermalproperties3(nz,z(:),T,porosity, &
       &                        ti(:),rhocv(:),icefrac(:),zdepthT)
  ell0 = meanfreepathinsoil(diam,porosity(1)) ! surface
  do j=1,nz
     ell(j) = meanfreepathinsoil(diam,porosity(j)) 
     if (z(j)>zdepthT .and. zdepthT>=0.) then
        ell(j) = meanfreepathinsoil(diam,porosity(j)+icefrac(j)) 
     endif
  end do
     
  ! run thermal model
  call ajsub_asteroid(latitude*d2r, z, ti, rhocv, Orbit, S0, &
       &     typeT, avSice, Tinit, Tmean(:), Tmini(:), Tmaxi(:))

  ! run ice evolution model
  if (typeT<=1) then
     elleff = ell0
  else
     deltaz = colint(spread(1d0,1,nz),z,nz,1,typeT-1)  ! for normalization
     elleff = deltaz/colint(1./ell(:),z(:),nz,1,typeT-1) 
     if (minval(ell(1:typeT-1))<=0.) then
        stop 'ELL_EFF PROBLEM'
     endif
  endif
  call icechanges3(nz,z(:),avSice,elleff,bigstep,zdepthT,porosity,icefrac)

  ! diagnose
  if (zdepthT>=0.) then
     jump = 0
     do j=1,nz
        if (zdepth_old<z(j) .and. zdepthT>z(j)) jump=jump+1
     end do
  else
     jump=-9
  endif
  write(34,'(f12.2,1x,f6.2,1x,f11.5,1x,g11.4,1x,i3,1x,g10.4)') &
       &        bigstep,latitude,zdepthT,avSice,jump,elleff
  zdepth_old = zdepthT

end subroutine icelayer_asteroid



subroutine ajsub_asteroid(latitude, z, ti, rhocv, Orbit, S0, &
     &     typeT, avSice, Tinit, Tmean, Tmini, Tmaxi)
!************************************************************************
!  A 1D thermal model that also returns various time-averaged quantities
!
!  Tinit = initalize if .true., otherwise use Tmean
!************************************************************************
  use body, only : pi, EQUILTIME, dt, Fgeotherm, nz, emiss, albedo, orbitp
  use allinterfaces
  implicit none
  real(8), intent(IN) :: latitude  ! [radians]
  real(8), intent(IN) :: z(nz)
  real(8), intent(IN) :: ti(nz), rhocv(nz)
  type(orbitp), intent(IN) :: Orbit
  real(8), intent(IN) :: S0
  integer, intent(IN) :: typeT
  real(8), intent(OUT) :: avSice  ! annual mean sublimation rate [kg/m^2/s]
  logical, intent(IN) :: Tinit
  real(8), intent(INOUT) :: Tmean(0:nz)
  real(8), intent(OUT) :: Tmini(0:nz), Tmaxi(0:nz)
  real(8), parameter :: sigSB=5.6704e-8
  real(8), parameter :: mmass = 18.  ! H2O
  !real(8), parameter :: mmass = 27.0  ! HCN
  real(8), parameter :: zero = 0.
  integer nsteps, n, nm
  real(8) tmax, time, Qn, Qnp1, tdays
  real(8) orbitR, orbitLs, orbitDec, HA
  real(8) Tsurf, Fsurf, T(nz)
  real(8) Tmean0, S1, coslat, solsperorbit
  
  ! initialize
  solsperorbit = sols_per_orbit( orbit%semia, orbit%solarDay)
  if (Tinit) then
     S1 = S0*1365./ (orbit%semia)**2  ! must match solar constant in flux_noatm
     coslat = max( &
          cos(latitude), cos(latitude+orbit%eps), cos(latitude-orbit%eps), 0.d0)
     Tmean0 = (S1*(1.-albedo)*coslat/(pi*emiss*sigSB))**0.25 ! estimate
     Tmean0 = Tmean0-5.  ! lower due to finite-amplitude effect
     if (Tmean0<20. .or. Tmean0/=Tmean0) Tmean0=20.
     print *,Tmean0,S1,latitude,cos(latitude)
     write(*,*) '# initialized with temperature estimate of',Tmean0,'K'
     write(34,*) '# initialized with temperature estimate of',Tmean0,'K'
     T(1:nz) = Tmean0 
     Tsurf = Tmean0
     tmax = 3*EQUILTIME*nint(solsperorbit)
  else
     !T(:) = ( Tmean1*(z(nz)-z(:)) + Tmean3*z(:) ) / z(nz)
     !Tsurf = Tmean1
     T(:) = Tmean(1:nz)
     Tsurf = Tmean(0)
     tmax = EQUILTIME*nint(solsperorbit)
  endif
  Fsurf=0.

  nsteps=int(tmax/dt)       ! calculate total number of timesteps

  nm=0
  Tmean(:) = 0.
  avSice = 0.
  Tmini(:)=+1e32; Tmaxi(:)=-9.

  time=0.
  call generalorbit( 0.d0, orbit%semia, orbit%ecc, orbit%omega, orbit%eps, &
       & orbitLs, orbitDec, orbitR)
  HA = 2.*pi*time            ! hour angle
  Qn = S0*(1-albedo)*flux_noatm(orbitR,orbitDec,latitude,HA,zero,zero)
  !----loop over time steps 
  do n=0,nsteps-1
     time = (n+1)*dt         !   time at n+1 
     tdays = time*( orbit%solarDay / 86400. ) ! parenthesis may improve roundoff
     call generalorbit(tdays,orbit%semia,orbit%ecc,orbit%omega,orbit%eps, &
          & orbitLs,orbitDec,orbitR)
     HA = 2.*pi*mod(time,1.d0)  ! hour angle
     Qnp1 = S0*(1-albedo)*flux_noatm(orbitR,orbitDec,latitude,HA,zero,zero)
     
     call conductionQ(nz,z,dt*orbit%solarDay,Qn,Qnp1,T,ti,rhocv,emiss, &
          &           Tsurf,Fgeotherm,Fsurf)
     Qn = Qnp1

     if ( time >= tmax - nint(solsperorbit) ) then
        Tmean(0) = Tmean(0) + Tsurf
        Tmean(1:nz) = Tmean(1:nz) + T(1:nz)
        if (typeT>0 .and. typeT<=nz) then
           !rhosatav = rhosatav+psv(T(typeT))/T(typeT)
           avSice = avSice + psv(T(typeT)) / sqrt(T(typeT))
        end if
        nm = nm+1

        if (Tsurf<Tmini(0)) Tmini(0)=Tsurf
        if (Tsurf>Tmaxi(0)) Tmaxi(0)=Tsurf
        where (T<Tmini(1:nz)) Tmini(1:nz)=T(:)
        where (T>Tmaxi(1:nz)) Tmaxi(1:nz)=T(:)
     endif

  end do  ! end of time loop
  
  Tmean(:) = Tmean(:)/nm
  !rhosatav = rhosatav/nm
  !rhosatav = rhosatav*mmass/8314.46
  avSice = avSice/nm
  avSice = avSice*sqrt(mmass/(2*pi*8314.46))
  
  if (typeT<=0 .or. typeT>nz) avSice = -9999.
end subroutine ajsub_asteroid



subroutine icechanges3(nz,z,avSice,elleff,bigstep,zdepthT,porosity,icefrac)
!***********************************************************************
! advances ice interface, no deflation
!***********************************************************************
  use allinterfaces, only : getfirst
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), avSice, elleff, bigstep
  real(8), intent(IN) :: porosity(nz), icefrac(nz)
  real(8), intent(INOUT) :: zdepthT
  real(8), parameter :: icedensity = 933.  ! 120K  [kg/m^3]
  integer typeT
  real(8) zdepthTnew, buf, bigdtsec, beta

  if (zdepthT<0.) return   ! no ice anywhere

  bigdtsec = bigstep*86400*365.24

  ! advance ice table
  if (any( icefrac > porosity )) then  ! deflation
     !beta = (1-icefrac)/(1-porosity)  ! still valid?
     stop 'icechanges3: Deflation not implemented'
  else
     beta = 1. ! no deflation
  end if
  if (zdepthT>z(nz)) then
     zdepthT=-9999.
  else
     typeT = getfirst(zdepthT,nz,z)
     if (typeT<1) stop 'typeT<1 in subroutine icechanges3'
     buf = elleff*avSice*beta/(icefrac(typeT)*icedensity)
     zdepthTnew = sqrt(2*buf*bigdtsec + (zdepthT+elleff)**2) - elleff
     print *,'# advance of ice table',zdepthT,zdepthTnew

     zdepthT = zdepthTnew
  end if
  
end subroutine icechanges3



pure function getfirst(zdepth,nz,z)
  ! returns index of shallowest grid point below ice table
  ! if zdepth<0 or zdepth>z(nz) (no ice table) then it returns -9
  implicit none
  integer getfirst
  integer, intent(IN) :: nz
  real(8), intent(IN) :: zdepth, z(nz)
  integer j
  getfirst = -9
  if (zdepth<0.) return
  do j=1,nz
     if (z(j)>zdepth) then
        getfirst = j  
        exit
     endif
  enddo
end function getfirst
