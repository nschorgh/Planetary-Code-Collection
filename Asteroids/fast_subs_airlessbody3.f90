!************************************************************************
! Subroutines for retreat of subsurface ice on airless bodies,
!      based on time-averaged vapor transport equations
!
! written by Norbert Schorghofer 2025-2026
!************************************************************************


subroutine icelayer_asteroid(bigstep,NP,z,porosity,icefrac,Tinit, &
     & zdepthT,Tmean1,Tmean3,Tmin,Tmax,latitude,Orbit,S0)
!************************************************************************
! bigstep = time step [Earth years]
! NP = number of geographic locations (sites)
! z = depths [m]
! zdetphT = depth of ice table or shallowest depth with perennial ice [m]  
! latitude  [degree]
! eps = axis tilt [radians]  
! S0 = solar constant relative to present
!************************************************************************
  use constants, only : d2r
  use body, only : Tnominal, nz, diam, orbitp
  use allinterfaces
  implicit none
  integer, intent(IN) :: NP
  real(8), intent(IN) :: bigstep
  real(8), intent(IN) :: z(nz), porosity(nz), icefrac(nz)
  logical, intent(IN) :: Tinit
  real(8), intent(INOUT) :: zdepthT(NP), Tmean1(NP), Tmean3(NP)
  real(8), intent(OUT) :: Tmin(NP), Tmax(NP)
  real(8), intent(IN) :: latitude(NP), S0
  type(orbitp), intent(IN) :: Orbit
  
  integer k, typeT, j, jump
  real(8) ti(nz), rhocv(nz), T(nz)
  real(8) ell0, elleff, deltaz, avSice
  real(8), dimension(nz) :: ell
  real(8), SAVE :: zdepth_old(100)  ! NP<=100

  do k=1,NP   ! big loop over sites

     typeT = getfirst(zdepthT(k),nz,z)

     ! assign/update property profiles
     if (Tinit) then
        T(:) = spread(Tnominal,1,nz)
     else
        T(:) = ( Tmean1(k)*(z(nz)-z(:)) + Tmean3(k)*z(:) ) / z(nz)
     end if
     call assignthermalproperties3(nz,z(:),T,porosity, &
          &                        ti(:),rhocv(:),icefrac(:),zdepthT(k))
     ell0 = meanfreepathinsoil(diam,porosity(1)) ! surface
     do j=1,nz
        ell(j) = meanfreepathinsoil(diam,porosity(j)) 
        if (z(j)>zdepthT(k) .and. zdepthT(k)>=0.) then
           ell(j) = meanfreepathinsoil(diam,porosity(j)+icefrac(j)) 
        endif
     end do
     
     ! run thermal model
     call ajsub_asteroid(latitude(k)*d2r, z, ti, rhocv, Orbit, S0, &
          &     typeT, avSice, Tinit, Tmean1(k), Tmean3(k), Tmin(k), Tmax(k))

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
     call icechanges3(nz,z(:),avSice,elleff,bigstep,zdepthT(k),porosity,icefrac)

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
          &        bigstep,latitude(k),zdepthT(k),avSice,jump,elleff
     zdepth_old(k) = zdepthT(k)

  end do  ! end of big loop
end subroutine icelayer_asteroid



subroutine ajsub_asteroid(latitude, z, ti, rhocv, Orbit, S0, &
     &     typeT, avSice, Tinit, Tmean1, Tmean3, Tmin, Tmaxi)
!************************************************************************
!  A 1D thermal model that also returns various time-averaged quantities
!
!  Tinit = initalize if .true., otherwise use Tmean1 and Tmean3
!************************************************************************
  use constants
  use body, only : EQUILTIME, dt, Fgeotherm, nz, emiss, albedo, orbitp
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
  real(8), intent(INOUT) :: Tmean1, Tmean3
  real(8), intent(OUT) :: Tmin, Tmaxi
  real(8), parameter :: mmass = 18.  ! H2O
  real(8), parameter :: zero = 0.
  integer nsteps, n, nm
  real(8) tmax, time, Qn, Qnp1, tdays
  real(8) orbitR, orbitLs, orbitDec, HA
  real(8) Tsurf, Fsurf, T(nz)
  real(8) Tmean0, S1, coslat, solsperyear
  
  ! initialize
  solsperyear = sols_per_year( orbit%semia, orbit%solarDay)
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
     tmax = 3*EQUILTIME*solsperyear
  else
     T(:) = ( Tmean1*(z(nz)-z(:)) + Tmean3*z(:) ) / z(nz)
     Tsurf = Tmean1
     tmax = EQUILTIME*solsperyear
  endif
  Fsurf=0.

  nsteps=int(tmax/dt)       ! calculate total number of timesteps

  nm=0
  Tmean1=0.; Tmean3=0.
  avSice = 0.
  Tmin=+1e32; Tmaxi=-9.

  time=0.
  call generalorbit( 0.d0, orbit%semia, orbit%ecc, orbit%omega, orbit%eps, &
       & orbitLs,orbitDec,orbitR)
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

     if (time>=tmax-solsperyear) then
        Tmean1 = Tmean1+Tsurf
        Tmean3 = Tmean3+T(nz)
        if (typeT>0 .and. typeT<=nz) then
           !rhosatav = rhosatav+psv(T(typeT))/T(typeT)
           avSice = avSice + psv(T(typeT)) / sqrt(T(typeT))
        end if
        nm = nm+1

        if (Tsurf<Tmin) Tmin=Tsurf
        if (Tsurf>Tmaxi) Tmaxi=Tsurf
     endif

  end do  ! end of time loop
  
  Tmean1 = Tmean1/nm; Tmean3 = Tmean3/nm
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
