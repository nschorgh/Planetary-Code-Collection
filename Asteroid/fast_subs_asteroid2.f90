!*****************************************************
! Subroutines for fast asteroid method
! written by Norbert Schorghofer 2014-2015
!*****************************************************


subroutine icelayer_asteroid(bigstep,NP,z,porosity,Tinit, &
     & zdepthP,sigma,Tmean1,Tmean3,Tmin,Tmax,latitude,albedo,ecc,omega,eps,S0)
!*************************************************************************
! bigstep = time step [Earth years]
! latitude  [degree]
! S0 = solar constant relative to present (as defined in flux_noatm.f90)
!*************************************************************************
  use constants, only : d2r, NMAX
  use body, only : icedensity, Tnominal, nz
  use allinterfaces, except_this_one => icelayer_asteroid
  implicit none
  integer, intent(IN) :: NP
  real(8), intent(IN) :: bigstep
  real(8), intent(IN) :: z(NMAX), porosity(nz)
  logical, intent(IN) :: Tinit
  real(8), intent(INOUT) :: sigma(nz,NP), zdepthP(NP), Tmean1(NP), Tmean3(NP)
  real(8), intent(OUT) :: Tmin(NP), Tmax(NP)
  real(8), intent(IN) :: latitude(NP), albedo(NP), ecc, omega, eps, S0

  integer k, typeP, j, jump
  real(8) ti(NMAX), rhocv(NMAX), maxsigma, diam
  real(8) Deff, deltaz, Diff0, Jp, avrhotmp
  real(8), dimension(nz) :: Diff, ypp, avrho, porefill
  real(8), SAVE :: zdepth_old(100)  ! NP<=100

  do k=1,NP   ! big loop over sites

     typeP = -9
     do j=1,nz
        if (z(j)>zdepthP(k)) then
           typeP = j 
           exit
        endif
     enddo

     ! assign/update property profiles
     porefill(:) = sigma(:,k)/(porosity(:)*icedensity)
     call assignthermalproperties(nz,z,Tnominal,porosity,ti,rhocv,porefill)
     diam = 100e-6  ! assumed mean-free-path in mantle
     Diff0 = vapordiffusivity(diam,porosity(1),Tnominal) ! surface
     do j=1,nz
        if (z(j)>0.5) diam=1e-3  ! coarser below 0.5m
        Diff(j) = vapordiffusivity(diam,porosity(j),Tnominal) 
        if (sigma(j,k)>0.) then
           porefill(j) = sigma(j,k)/(porosity(j)*icedensity)
           Diff(j) = constriction(porefill(j))*Diff(j)
        endif
     enddo
     
     ! run thermal model
     call ajsub_asteroid(latitude(k)*d2r, albedo(k), z, ti, rhocv, & 
          &     ecc, omega, eps, S0, typeP, Diff(:), Diff0, avrho(:), &
          &     Tinit, ypp, Jp, Tmean1(k), Tmean3(k), Tmin(k), Tmax(k))

     ! run ice evolution model
     if (typeP<=1) then
        Deff = Diff0
     else
        deltaz = colint(spread(1d0,1,nz),z,nz,1,typeP-1)  ! for normalization
        Deff = deltaz/colint(1./Diff,z,nz,1,typeP-1) 
        if (minval(Diff(1:typeP-1))<=0.) then
           print *,typeP,porefill(1:typeP-1),Diff(1:typeP-1)
           stop 'D_EFF PROBLEM'
        endif
     endif
     call impactstirring(nz,z(:),bigstep,sigma(:,k))  ! turn impact stiring on and off here
     call icechanges(nz,z(:),typeP,avrho(:),ypp(:),Deff,bigstep,Jp,zdepthP(k),sigma(:,k))
     where(sigma<0.) sigma=0.
     do j=1,nz
        maxsigma = porosity(j)*icedensity
        if (sigma(j,k)>maxsigma) sigma(j,k)=maxsigma
     enddo

     ! diagnose
     if (zdepthP(k)>=0.) then
        jump = 0
        do j=1,nz
           if (zdepth_old(k)<z(j).and.zdepthP(k)>z(j)) jump=jump+1
        enddo
     else
        jump=-9
     endif
     if (typeP>0) then
        avrhotmp = avrho(typeP)
     else
        avrhotmp = -9999.
     end if
     write(34,'(f12.2,1x,f6.2,1x,f11.5,1x,g11.4,1x,i3,1x,g10.4)') &
          &        bigstep,latitude(k),zdepthP(k),avrhotmp,jump,Deff
     zdepth_old(k) = zdepthP(k)

  enddo  ! end of big loop
end subroutine icelayer_asteroid



subroutine ajsub_asteroid(latitude, albedo, z, ti, rhocv, ecc, omega, eps, &
     &     S0, typeP, Diff, Diff0, rhosatav, Tinit, ypp, Jp, Tmean1, Tmean3, Tmin, Tmaxi)
!***********************************************************************
!  A 1D thermal model that also returns various time-averaged quantities
!
!  Tinit = initalize if .true., otherwise use Tmean1 and Tmean3
!***********************************************************************
  use constants
  use body, only : EQUILTIME, dt, solsperyear, Fgeotherm, semia, nz, emiss, solarDay
  use allinterfaces, except_this_one => ajsub_asteroid
  implicit none
  real(8), intent(IN) :: latitude  ! in radians
  real(8), intent(IN) :: albedo, z(NMAX)
  real(8), intent(IN) :: ti(NMAX), rhocv(NMAX)
  real(8), intent(IN) :: ecc, omega, eps, Diff(nz), Diff0, S0
  integer, intent(IN) :: typeP
  real(8), intent(OUT) :: rhosatav(nz)   ! annual mean vapor density
  logical, intent(IN) :: Tinit
  real(8), intent(OUT) :: ypp(nz), Jp
  real(8), intent(INOUT) :: Tmean1, Tmean3
  real(8), intent(OUT) :: Tmin, Tmaxi
  integer nsteps, n, j, nm
  real(8) tmax, time, Qn, Qnp1, tdays
  real(8) orbitR, orbitLs, orbitDec, HA
  real(8) Tsurf, Fsurf, T(NMAX)
  real(8) Tmean0, rhosatav0, rlow , S1, coslat
  real(8), external :: psv
  
  ! initialize
  if (Tinit) then 
     S1=S0*1365./semia**2  ! must match solar constant defined in flux_noatm
     coslat = max(cos(latitude),cos(latitude+eps),cos(latitude-eps))
     Tmean0 = (S1*(1.-albedo)*coslat/(pi*emiss*sigSB))**0.25 ! estimate
     Tmean0 = Tmean0-5.
     if (Tmean0<50.) Tmean0=50.
     print *,Tmean0,S1,latitude,cos(latitude)
     write(6,*) '# initialized with temperature estimate of',Tmean0,'K'
     write(34,*) '# initialized with temperature estimate of',Tmean0,'K'
     T(1:nz) = Tmean0 
     Tsurf = Tmean0
     tmax = 3*EQUILTIME*solsperyear
  else
     forall(j=1:nz) T(j) = (Tmean1*(z(nz)-z(j))+Tmean3*z(j))/z(nz)
     Tsurf = Tmean1
     tmax = EQUILTIME*solsperyear
  endif
  Fsurf=0.

  nsteps=int(tmax/dt)       ! calculate total number of timesteps

  nm=0
  Tmean1=0.; Tmean3=0.
  rhosatav0 = 0.; rhosatav(:) = 0.
  Tmin=+1e32; Tmaxi=-9.

  time=0.
  call generalorbit(0.d0,semia,ecc,omega,eps,orbitLs,orbitDec,orbitR)
  HA=2.*pi*time             ! hour angle
  Qn=S0*(1-albedo)*flux_noatm(orbitR,orbitDec,latitude,HA,0.d0,0.d0)
  !----loop over time steps 
  do n=0,nsteps-1
     time =(n+1)*dt         !   time at n+1 
     tdays = time*(solarDay/86400.) ! parenthesis may improve roundoff
     call generalorbit(tdays,semia,ecc,omega,eps,orbitLs,orbitDec,orbitR)
     HA=2.*pi*mod(time,1.d0)  ! hour angle
     Qnp1=S0*(1-albedo)*flux_noatm(orbitR,orbitDec,latitude,HA,0.d0,0.d0)
     
     call conductionQ(nz,z,dt*solarDay,Qn,Qnp1,T,ti,rhocv,emiss, &
          &           Tsurf,Fgeotherm,Fsurf)
     Qn=Qnp1
     
     if (time>=tmax-solsperyear) then
        Tmean1 = Tmean1+Tsurf
        Tmean3 = Tmean3+T(nz)
        rhosatav0 = rhosatav0+psv(Tsurf)/Tsurf
        do j=1,nz
           rhosatav(j) = rhosatav(j)+psv(T(j))/T(j)
        enddo
        nm=nm+1

        if (Tsurf<Tmin) Tmin=Tsurf
        if (Tsurf>Tmaxi) Tmaxi=Tsurf
     endif

  enddo  ! end of time loop
  
  Tmean1 = Tmean1/nm; Tmean3 = Tmean3/nm
  rhosatav0 = rhosatav0/nm; rhosatav(:)=rhosatav(:)/nm
  rhosatav0 = rhosatav0*18./8314.; rhosatav(:)=rhosatav(:)*18./8314.

  rlow=rhosatav(nz-1)
  call avmeth(nz, z, rhosatav(:), rhosatav0, rlow, typeP, Diff(:), Diff0, ypp(:), Jp)

  if (typeP<=0) rhosatav = -9999.
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



subroutine avmeth(nz, z, rhosatav, rhosatav0, rlow, typeP, Diff, Diff0, ypp, Jpump1)
!***********************************************************************
!  returns 2nd derivative ypp and pumping flux
!***********************************************************************
  use allinterfaces, except_this_one => avmeth
  implicit none
  integer, intent(IN) :: nz, typeP
  real(8), intent(IN), dimension(nz) :: z, rhosatav, Diff
  real(8), intent(IN) :: rhosatav0, rlow, Diff0
  real(8), intent(OUT) :: ypp(nz), Jpump1
  real(8) yp(nz), ap_one, ap(nz)
  !real(8), external :: deriv1_onesided

!-calculate pumping flux at interface
  call deriv1(z,nz,rhosatav,rhosatav0,rlow,yp)  ! yp also used below
  Jpump1 = -Diff(typeP)*yp(typeP)
  ! yp is always <0

!-calculate ypp
  call deriv1(z,nz,Diff(:),Diff0,Diff(nz-1),ap)
  if (typeP>0 .and. typeP<nz-2) then
     ap_one=deriv1_onesided(typeP,z(:),nz,Diff(:))
     ! print *,typeP,ap(typeP),ap_one
     ap(typeP)=ap_one
  endif
  call deriv2_simple(z,nz,rhosatav(1:nz),rhosatav0,rlow,ypp(:))
  ypp(:) = ap(:)*yp(1:)+Diff(:)*ypp(:)

end subroutine avmeth



subroutine icechanges(nz,z,typeP,avrho,ypp,Deff,bigstep,Jp,zdepthP,sigma)
!***********************************************************
! advances ice interface and grows pore ice
!***********************************************************
  use allinterfaces, except_this_one => icechanges
  implicit none
  integer, intent(IN) :: nz, typeP
  real(8), intent(IN) :: z(nz), ypp(nz), avrho(nz)
  real(8), intent(IN) :: Deff, bigstep, Jp
  real(8), intent(INOUT) :: zdepthP, sigma(nz)
  integer j, newtypeP
  real(8) zdepthPnew, buf, bigdtsec, dtcorr, dtstep, dz(nz)

  if (typeP<0) return   ! no ice anywhere
  if (zdepthP<0.) print *,'Error: No ice in icechanges'
  bigdtsec = bigstep*86400*365.24

  ! loss of gardened-up ice
  dtcorr = 0.
  if (typeP-1<1) goto 20
  if (maxval(sigma(1:typeP-1))==0.) goto 20
  call dzvector(nz,z,dz) 
  do j=1,typeP-1
     dtcorr = dtcorr + sigma(j)*z(j)/(Deff*avrho(j))*dz(j)
     sigma(j) = 0.
     if (dtcorr>bigdtsec) then
        dtcorr=bigdtsec
        write(6,*) '# icechanges: early return',j,typeP-1
        newtypeP = j+1
        goto 30
     endif
  enddo
  print *,'# correction time ratio',dtcorr/bigdtsec
20 continue

  ! advance ice table
  buf = (Deff*avrho(typeP)+zdepthP*Jp)/sigma(typeP)
  !zdepthPnew = sqrt(2*buf*bigdtsec + zdepthP**2)
  zdepthPnew = sqrt(2*buf*(bigdtsec-dtcorr) + zdepthP**2)
  newtypeP = gettypeP(zdepthPnew,nz,z)
  if (newtypeP>typeP+1) then  ! take two half steps
     print *,'# icechanges: half step',typeP,newtypeP,sigma(typeP),sigma(newtypeP),zdepthPnew
     !dtstep = bigdtsec/2.  ! half the time step
     dtstep = (bigdtsec-dtcorr)/2.  ! half the time step

     buf = (Deff*avrho(typeP)+zdepthP*Jp)/sigma(typeP)
     zdepthPnew = sqrt(2*buf*dtstep + zdepthP**2)  ! 1st half
     newtypeP = gettypeP(zdepthPnew,nz,z)
  
     buf = (Deff*avrho(newtypeP)+zdepthPnew*Jp)/sigma(newtypeP)
     zdepthPnew = sqrt(2*buf*dtstep + zdepthPnew**2) ! 2nd half
     newtypeP = gettypeP(zdepthPnew,nz,z)
  endif
  write(6,*) '# advance of ice table',typeP,zdepthP,newtypeP,zdepthPnew

  zdepthP = zdepthPnew
  if (zdepthP>z(nz)) zdepthP=-9999.
  if (newtypeP>1) sigma(1:newtypeP-1)=0.
  
  ! diffusive filling
30 continue
  if (newtypeP>0) then  
     do j=newtypeP,nz
        sigma(j) = sigma(j) + ypp(j)*bigdtsec
     enddo
  end if
  where(sigma<0.) sigma=0.

end subroutine icechanges



pure function zint(y1,y2,z1,z2)
  ! interpolate between two values, y1*y2<0
  implicit none
  real(8), intent(IN) :: y1,y2,z1,z2
  real(8) zint
  zint = (y1*z2 - y2*z1)/(y1-y2)
end function zint



pure function colint(y,z,nz,i1,i2)
  ! column integrates y
  implicit none
  integer, intent(IN) :: nz, i1, i2
  real(8), intent(IN) :: y(nz),z(nz)
  real(8) colint
  integer i
  real(8) dz(nz)
  dz(1)=(z(2)-0.)/2
  do i=2,nz-1
     dz(i) = (z(i+1)-z(i-1))/2.
  enddo
  dz(nz) = z(nz)-z(nz-1)
  colint= sum(y(i1:i2)*dz(i1:i2))
end function colint



subroutine compactoutput(unit,sigma,nz)
  implicit none
  integer, intent(IN) :: unit,nz
  real(8), intent(IN) :: sigma(nz)
  integer j
  do j=1,nz
     if (sigma(j)==0.) then
        write(unit,'(1x,f2.0)',advance='no') sigma(j)
     else
        write(unit,'(1x,f7.3)',advance='no') sigma(j)
     endif
  enddo
  write(unit,"('')")
end subroutine compactoutput



subroutine assignthermalproperties(nz,z,Tnom,porosity,ti,rhocv,porefill)
!*********************************************************
! assign thermal properties of soil
! specify thermal interia profile here
!*********************************************************
  use body, only : icedensity
  use allinterfaces, only : heatcapacity
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), Tnom, porosity(nz)
  real(8), intent(OUT) :: ti(nz), rhocv(nz)
  real(8), intent(IN), optional :: porefill(nz)
  real(8), parameter :: rhodry = 2400  ! bulk density
  real(8), parameter :: kbulk = 2. ! conductivity for zero porosity dry rock
  real(8), parameter :: kice=4.6, cice=1145   ! 140K
  !real(8), parameter :: kice=4.3, cice=1210   ! 150K
  integer j
  real(8) cdry  ! heat capacity of dry regolith
  real(8) k(nz)  ! thermal conductivity
  real(8) thIn

  if (minval(porosity)<0. .or. maxval(porosity)>0.8) then
     print *,'Error: unreasonable porosity',minval(porosity),maxval(porosity)
     stop
  endif

  cdry = heatcapacity(Tnom)
  thIn = 15.
  do j=1,nz
     rhocv(j) = (1.-porosity(j))*rhodry*cdry
     if (z(j)>0.5) thIn=50.
     k(j) = thIn**2/rhocv(j) 
     !if (porosity(j)<phi0) then  ! weighted harmonic mean of k and kbulk
     !   k0 = thIn**2/rhocv(1) 
     !   k(j) = 1./( (porosity(j)/phi0)/k(j) + (1-porosity(j)/phi0)/kbulk)
     !endif
  enddo
  if (present(porefill)) then
     do j=1,nz
        if (porefill(j)>0) then
           k(j) = k(j) + porosity(j)*kice*porefill(j)
           rhocv(j) = rhocv(j) + icedensity*cice*porefill(j)*porosity(j)
        endif
     enddo
  end if

  ti(1:nz) = sqrt(k(1:nz)*rhocv(1:nz))
end subroutine assignthermalproperties



elemental function heatcapacity(T)
  implicit none
  real(8), intent(IN) :: T
  real(8) heatcapacity
  real(8) c
  
  ! heat capacity from Ledlow et al. (1992), <350K
  !c = 0.1812 + 0.1191*(T/300.-1) + 0.0176*(T/300.-1)**2 + &
  !     0.2721*(T/300.-1)**3 + 0.1869*(T/300.-1)**4
  !c = c*1000*4.184  ! cal/(g K) -> J/(kg K)

  ! heat capacity from Winter & Saari (1969),  20K<T<500K
  c = -0.034*T**0.5 + 0.008*T - 0.0002*T**1.5
  heatcapacity = c*1000   ! J/(g K) -> J/(kg K)
end function heatcapacity



elemental function conductivity(T)
  implicit none
  real(8), intent(IN) :: T
  real(8) conductivity
  real(8) A, B

  A = 0.001; B = 3e-11     ! D= 100um based on Sakatani et al. (2012)
  !A = 0.003; B = 1.5e-10   ! D= 1mm based on Sakatani et al. (2012)
  !A = 0.01; B = 2e-9       ! D= 1cm

  conductivity = A + B*T**3
end function conductivity



function vapordiffusivity(diam,porosity,T)
  ! diam = rms grain diameter
  use constants, only : pi, kB
  implicit none
  real(8) vapordiffusivity
  real(8), intent(IN) :: diam,porosity,T
  real(8) vbar, r
  real(8), parameter :: tau = 2.  ! tortuosity

  r = diam/2.
  vbar = sqrt(8*kB*T/(pi*18*1.66e-27))
  ! for 0<=porosity<=0.5
  vapordiffusivity = pi/(8+pi)*porosity/(1-porosity)*vbar*r/tau
end function vapordiffusivity



elemental function constriction(porefill)
! vapor diffusivity constriction function, 0<=eta<=1
  implicit none
  real(8), intent(IN) :: porefill
  real(8) eta, constriction
  if (porefill<=0.) eta = 1.
  if (porefill>0. .and. porefill<1.) then
     ! eta = 1.
     ! eta = 1-porefill
     eta = (1-porefill)**2  ! Hudson et al., JGR, 2009
  endif
  if (porefill>=1.) eta = 0.
  constriction = eta
end function constriction



elemental function faintsun(t)
  implicit none
  real(8) faintsun
  real(8), intent(IN) :: t   ! time before present [years]
  ! Gough, D. O. (1981), Sol. Phys., 74, 21–34
  faintsun = 1./(1+0.4*abs(t)/4.57e9)
end function faintsun



subroutine dzvector(nz,z,dz) 
  ! matches colint
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz)
  real(8), intent(OUT) :: dz(nz)
  integer i
  dz(1)=(z(2)-0.)/2
  do i=2,nz-1
     dz(i) = (z(i+1)-z(i-1))/2.
  enddo
  dz(nz) = z(nz)-z(nz-1)
end subroutine dzvector


integer function gettypeP(zdepthP,nz,z)
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: zdepthP, z(nz)
  integer j
  gettypeP = -9 
  do j=1,nz
     if (z(j)>zdepthP) then
        gettypeP = j  
        exit
     endif
  enddo
end function gettypeP
