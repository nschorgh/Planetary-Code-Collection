subroutine jsubv(NS, zdepth, latitude, albedo0, thIn, pfrost, nz, &
     &     rhoc, fracIR, fracDust, Fgeotherm, dt, zfac, icefrac, & 
     &     surfaceSlope, azFac, mode, avdrho, Tb)
!***********************************************************************
!  jsubv: runs thermal model for Mars at several sites simultaneously
!         allows for thermal emission of slopes to other slopes
!
!   Eqn: T_t = (D*T_z)_z   where D=k/(rho*c)
!   BC (z=0): Q(t) + k*T_z = em*sig*T^4 + L*dm/dt
!   BC (z=L): fixed flux, Fgeotherm
!
!  OUTPUT: 
!           avdrho = difference in annual mean vapor density
!
!  IN- or OUTPUT: 
!           Tb = temperature at the bottom of the domain at
!                a particular time of year
!           when positive, Tb is used for initialization
!
!  mode = 0  ends after fixed amount of time and writes all data
!  mode = 1  ends when certain accuracy is reached and writes no data
!
!  latitude, surfaceSlope, and azFac  must be in RADIANS
!
!  Grid: surface is at z=0; T(i) is at z(i)
!***********************************************************************

  implicit none
  integer, parameter :: NMAX=1000
  real*8, parameter :: pi=3.1415926535897932, d2r=pi/180., zero=0.
  real*8, parameter :: earthDay=86400., marsDay=88775.244, solsperyear=668.60
  real*8, parameter :: sigSB=5.6704d-8, Lco2frost=6.0e5, co2albedo=0.65, co2emiss=1.

  integer, intent(IN) :: NS, nz, mode
  real*8, intent(IN) :: zdepth(NS), latitude, albedo0, thIn, pfrost, rhoc
  real*8, intent(IN) :: fracIR, fracDust, Fgeotherm, dt, zfac, icefrac
  real*8, intent(IN) :: surfaceSlope(NS), azFac(NS)
  real*8, intent(OUT) :: avdrho(NS)
  real*8, intent(INOUT) :: Tb(NS)

  logical outf
  integer nsteps, n, i, nm, k, i0(NS)
  integer julday, iyr, imm, iday

  real*8 T(NMAX,NS),tmax, time, zmax, emiss(NS), albedo(NS)
  real*8 Qn(NS), Qnp1(NS), tdays, dtsec
  real*8 marsR, marsLs, marsDec, HA
  real*8 jd, temp1, dcor, dt0_j2000
  real*8 Tsurf(NS), z(NMAX,NS), ti(NMAX,NS), rhocv(NMAX,NS), Told(NMAX)
  real*8 Tsurfold, Fsurf(NS), Fsurfold, m(NS), dE, Tco2frost
  real*8 Tmean1(NS), Tmean2(NS), rhoavs(NS), rhoavb(NS), Tbold(NS)
  real*8 Qmean(NS), Qland(NS)
  real*8 marsLsold, psv, tfrostco2, emiss0
  real*8 Qdir, Qscat, Qlw, skyviewfactor(NS), Qdir1
  external julday, psv, tfrostco2

  outf = .false.  ! additional output

  select case (mode)
  case (0) ! full mode
     tmax = 7020.
  case (1) ! equilibrate fast
     tmax = 50.*solsperyear  ! use integer multiples of solsperyear
  case default
     stop 'no valid mode'
  end select
  
  iyr=1996; imm=10; iday=5  ! starting year and date
  nsteps = int(tmax/dt)     ! calculate total number of timesteps
  emiss0 = 1.  ! frost-free emissivity
  emiss(:) = emiss0

  zmax = 5.*26.*thIn/rhoc*sqrt(marsDay/pi)
  ! print *,25.86*thIn/rhoc*sqrt(marsDay/pi)
  if (latitude>=0.) then    ! north
     Tco2frost = 147. 
  else                      ! south
     Tco2frost = 143.
  endif
  !Tco2frost = tfrostco2(patm)
  dtsec = dt*marsDay

  if (minval(Tb(:))<=0.) then
     !Tmean1=210.15              ! black-body temperature of planet
     Tmean1 = (589.*(1.-albedo0)*cos(latitude)/pi/5.67e-8)**0.25 ! better estimate
  else
     Tmean1(:) = Tb(:)
  endif
  
  jd = dble(julday(imm,iday,iyr)) ! JD for noon UTC on iyear/imm/iday
  temp1 = (jd-2451545.d0)/36525.d0
  dcor = (64.184d0 + 95.*temp1 + 35.*temp1**2) ! correction in sec
! All time is referenced to dt0_j2000
  dt0_j2000 = jd + dcor/earthDay - 2451545.d0 

  albedo(:) = albedo0
  ti(1:nz,:) = thIn
  rhocv(1:nz,:) = rhoc
  do k=1,NS
     T(1:nz,k) = Tmean1(k)
  enddo
  where (T(1:nz,:) < Tco2frost) T=Tco2frost
  Tsurf(:) = Tmean1(:)
  m(:)=0.; Fsurf(:)=0.
  marsLsold=-1.e32
  Tbold(:)=-1.e32

  Tmean1(:)=0.; Tmean2(:)=0.; nm=0   
  rhoavs(:)=0.; rhoavb(:)=0.
  Qmean(:)=0.; Qland(:)=0.

  do k=1,NS
     call setgrid(nz,z(:,k),zmax,zfac)
     call smartgrid(nz,z(:,k),zdepth(k),thIn,rhoc,icefrac,ti(:,k),rhocv(:,k),1,zero)
     if (mode==0 .and. outf) write(30,'(999(f8.5,1x))') z(1:nz,k)
     if (zdepth(k)<=0 .or. zdepth(k)>=z(nz,k)) then  ! no ice
        i0(k) = nz
     else  ! find first grid point in ice
        do i=1,nz
           if (z(i,k)>zdepth(k)) then ! ice
              i0(k) = i
              exit
           endif
        enddo
     endif
  enddo

  time=0.
  tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
  call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR); 
  HA = 2.*pi*time           ! hour angle

  do k=1,NS
     !Qn(k) = flux(marsR,marsDec,latitude,HA,albedo(k),fracir,fracdust, &
     !     & surfaceSlope(k),azFac(k))
     call flux_mars2(marsR,marsDec,latitude,HA,fracir,fracdust, &
          &          surfaceSlope(k),azFac(k),zero,Qdir,Qscat,Qlw)
     skyviewfactor(k) = cos(surfaceSlope(k)/2.)**2
     Qn(k) = (1-albedo(k))*(Qdir+Qscat*skyviewfactor(k)) + emiss(k)*Qlw*skyviewfactor(k)
     ! neglect terrain irradiances at initialization
  enddo

!-loop over time steps 
  do n=0,nsteps-1
     time = (n+1)*dt        !   time at n+1 
     tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
     call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR)
     HA = 2.*pi*mod(time,1.d0) ! hour angle

     ! k=1 (horizontal unobstructed slope)
     call flux_mars2(marsR,marsDec,latitude,HA,fracir,fracdust, &
          &          surfaceSlope(1),azFac(1),zero,Qdir1,Qscat,Qlw)
     Qnp1(1) = (1-albedo(1))*(Qdir1+Qscat) + emiss(1)*Qlw
     do k=2,NS
        ! direct and sky irradiance
        !Qnp1(k) = flux(marsR,marsDec,latitude,HA,albedo(k),fracir,fracdust,surfaceSlope(k),azFac(k))
        call flux_mars2(marsR,marsDec,latitude,HA,fracir,fracdust, &
             &          surfaceSlope(k),azFac(k),zero,Qdir,Qscat,Qlw)
        Qnp1(k) = (1-albedo(k))*(Qdir+Qscat*skyviewfactor(k)) + emiss(k)*Qlw*skyviewfactor(k)

        ! terrain irradiance
        Qnp1(k) = Qnp1(k) + (1.-skyviewfactor(k))*sigSB*emiss(1)*Tsurf(1)**4
        !Qnp1(k) = Qnp1(k) + (1.-skyviewfactor(k))*(1.-albedo(k))*albedo(1)*Qdir1
     enddo

     do k=1,NS
        Tsurfold = Tsurf(k)
        Fsurfold = Fsurf(k)     
        Told(1:nz) = T(1:nz,k)
        if (m(k)<=0.) then
           call conductionQ(nz,z(:,k),dtsec,Qn(k),Qnp1(k),T(:,k),ti(:,k), &
                & rhocv(:,k),emiss(k),Tsurf(k),Fgeotherm,Fsurf(k))
        endif
        if (Tsurf(k)<Tco2frost .or. m(k)>0.) then   ! CO2 condensation
           T(1:nz,k) = Told(1:nz)
           call conductionT(nz,z(:,k),dtsec,T(:,k),Tsurfold,Tco2frost, &
                & ti(:,k),rhocv(:,k),Fgeotherm,Fsurf(k))
           Tsurf(k) = Tco2frost
           dE = (- Qn(k) - Qnp1(k) + Fsurfold + Fsurf(k) + &
                & emiss(k)*sigSB*(Tsurfold**4+Tsurf(k)**4))/2.
           m(k) = m(k) + dtsec*dE/Lco2frost
        endif
        if (Tsurf(k)>Tco2frost .or. m(k)<=0.) then
           albedo(k) = albedo0
           emiss(k) = emiss0
        else
           albedo(k) = co2albedo
           emiss(k) = co2emiss
        endif
     enddo

     Qn(:) = Qnp1(:)

     if (mode==0 .and. time>=tmax-solsperyear) then
        Tmean1(:) = Tmean1(:)+Tsurf(:)
        Tmean2(:) = Tmean2(:)+T(nz,:)
        do k=1,NS
           rhoavs(k) = rhoavs(k) + min(psv(Tsurf(k)),pfrost) / Tsurf(k)
           rhoavb(k) = rhoavb(k) + psv(T(i0(k),k)) / T(i0(k),k)
        enddo
        Qmean(:) = Qmean(:) + Qn(:)
        Qland(:) = Qland(:) + (1-skyviewfactor(:))*sigSB*emiss(1)*Tsurf(1)**4
        !Qland(:) = Qland(:) + (1-skyviewfactor(:))*(1-albedo(:))*albedo(1)*Qdir1
        nm=nm+1
     endif

     if (mode==1) then
        if (marsLs<marsLsold) then
           Tmean2(:) = Tmean2(:)/nm
           if (maxval(abs(Tmean2(:)-Tbold(:)))<0.05) exit
           Tbold(:) = Tmean2(:)
           Tmean2=0.; nm=0
        else
           Tmean2(:) = Tmean2(:)+T(nz,:)
           nm = nm+1
        endif
     endif
     marsLsold = marsLs
     
  enddo  ! end of time loop
  
  if (mode==0) then
     rhoavs(:) = rhoavs(:)/nm
     rhoavb(:) = rhoavb(:)/nm
     Tmean1(:) = Tmean1(:)/nm
     Tmean2(:) = Tmean2(:)/nm
     avdrho(:) = rhoavb(:)-rhoavs(:)
     Qmean(:) = Qmean(:)/nm
     Qland(:) = Qland(:)/nm

     if (outf) then
        do k=1,NS
           write(34,'(2(f6.2,1x),2(f6.2,1x),2(g10.4,1x),2(g10.4,1x))')  &
                & surfaceSlope(k)/d2r,azFac(k)/d2r,Tmean1(k),Tmean2(k),rhoavb(k),rhoavs(k) &
                & ,Qmean(k),Qland(k)
        enddo
     endif
  endif
  Tb(:) = T(nz,:)
end subroutine jsubv

