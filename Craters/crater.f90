program crater
!***********************************************************************
!   crater: program to calculate the conduction of temperature 
!           into the ground for Mars orbit (as mars_thermal1d) at
!           several points simultaneously
!   Eqn: T_t = D*T_zz   where D=k/(rho*c)
!   BC (z=0): Q(t) + kT_z = em*sig*T^4 + L*dm/dt
!   BC (z=L): no flux
!
!   Skin depth delta =(I/(rho*c))*sqrt(P/pi)
!   where I=thermal inertia [J m^-2 K^-1 s^-1/2], rho=density [kg m^-3],
!   c=specific heat [J K^-1 kg^-1], P=period [s]
!   Grid: surface is at z=0; T(i) is at z(i)
!***********************************************************************

  implicit none
  integer, parameter :: NMAX=3500, NS=2
  real*8 pi, d2r, earthDay, marsDay
  real*8 solsy   ! number of sols in a Mars year
  parameter (pi=3.1415926535897932, d2r=pi/180., &
       &     earthDay=86400., marsDay=88775.244, solsy=668.60)
  real*8, parameter :: sigSB=5.67051d-8, Lco2frost=6.0e5
  
  integer nz, nsteps, n, i, nm, k
  real*8 T(NMAX,NS),tmax, time, dt, zmax, dz
  real*8 latitude, thermalInertia, albedo(NS), emiss(NS)
  real*8 fracIR, fracDust, rhoc, delta, fracIR0
  real*8 Qn(NS), Qnp1(NS), tdays, dtsec
  real*8 marsR, marsLs, marsDec, HA
  real*8 surfaceSlope(NS), azFac(NS), zfac
  real*8 jd, temp1, dcor, dt0_j2000, flux
  real*8 ti(NMAX), Tsurf(NS), rhocv(NMAX), z(NMAX), albedo0, Told(NMAX)
  real*8 co2albedo, Fgeotherm, Tsurfold(NS), Fsurf(NS), Fsurfold(NS), m(NS), dE
  real*8 co2emiss, Tmean, Tmean1(NS), Tmean2(NS), geof, Tco2frost !,a(NS,NS)
  real*8 :: psurf, tfrostco2, rhos(NS), rhob(NS), psv, ps=520
  integer julday, iyr, imm, iday
  external julday, flux, psurf, tfrostco2, psv !, iratmosphere

  !-------inputs
!  dt=0.02; tmax=10000.00
!  nz=80; zmax=4.0; zfac=1.05
!  fracIR0=0.04;  fracDust=0.02
!  thermalInertia=250.;  albedo0=0.20
!  latitude=+38; rhoc=1000000.
!  Fgeotherm=0.028; co2albedo=0.65; co2emiss=1.

  dt=0.02; tmax=7000.
  nz=80; zmax=6.; zfac=1.05
  fracIR0=0.04;  fracDust=0.0
  latitude = -1.95; thermalInertia=200.; albedo0=0.12
  !latitude = -14.57; thermalInertia=315.; albedo0=0.23
  rhoc=1000000.
  Fgeotherm=0.0; co2albedo=0.65; co2emiss=1.

  iyr=1996; imm=10; iday=5  ! starting year and date
  ! azimuth in degrees east of north, 0=north facing
  surfaceSlope(1)= 0.0; azFac(1)=0.0
  !surfaceSlope(2)=14.; azFac(2)= 180.0
  !surfaceSlope(3)=14.; azFac(3)= 200.0
  !surfaceSlope(2)=16; azFac(2)=180.+40.
  !surfaceSlope(3)=16; azFac(3)=180.-41.
  !forall(i=2:NS) azFac(i)=180.0+(i-2)*5.
  !forall(i=2:NS) surfaceSlope(i)=(i-1)*10; azFac(2:NS)=0.
  !surfaceSlope(2:12)=14.0; 
  !azFac(3)=190.; azFac(4)=200.; azFac(5)=210.; azFac(6)=220.; azFac(7)=225.;
  !azFac(8)=170.; azFac(9)=160.; azFac(10)=150.; azFac(11)=140.; azFac(12)=135.;
  surfaceSlope(2)=20.; azFac(2)= 180.0 

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  emiss(:) = 1. ! emissivity
  delta = thermalInertia/rhoc*sqrt(marsDay/pi)  ! skin depth
  albedo(:) = albedo0
  dtsec = dt*marsDay
  
  jd=dble(julday(imm,iday,iyr))  !  JD for noon UTC on iyear/imm/iday
  temp1 = (jd-2451545.d0)/36525.d0
  dcor = (64.184d0 + 95.*temp1 + 35.*temp1**2) ! correction in sec
  !     All time is referenced to dt0_j2000
  dt0_j2000 = jd + dcor/earthDay - 2451545.d0 
  
  call marsorbit(dt0_j2000,0.d0,marsLs,marsDec,marsR)

  write(*,*) 'Model parameters'
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'zmax=',zmax,' zfac=',zfac
  write(*,*) 'Thermal inertia=',thermalInertia,' rho*c=',rhoc
  print *,'Diurnal skin depth=',delta,' Geothermal flux=',Fgeotherm
  write(*,*) 'albedo=',albedo0,' fracIR=',fracIR0,' fracDust=',fracDust
  write(*,*) 'Reference time is noon UTC ',iyr,'/',imm,'/',iday
  write(*,*) 'Julian day ',jd
  write(*,*) 'Initial Mars orbit parameters: '
  write(*,*) 'Ls(deg)=',marsLs/d2r,' declination(deg)=',marsDec/d2r
  write(*,*) 'r (AU)=',marsR
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Surface slope: ',surfaceSlope,' surface azimuth: ',azFac
  write(*,*) 'CO2 albedo=',co2albedo,' CO2 emissivity=',co2emiss
!  write(*,*) 'CO2 frost point',Tco2frost

!-----Initialize
  Tmean = 210.15     ! black-body temperature of planet
  geof = cos(latitude*d2r)/pi
  Tmean=(589.*(1.-albedo0)*geof/5.67e-8)**0.25 ! better estimate
  T(1:nz,:) = Tmean
  ti(1:nz) = thermalInertia
  rhocv(1:nz) = rhoc
  Tsurf(:) = Tmean
  latitude=latitude*d2r
  surfaceSlope(:)=surfaceSlope(:)*d2r
  azFac(:)=azFac(:)*d2r
  Tmean1(:)=0.; Tmean2(:)=0.; nm=0   
  rhos(:)=0.; rhob(:)=0.
  fracIR = fracIR0

! set grid
  do i=1,nz
     z(i) = (i-0.5)*zmax/nz
  enddo
  if (zfac>1.) then
     dz=zmax/(3.+2.*zfac*(zfac**(nz-2)-1.)/(zfac-1.))
     z(1)=dz; z(2)=3*z(1)
     do i=3,nz
        z(i)=(1.+zfac)*z(i-1)-zfac*z(i-2)
     enddo
  endif
  do i=1,nz
     if (z(i)<delta) cycle
     print *,i-1,' grid points within diurnal skin depth'
     exit
  enddo
  if (z(1)<1.e-5) print *,'WARNING: first grid point is too shallow'
  open(unit=30,file='z',status='unknown',action='write');
  write(30,*) (z(i),i=1,nz)
  close(30)
  
  time=0.
  m(:)=0.; Fsurf(:)=0.
  !aa(:,:)=0.
  !aa(2:NS,1)=sin(surfaceSlope(2:NS)/2.)**2
  open(unit=21,file='Tsurface',status='unknown',action='write')
  open(unit=22,file='Tprofile',status='unknown',action='write')
  HA=2.*pi*time   ! hour angle
  ! net flux: solar insolation + IR
  do k=1,NS
     Qn(k)=flux(marsR,marsDec,latitude,HA,albedo(k),fracir,fracdust, &
          & surfaceSlope(k),azFac(k))
  enddo

!-------loop over time steps 
  do n=0,nsteps-1
     time =(n+1)*dt   !   time at n+1 
     ! Solar insolation and IR at future time step
     tdays = time*(marsDay/earthDay)  ! parenthesis may improve roundoff
     call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR); 
     HA=2.*pi*mod(time,1.d0)    ! hour angle
     do k=1,NS
        Qnp1(k)=flux(marsR,marsDec,latitude,HA,albedo(k),fracir,fracdust, &
             & surfaceSlope(k),azFac(k))
        !Qnp1(k)=Qnp1(k) + sigSB*sum(aa(k,:)*emiss(:)*Tsurf(:)**4)
        Qnp1(k)=Qnp1(k) + sigSB*sin(surfaceSlope(k)/2.)**2*emiss(1)*Tsurf(1)**4
     enddo
     Tsurfold(:)=Tsurf(:)
     Fsurfold(:)=Fsurf(:)
     !ps=psurf(marsLs/d2r)
     Tco2frost=tfrostco2(ps)
     fracIR = fracIR0*ps/520.
     do k=1,NS
        Told(1:nz)=T(1:nz,k)
        call conductionQ(nz,z,dtsec,Qn(k),Qnp1(k),T(:,k),ti,rhocv, &
             & emiss(k),Tsurf(k),Fgeotherm,Fsurf(k))
        if (Tsurf(k)<Tco2frost.or.m(k)>0.) then   ! CO2 condensation
           T(1:nz,k)=Told(1:nz)
           call conductionT(nz,z,dtsec,T(:,k),Tsurfold(k),Tco2frost,ti,rhocv,Fgeotherm,Fsurf(k)) 
           Tsurf(k)=Tco2frost
           dE = (- Qn(k) - Qnp1(k) + Fsurfold(k) + Fsurf(k) + &
                & emiss(k)*sigSB*(Tsurfold(k)**4+Tsurf(k)**4))/2.
           m(k) = m(k) + dtsec*dE/Lco2frost;
        endif
        if (Tsurf(k)>Tco2frost.or.m(k)<=0.) then
           albedo(k) = albedo0
           emiss(k) = 1.
        else
           albedo(k) = co2albedo
           emiss(k) = co2emiss
        endif
     enddo
     
     Qn(:)=Qnp1(:)

!--------only output and diagnostics below this line
     if (time>=tmax-solsy) then
        Tmean1(:)=Tmean1(:)+Tsurf(:)
        Tmean2(:)=Tmean2(:)+T(nz,:)
        nm=nm+1
        do k=1,NS
           rhos(k) = rhos(k) + min(psv(Tsurf(k)),0.1202d0)/Tsurf(k)
           rhob(k) = rhob(k) + psv(T(nz,k))/T(nz,k)
        enddo
     endif
     if (time>=tmax-solsy.and.mod(n,2)==0) then
        ! if (mod(n,max(1,nsteps/100000))==0) then
        ! if ((HA>1.*2.*pi/24.and.HA<3.*2.*pi/24.).or.(HA>13.*2.*pi/24.and.HA<15.*2.*pi/24)) then
        write(21,'(f12.6,999f10.4)') time,marsLs/d2r,Tsurf(:),m(:) !,ps
        write(40,'(f12.6,999f10.4)') time,marsLs/d2r,HA,marsDec/d2r,Qn(1),Tsurf(1)
        !Qir=iratmosphere(marsR,marsDec,latitude,fracir,fracdust,surfaceSlope(1))
        !print *,time,Qnp1(1),Qir,0.04*emiss(1)*sigSB*Tsurf(1)**4
     endif
     if (time>=tmax-solsy.and.mod(n,nint(solsy/dt/32.))==0) then
        write(22,'(f12.4,1000(1x,f6.2))') time,(T(i,1),i=1,nz)
        !write(23,'(f12.4,1000(1x,f6.2))') time,(T(i,2),i=1,nz)
     endif
     
  enddo                     ! end the loop

  close(21)
  close(22)
  
  Tmean1=Tmean1/nm; Tmean2=Tmean2/nm; rhob=rhob/nm; rhos=rhos/nm

  !print *,Tmean1(:),Tmean2(:)
  do k=1,NS
     print *,surfaceSlope(k)/d2r,Tmean1(k),Tmean2(k),rhob(k),rhos(k)
  enddo
end program
 

