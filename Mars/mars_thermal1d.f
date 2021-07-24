      program mars_thermal1d
C***********************************************************************
C   mars_thermal1d:  program to calculate the diffusion of temperature 
C                    into the ground
C   Eqn: rho*c*T_t = (k*T_z)_z   where k is the thermal conductivity
C   BC (z=0): Q(t) + k*T_z = emiss*sig*T^4 + L*dm/dt
C   BC (zmax): geothermal flux
C
C   Grid: surface is at z=0; T(i) is at z(i)
C***********************************************************************
      implicit none
      integer NMAX
      real*8 pi, d2r, zero, earthDay, marsDay
      real*8 sigSB, Tco2frost, Lco2frost
      real*8 solsy   ! number of sols in a Mars year
      parameter (NMAX=1000, sigSB=5.6704d-8)
      parameter (pi=3.1415926535897932, d2r=pi/180., zero=0.)
      parameter (earthDay=86400., marsDay=88775.244, solsy=668.60)
      parameter (Lco2frost=6.0e5)

      integer nz, nsteps, n, i, nm
      integer julday, iyr, imm, iday
      real*8 T(NMAX), tmax, time, dt, zmax, dz, zfac
      real*8 latitude, thermalInertia, albedo, albedo0, emiss, emiss0
      real*8 fracIR, fracDust, rhoc, delta
      real*8 Qn, Qnp1, tdays, dtsec
      real*8 marsR, marsLs, marsDec, HA
      real*8 jd, temp1, dcor, dt0_j2000, flux_mars77
      real*8 ti(NMAX), Tsurf, rhocv(NMAX), z(NMAX), Told(NMAX)
      real*8 co2albedo, Fgeotherm, Tsurfold, Fsurf, Fsurfold, m, dE
      real*8 co2emiss, Tinit, Tmean0, Tmean2, geof, ps, pb, psv
!     real*8 psurf, psurf_season, tfrostco2
      character*100 dum1
      character*40 fileout1, fileout2  ! character arrays for output filenames
      external julday, flux_mars77, psv

C-----read input       
      open(unit=20,file='input.par',status='old')
      read(20,'(a)')dum1
      read(20,*)dt,tmax
      read(20,'(a)')dum1
      read(20,*)nz,zmax,zfac
      read(20,'(a)')dum1
      read(20,*)latitude,albedo0,emiss0
      read(20,'(a)')dum1
      read(20,*)thermalInertia,rhoc
      read(20,'(a)')dum1
      read(20,*)fracIR,fracDust,Fgeotherm
      read(20,'(a)')dum1
      read(20,*)Tco2frost,co2albedo,co2emiss
      read(20,'(a)')dum1
      read(20,*)iyr,imm,iday
      fileout1='Tsurface'; fileout2='Tprofile'
      close(20)

C     set some constants
      dz = zmax/nz              ! set the space step
      nsteps = nint(tmax/dt)     ! calculate total number of timesteps
      emiss = emiss0     ! frost-free emissivity
      delta = thermalInertia/rhoc*sqrt(marsDay/pi)  ! skin depth
      albedo = albedo0
      dtsec = dt*marsDay

      jd = dble(julday(imm,iday,iyr))  !  JD for noon UTC on iyear/imm/iday
      temp1 = (jd-2451545.d0)/36525.d0
      dcor = (64.184d0 + 95.*temp1 + 35.*temp1**2) ! correction in sec
C     All time is referenced to dt0_j2000
      dt0_j2000 = jd + dcor/earthDay - 2451545.d0 
      
      call marsorbit(dt0_j2000,0.d0,marsLs,marsDec,marsR)

      write(*,*) 'Model parameters'
      write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
      write(*,*) 'zmax=',zmax,' zfac=',zfac
      write(*,*) 'Thermal inertia=',thermalInertia,' rho*c=',rhoc
      print *,'Diurnal skin depth=',delta,' Geothermal flux=',Fgeotherm
      write(*,*) 'albedo=',albedo,' emiss=',emiss
      write(*,*) 'fracIR=',fracIR,' fracDust=',fracDust
      write(*,*) 'CO2 frostpoint=',Tco2frost
      write(*,*) 'CO2 albedo=',co2albedo,' CO2 emissivity=',co2emiss
      write(*,*) 'Reference time is noon UTC ',iyr,'/',imm,'/',iday
      write(*,*) 'Julian day ',jd
      write(*,*) 'Initial Mars orbit parameters: '
      write(*,*) 'Ls(deg)=',marsLs/d2r,' declination(deg)=',marsDec/d2r
      write(*,*) 'r (AU)=',marsR
      write(*,*) 'Calculations performed for latitude=',latitude

C-----initialize
      Tinit = 210.15     ! black-body temperature of planet
      geof = cos(latitude*d2r)/pi
      Tinit = (589.*(1.-albedo0)*geof/sigSB)**0.25  ! better estimate
      do i=1,nz
         T(i) = Tinit
         ! z(i) = (i-0.5)*dz
      enddo
      Tsurf = Tinit
      latitude = latitude*d2r
      Tmean0=0.; Tmean2=0.; nm=0
      ps=0.; pb=0.

      call setgrid(nz,z,zmax,zfac)
      if (z(6)>delta) then
         print *,'WARNING: less than 6 points within diurnal skin depth'
      endif
      do i=1,nz
         if (z(i)<delta) cycle
         print *,i-1,' grid points within diurnal skin depth'
         exit
      enddo
      if (z(1)<1.e-5) print *,'WARNING: first grid point is too shallow'
      call smartgrid(nz,z,0.1d0,thermalInertia,rhoc,0.4d0,ti,rhocv,
     &     1,zero)
!      call smartgrid(nz,z,0.1d0,thermalInertia,rhoc,zero,ti,rhocv,
!     &     3,zero)
      open(unit=30,file='z',status='unknown')
      write(30,*) (z(i),i=1,nz)
      close(30)
      do i=1,nz ! assign thermal properties, if not already set by smartgrid
!         ti(i) = thermalInertia
!         rhocv(i) = rhoc
!         if (z(i)>0.1) then
!            call soilthprop(0.4d0,1.d0,rhoc,thermalInertia,1,
!     &           rhocv(i),ti(i))
!            call soilthprop(zero,zero,zero,zero,3,rhocv(i),ti(i),zero)
!         endif
      enddo

      time=0.
      m=0.; Fsurf=0.
      open(unit=21,file=fileout1,status='unknown') ! surface temperature
      open(unit=22,file=fileout2,status='unknown') ! temperature profile
      HA = 2.*pi*time   ! hour angle
C     net flux: solar insolation + IR
      Qn = flux_mars77(marsR,marsDec,latitude,HA,albedo,fracir,fracdust)

C-----loop over time steps 
      do n=0,nsteps-1
         time = (n+1)*dt   !   time at n+1 
C        Solar insolation and IR at future time step
         tdays = time*(marsDay/earthDay)  ! parenthesis may improve roundoff
         call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR) 
         HA = 2.*pi*mod(time,1.d0)    ! hour angle
!        psurf = psurf_season(520d0,marsLs) ! seasonally varying surf. pressure
!        Tco2frost = tfrostco2(psurf)
         Qnp1 = flux_mars77(marsR,marsDec,latitude,HA,albedo,
     &        fracir,fracdust)
         Tsurfold = Tsurf
         Fsurfold = Fsurf

         do i=1,nz; Told(i)=T(i); enddo
         if (Tsurf>Tco2frost .or. m<=0.) then   
            call conductionQ(nz,z,dtsec,Qn,Qnp1,T,ti,rhocv,emiss,
     &           Tsurf,Fgeotherm,Fsurf)
         endif

         if (Tsurf<Tco2frost .or. m>0.) then   ! CO2 condensation
            do i=1,nz; T(i)=Told(i); enddo
            call conductionT(nz,z,dtsec,T,Tsurfold,Tco2frost,ti,
     &              rhocv,Fgeotherm,Fsurf) 
            Tsurf = Tco2frost
            dE = (- Qn - Qnp1 + Fsurfold + Fsurf +
     &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
            m = m + dtsec*dE/Lco2frost
         endif
         if (Tsurf>Tco2frost .or. m<=0.) then
            albedo = albedo0
            emiss = emiss0
         else
            albedo = co2albedo
            emiss = co2emiss
         endif

         Qn = Qnp1

c--------only output and diagnostics below this line
         if (time>=tmax-solsy) then
            Tmean0 = Tmean0 + Tsurf
            Tmean2 = Tmean2 + T(nz)
            ps = ps + min(psv(Tsurf),0.12d0)/Tsurf
            pb = pb + psv(T(nz))/T(nz)
            nm = nm+1
         endif
         if (time>=tmax-solsy .and. mod(n,5)==0) then
            write(21,'(f12.6,9f10.4)') time,marsLs/d2r,Tsurf,m
         endif
         if (time>=tmax-solsy .and. mod(n,nint(solsy/dt/40.))==0) then
            write(22,'(f12.4,1000(1x,f6.2))') time,(T(i),i=1,nz)
         endif

      enddo                     ! end the time loop

      close(21)
      close(22)

      print *,Tmean0/nm,Tmean2/nm,ps/nm,pb/nm

      end
 
