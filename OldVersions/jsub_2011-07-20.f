      subroutine jsub(zdepth, latitude, albedo0, thIn, pfrost, nz,
     &     rhoc, fracIR, fracDust, Fgeotherm, dt, zfac, icefrac, 
     &     surfaceSlope, azFac, mode, avdrho, Tb, patm)
C***********************************************************************
C  jsub: runs thermal model for Mars and returns difference in mean 
C        annual vapor density between surface and ice at depth zdepth
C
C  OUTPUTS: 
C           avdrho = difference in annual mean vapor density
C
C  IN- or OUTPUT: 
C                 Tb = temperature at the bottom of the domain at
C                      a particular time of year
C
C  mode = 0  ends after fixed amount of time and then outputs all data
C  mode = 1  ends when certain accuracy is reached and outputs no data
C
C  latitude, surfaceSlope, and azFac  must be in RADIANS
C  thermal emission from slope is implemented by changing emissivity
C***********************************************************************

      implicit none
      integer NMAX
      real*8 pi, earthDay, marsDay, sigSB, Lco2frost, solsperyear
      real*8 co2albedo, co2emiss
      parameter (NMAX=2000)
      parameter (pi=3.1415926535897932, 
     &     earthDay=86400.,marsDay=88775.244,solsperyear=668.60)
      parameter (sigSB=5.67051d-8, Lco2frost=6.0e5)
      parameter (co2albedo=0.65, co2emiss=1.)

      integer nz, nsteps, n, i, nm, mode, din
      integer julday, iyr, imm, iday
      real*8 T(NMAX),tmax, time, dt, zmax, surfaceSlope, azFac
      real*8 latitude,thIn,albedo,emiss0,emiss
      real*8 fracIR,fracDust,rhoc,Qn,Qnp1,tdays
      real*8 marsR,marsLs,marsDec,HA
      real*8 jd,temp1,dcor,dt0_j2000,flux
      real*8 ti(NMAX),Tsurf,rhocv(NMAX),z(NMAX),albedo0,zfac
      real*8 Fsurf, m, dE, Told(NMAX), Fgeotherm, Fsurfold, Tsurfold
      real*8 Tpeak(NMAX), Tmean1, Tlow(NMAX), zdepth, Tco2frost, pfrost
      real*8 rhosatav(NMAX), rhoavs, psv, icefrac, avdrho, Tb, Tmean2
      real*8 Tbold, marsLsold, oldtime
      real*8 dirho(NMAX), dirholowz, dirholowLs(NMAX)
      real*8 patm, NULL, tfrostco2
      external julday, flux, psv, tfrostco2
      parameter (NULL=0.)

      if (mode==0) then  ! full mode
         tmax = 7020.
      elseif (mode==1) then  ! equilibrate fast
         tmax = 40.*solsperyear  ! use integer multiple of solsperyear
      else
         stop 'no valid mode'
      endif

      iyr=1996; imm=10; iday=5
      nsteps=int(tmax/dt)       ! calculate total number of timesteps
      emiss0 = 1.*cos(surfaceSlope/2.)**2   ! emissivity
      zmax=5.*26.*thIn/rhoc*sqrt(marsDay/pi)
      print *,thIn/rhoc*sqrt(marsDay/pi)
      print *,thIn/rhoc*sqrt(solsperyear*marsDay/pi)

      !if (latitude>=0.) then    ! north
      !   Tco2frost=147.
      !else                      ! south
      !   Tco2frost=143.
      !endif
      Tco2frost = tfrostco2(patm)

      if (Tb<=0.) then
         !Tmean2=210.15              ! black-body temperature of planet
         Tmean2=(589.*(1.-albedo0)*cos(latitude)/pi/5.67e-8)**0.25 ! better estimate
      else
         Tmean2=Tb
      endif

      jd=dble(julday(imm,iday,iyr)) !  JD for noon UTC on iyear/imm/iday
      temp1 = (jd-2451545.d0)/36525.d0
      dcor = (64.184d0 + 95.*temp1 + 35.*temp1**2) ! correction in sec
C     All time is referenced to dt0_j2000
      dt0_j2000 = jd + dcor/earthDay - 2451545.d0 
      !call marsorbit(dt0_j2000,0.d0,marsLs,marsDec,marsR)

      albedo=albedo0
      emiss = emiss0
      do i=1,nz
         T(i) = Tmean2
         if (T(i)<Tco2frost) T(i)=Tco2frost
         ti(i) = thIn
         rhocv(i) = rhoc
         Tpeak(i) = -1.e32
         Tlow(i) = +1.e32
         rhosatav(i) = 0.
         dirholowLs(i) = 1.e32
      enddo
      Tsurf = T(1)
      m=0.; Fsurf=0.
      nm=0; Tmean1=0.
      rhoavs=0.
      marsLsold=-1.e32; Tmean2=0.
      Tbold=-1.e32
      oldtime=1.e32
      din=-1
      
      call setgrid(nz,z,zmax,zfac)
      call smartgrid(nz,z,zdepth,thIn,rhoc,icefrac,ti,rhocv,1,NULL)
      !call smartgrid(nz,z,zdepth,thIn,rhoc,0.5d0,ti,rhocv,2,icefrac)
      !call smartgrid(nz,z,zdepth,thIn,rhoc,NULL,ti,rhocv,3,NULL)
      if (mode==0) write(30,'(999(f8.5,1x))') (z(i),i=1,nz)

      time=0.
      tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
      call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR)
      HA=2.*pi*time             ! hour angle
      Qn=flux(marsR,marsDec,latitude,HA,albedo,
     &     fracir,fracdust,surfaceSlope,azFac)
C-----loop over time steps 
      do n=0,nsteps-1
         time =(n+1)*dt         !   time at n+1 
         tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
         call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR) 
         HA=2.*pi*mod(time,1.d0)  ! hour angle
         Qnp1=flux(marsR,marsDec,latitude,HA,albedo,
     &        fracir,fracdust,surfaceSlope,azFac)

         Tsurfold=Tsurf
         Fsurfold=Fsurf
         do i=1,nz; Told(i)=T(i); enddo
         if (Tsurf>Tco2frost.or.m<=0.) then
            call conductionQ(nz,z,dt*marsDay,Qn,Qnp1,T,ti,rhocv,emiss,
     &           Tsurf,Fgeotherm,Fsurf)
         endif
         if (Tsurf<Tco2frost.or.m>0.) then ! CO2 condensation
            do i=1,nz; T(i)=Told(i); enddo
            call conductionT(nz,z,dt*marsDay,T,Tsurfold,Tco2frost,ti,
     &              rhocv,Fgeotherm,Fsurf) 
            Tsurf=Tco2frost
            dE = (- Qn - Qnp1 + Fsurfold + Fsurf +
     &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
            m = m + dt*marsDay*dE/Lco2frost
         endif
         if (Tsurf>Tco2frost.or.m<=0.) then
            albedo=albedo0
            emiss =emiss0
         else
            albedo=co2albedo
            emiss =co2emiss*cos(surfaceSlope/2.)**2
         endif
         Qn=Qnp1

         if (time>=tmax-solsperyear.and.mode==0) then
            do i=1,nz
               if (T(i)>Tpeak(i)) Tpeak(i)=T(i)
               if (T(i)<Tlow(i)) Tlow(i)=T(i)
               rhosatav(i)=rhosatav(i)+psv(T(i))/T(i)
            enddo
            Tmean1=Tmean1+Tsurf
            rhoavs=rhoavs+min(psv(Tsurf),pfrost)/Tsurf
            nm=nm+1
         endif

         if (time>=tmax-solsperyear-2.and.mode==0) then
            if (mod(time,1.d0)-oldtime<0.) then
               if (din==100) then  ! 100 steps in a day
                  dirholowz=1.e32
                  do i=1,nz
                     dirho(i)=dirho(i)/100.
                     if (dirho(i)<dirholowz) dirholowz=dirho(i)
                     if (dirho(i)<dirholowLs(i)) dirholowLs(i)=dirho(i)
                  enddo
                  !write(38,'(f7.3,1x,g9.3)') marsLs*180./pi,dirholowz
               endif
               oldtime=mod(time,1.d0)
               do i=1,nz
                  dirho(i)=0.
                  din=0
               enddo
            endif
            do i=1,nz
               dirho(i) = dirho(i)+psv(T(i))/T(i)
            enddo
            din = din+1
         endif

         if (mode==1) then
            if (marsLs<marsLsold) then
               Tmean2=Tmean2/nm
               if (abs(Tmean2-Tbold)<0.05) exit
               Tbold=Tmean2
               Tmean2=0.; nm=0
            else
               Tmean2=Tmean2+T(nz)
               nm=nm+1
            endif
         endif
         marsLsold = marsLs

      enddo 

      if (mode==0) then
         rhoavs=rhoavs/nm
         do i=1,nz
            rhosatav(i)=rhosatav(i)/nm
         enddo
         Tmean1 = Tmean1/nm
         !write(31,'(999(f6.2,1x))') (Tpeak(i),i=1,nz)
         !write(32,'(999(f6.2,1x))') (Tlow(i),i=1,nz)
         !write(35,'(999(g10.4,1x))') (rhosatav(i),i=1,nz)
         if (zdepth<=0.or.zdepth>=z(nz)) then
            avdrho = rhosatav(nz)-rhoavs ! no ice
         else
            do i=1,nz
               if (z(i)>zdepth) then ! ice
                  avdrho = rhosatav(i)-rhoavs
                  exit
               endif
            enddo
         endif
         write(34,'(f7.3,1x,2(f6.2,1x),2(g10.4,1x))') 
     &        Tmean1,Tlow(nz),Tpeak(nz),avdrho+rhoavs,rhoavs
         !write(39,'(999(1x,g9.3))') (dirholowLs(i)-rhoavs,i=1,nz)
      endif
      Tb = T(nz)
      end

