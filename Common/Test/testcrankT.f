      PROGRAM testcrankT
C***********************************************************************
C     test Crank Nicolson subroutine 
C***********************************************************************
      !use conductionT   ! for test of conductionT2.f90
      implicit none
      integer NMAX
      real*8 pi, Fgeo, Ta, Tm, Period
      parameter (NMAX=1000)
      parameter (pi=3.1415926535897931, Period=88775.244*670)
      parameter (Fgeo = 0., Ta=30., Tm=190.)

      integer nz, n, i, STEPSPERSOL
      real*8 T(NMAX), time, dt
      real*8 thIn  ! thermalInertia
      real*8 rhocv(NMAX)  ! volumetric heat capacity, rho*c
      real*8 Tsurf, Tsurfp1
      real*8 delta, ti(NMAX), Fsurf, z(NMAX)
      real*8 zmax, zfac

      STEPSPERSOL = 120
      dt = Period/STEPSPERSOL
      nz = 70
      zmax = 2.5; zfac=1.02
      thIn = 120.

      rhocv(:) = 1200.*800.
      delta = thIn/rhocv(1)*sqrt(Period/pi)
      print *,'Skin depth= ',delta
      print *,'zmax=',zmax
      print *,'Time step=',dt
      print *,'Thermal inertia=',thIn,' Period=',Period

C     Initialize
      open(unit=22,file='Tprofile',status='unknown') ! temperature profile

      do i = 1,nz
         T(i) = 0.
         ti(i) = thIn
      enddo

      !do i = 1,nz
      !   z(i) = (i-0.5)*dz
      !enddo
      call setgrid(nz,z,zmax,zfac)
      open(unit=30,file='z',status='unknown');
      write(30,*) (z(i),i=1,nz)
      close(30)
      
      time = 0.

      !call conductionT2_init(nz,z,dt,ti,rhocv,Fgeo)
      Tsurf = Tm + Ta*sin(2*pi*time/Period)
      
      do n=0,50000
c         print *,time/Period,Qn
         time = (n+1)*dt          !   time at n+1;
         Tsurfp1 = Tm + Ta*sin(2*pi*time/Period)
         call conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhocv,Fgeo,Fsurf)
         !call conductionT2(nz,Tsurf,Tsurfp1,T,Fsurf)
         Tsurf=Tsurfp1

         !write(21,'(f12.6,2f10.4)') time/Period,Tsurf

         ! write 12 profiles from the last sol
         if (n>50000-STEPSPERSOL) then
            if (mod(n,10)==9) then  ! synched with test_Tprofile.m
               print *,time/Period,Tsurf
               write(22,'(1000(1x,f7.2))') Tsurf,(T(i),i=1,nz)
            endif
         endif

      enddo                     ! end the loop

      close(22)
      end
