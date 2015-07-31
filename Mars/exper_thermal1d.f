      program exper_thermal1d
C***********************************************************************
C   exper_thermal1d:  program to calculate the diffusion of temperature 
C                     and vapor with prescribed surface temperature
C
C   Grid: surface is at z=0
C         see subroutine vapordiffusion for restrictions in grid spacing
C***********************************************************************

      implicit none
      integer NMAX
      real*8 pi, marsDay, icedensity, rescale
      parameter (NMAX=1000)
      parameter (pi=3.1415926535897932, marsDay=88775.244, rescale=1.)
      parameter (icedensity=926.*rescale) ! must be the same as in vapordiffusion

      integer nz, nsteps, n, i, outstep, nm, b1, b2
      real*8 T(NMAX),tmax, time, dt, zmax, dz, p(NMAX)
      real*8 thIn, dtsec
      real*8 rhoc, delta, rhof(NMAX)
      real*8 Tsurf,Tsurfp1,p0,Diff,D(NMAX),psv,psurf,psurfp1
      real*8 Tmean, Tampl, Told(NMAX), period, truetime
      real*8 z(NMAX), ti(NMAX), rhocv(NMAX), zfac
      real*8 cumice, cumvapor, fill, f(NMAX), f0, cumads, Fsurf
      real*8 rhoav(NMAX), rhoavs, rhoavb1, rhoavb2, zdepth, D0
      real*8 buf,minrho
      character*100 dum1
      real*8 Tsurface, Tsurface_some, rhofold(NMAX), porosity
      external Tsurface,psv,Tsurface_some

C-------read input       
      open(unit=20,file='exper.par',status='old')
      read(20,'(a)')dum1
      read(20,*)dt,tmax,period
      read(20,'(a)')dum1
      read(20,*)nz,zmax,zfac
      read(20,'(a)')dum1
      read(20,*)thIn,rhoc
      read(20,'(a)')dum1
      read(20,*)p0,Tmean,Tampl,Diff,fill
      close(20)

C     set some constants
      dz = zmax/nz    ! set the space step
      nsteps = int(tmax/dt)       ! calculate total number of timesteps
      delta = thIn/rhoc*sqrt(marsDay*669./pi)  ! skin depth

      open(unit=30,file='param.dat',status='unknown')
      write(30,*) 'Model parameters'
      write(30,*) 'Time step=',dt,' Max number of steps=',nsteps
      write(30,*) 'dz=',dz,' zmax=',zmax,' zfac=',zfac
      write(30,*) 'Thermal inertia=',thIn,' rho*c=',rhoc
      write(30,*) 'Mean Temperature=',Tmean, 
     &     ' Temperature amplitude=',Tampl
      write(30,*) 'Diffusivity=',Diff,' initial filling=',fill
      write(30,*) 'Skin depth= ',delta
      write(30,*) 'max. atm. H2O pressure= ',p0
      write(30,*) 'density of solid ice= ',icedensity
      close(30)

      if (nsteps<0) stop 'integer overflow'
      dtsec=dt*period
      outstep=max(5000,nint(669./dt/12.),nsteps/5000)

C     Initialize and write to file
      open(unit=22,file='Tprofile',status='unknown') ! temperature profile
      open(unit=23,file='pprofile',status='unknown') ! H2O pressure profile
      open(unit=24,file='iprofile',status='unknown') ! ice profile
      open(unit=25,file='ideep',status='unknown') ! depth-integrated contents
      open(unit=26,file='aprofile',status='unknown') ! adsorbate profile
      open(unit=27,file='rhoav',status='unknown') ! mean vapor density profile

      do i=1,nz
         T(i) = Tmean  ! initial temperature profile
         p(i) = min(p0,psv(T(i)))  ! initial pressure profile
         rhof(i)= fill*icedensity  ! initial ice profile
         D(i) = Diff
         z(i) = (i-0.5)*dz
         ti(i) = thIn
         rhocv(i) = rhoc
         rhoav(i) = 0.
      enddo
      D0=Diff
      nm=0; rhoavs=0.; 
      rhoavb1=0.; rhoavb2=0.
      if (zfac>1.) then
         dz=zmax/(3.+2.*zfac*(zfac**(nz-2)-1.)/(zfac-1.))
         print *,'allowed explicit stepsize=',(2*dz)**2/(2.*Diff*period)
         z(1)=dz; z(2)=3*z(1)
         do i=3,nz
            z(i)=(1.+zfac)*z(i-1)-zfac*z(i-2)
         enddo
      else
         print *,'using regularly space grid'
         print *,'allowed explicit stepsize=',dz**2/(2.*Diff*period)
      endif

      zdepth=-9999.
      porosity=0.4   ! only matters when ice is present
      b1=nz-1; b2=nz
      do i=1,nz
         if (z(i)<zdepth) then
            rhof(i)=0.
            b1=i; b2=i+1
         endif
      enddo
      if (b1==nz) then; b1=nz-1; b2=nz; endif
      call smartgrid(nz,z,zdepth,thIn,rhoc,porosity,ti,rhocv,1,0.d0)

      open(unit=30,file='z',status='unknown');
      write(30,*) (z(i),i=1,nz)
      close(30)
!      call feedback(nz,porosity,rhof,icedensity,rhoc,thIn,rhocv,ti,
!     &     T,rhof)
      time=0.

!      Tsurf=253.; Tsurfp1=Tsurf; 
!      call assignTprofile(z,nz,253.d0,233.d0,T)
!      write(22,'(f8.4,1000(x,f6.2))') time,(Th+T(i),i=1,nz)
      Tsurf=Tsurface(time,Tmean,Tampl)  ! prescribed surface temperature
!      Tsurf=Tsurface_some(time) 
C-------loop over time steps 
      do n=0,nsteps-1
         truetime =(n+1)*dtsec   !   time at n+1 
         time = (n+1)*dt
         Tsurfp1=Tsurface(truetime,Tmean,Tampl)
!         Tsurfp1=Tsurface_some(truetime) 
         do i=1,nz; Told(i)=T(i); enddo
!         call cranknT(nz,dz,dtsec,T,Tsurf,Tsurfp1,thIn,rhoc)
!        T(nz)=233.
        call conductionT(nz,z,dtsec,T,Tsurf,Tsurfp1,ti,rhocv,0.d0,Fsurf)
         psurf=min(p0,psv(Tsurf))
         do i=1,nz
!            D(i)=Diff*(T(i)/200.)**1.5
            D(i)=Diff
         enddo
!         D0=Diff*(Tsurf/200.)**1.5
         D0=Diff
         do i=1,nz; rhofold(i)=rhof(i); enddo
!         call vapordiffusion(nz,dz,dtsec,D,p,psurf,Told,T,rhof,Tsurf,D0)
         call vapordiffusioni(nz,z,dtsec,D,p,psurf,Told,T,rhof,Tsurf,D0)
         psurfp1=min(p0,psv(Tsurfp1))
!         call diffusionadsorption(nz,z,dtsec,D,p,Told,T,
!     &        psurf,psurfp1,Tsurf,Tsurfp1)
         Tsurf=Tsurfp1
!         call feedback(nz,porosity,rhof,icedensity,rhoc,thIn,rhocv,ti,
!     &        T,rhofold)

c--------only output and diagnostics below this line
         if (time>=0.and.mod(n,outstep/5)==0) then
!         if (time>=tmax-669.and.mod(n,nint(669./dt/512.))==0) then
!            call adsorption2(p,T,f,nz,psurf,Tsurf,f0)
            cumice=0.
            cumvapor=psurf/Tsurf*z(1)/2.
            cumads=f0*z(1)/2.
            minrho=1.e32
            do i=1,nz
               if (i==1) dz=z(2)/2
               if (i>1.and.i<nz) dz=(z(i+1)-z(i-1))/2.
               if (i==nz) dz=(z(nz)-z(nz-1))/2.
               cumice = cumice + rhof(i)*dz
               cumvapor= cumvapor + p(i)/T(i)*(1.-rhof(i)/icedensity)*dz
               cumads = cumads + f(i)*dz
               buf=psv(T(i))/T(i)
               if (buf<minrho) minrho=buf
            enddo
            cumvapor = cumvapor*18./8314.5
            cumice = cumice/icedensity
            write(25,'(f12.0,x,g17.11,2(x,g11.5))') 
     &           time,cumice,cumvapor,cumads
            buf = psurf/Tsurf*18./8314.5
            minrho = minrho*18./8314.5
!            write(25,'(f12.0,x,g17.11,5(x,g11.5))') 
!     &           time,cumice,cumvapor,cumads,buf,minrho
         endif
c        output from last year or last day
         if (time>=tmax-669.and.mod(n,nint(669./dt/512.))==0) then
!         if (time>=tmax-1..and.mod(n,nint(1./dt/32.))==0) then
!            write(22,'(f12.0,999(x,f6.2))') time,Tsurf,(T(i),i=1,nz)
!            write(23,'(f12.0,999(x,f7.4))') time,psurf,(p(i),i=1,nz)
!            write(24,'(f12.0,999(x,g9.3))') time,
!     &           (rhof(i)/icedensity,i=1,nz)
!            call adsorption2(p,T,f,nz,psurf,Tsurf,f0)
!            write(26,'(f12.0,999(x,g11.5))') time,f0,(f(i),i=1,nz)
!            write(27,'(f12.0,999(x,g11.5))') 
!     &           time,psurf/Tsurf,(psv(T(i))/T(i),i=1,nz)
         endif
c        output in equally spaced time intervals
         if (time>=0..and.mod(n,max(1,nsteps/60))==0) then
!         if (time>=0..and.mod(n,max(1,nsteps/897))==0) then
            write(22,'(f12.0,999(x,f6.2))') time,Tsurf,(T(i),i=1,nz)
            write(23,'(f12.0,999(x,f7.4))') time,psurf,(p(i),i=1,nz)
            write(24,'(f12.0,999(x,g10.4))') 
     &           time,(rhof(i)/icedensity,i=1,nz)
!            call adsorption2(p,T,f,nz,psurf,Tsurf,f0)
!            write(26,'(f12.0,999(x,g11.5))') time,f0,(f(i),i=1,nz)
         endif
c        vapor densities in last year
         if (time>=tmax-669.) then
            do i=1,nz
               rhoav(i) = rhoav(i)+p(i)/T(i)
!               rhoav(i) = rhoav(i)+psv(T(i))/T(i) 
            enddo
            nm=nm+1
            rhoavs = rhoavs+min(psv(Tsurf),p0)/Tsurf
            rhoavb1= rhoavb1+psv(T(b1))/T(b1)
            rhoavb2= rhoavb2+psv(T(b2))/T(b2) 
         endif
      enddo                     ! end the loop

      rhoavs=rhoavs/nm; 
      rhoavb1=rhoavb1/nm; rhoavb2=rhoavb2/nm
      do i=1,nz
         write(27,*) z(i),rhoav(i)/nm
      enddo
!      print *,rhoavs,rhoavb1

      close(22)
      close(23)
      close(24)
      close(25)
      close(26)

      end
 



      real*8 function Tsurface(time,Tmean,Tampl)
C     sinusoidal surface temperature
      implicit none
      real*8 time, omegaMars, pi, Tmean, Tampl
      parameter (pi=3.1415926535897932)
      parameter (omegaMars=2.*pi/(669.*88775.244))
      Tsurface=Tmean + Tampl*sin(omegaMars*time)
      return
      end
      



      real*8 function Tsurface_some(time)
C     periodic surface temperature
      implicit none
      real*8 time, omegaMars, pi, omegaMars2
      parameter (pi=3.1415926535897932)
      parameter (omegaMars=2.*pi/(668.6*88775.244))   ! annual
      parameter (omegaMars2=2.*pi/88775.244)   ! diurnal
      real*8 a0,a1,a2,b1,b2  !,a3,b3,a4,b4
!      parameter (a0= 198., a1= 28., b1=-4.,  a2=0.,  b2=5.)  ! VL2
      parameter (a0= 210., a1= -39.65, b1=-1.46,  a2=0.65,  b2=3.5) ! ~merA
      Tsurface_some = a0 +  
     &     a1*sin(omegaMars*time) + b1*cos(omegaMars*time) +
     &     a2*sin(2.*omegaMars*time) + b2*cos(2.*omegaMars*time) 
!     &     + a3*sin(omegaMars2*time) + b3*cos(omegaMars2*time) +
!     &     a4*sin(2.*omegaMars2*time) + b4*cos(2.*omegaMars2*time)
      return
      end



      subroutine assignTprofile(z,nz,Tsurf,Tbottom,T)
C***********************************************************************
C   assignTprofile:  assign time-independent temperature profile
C
C   INPUTS: z = grid point coordinates
C           nz = number of grid points
C           Tsurf = T(z=0)
C           Tbottom 
C   OUTPUT: T = vertical temperature profile [Kelvin]
C***********************************************************************
      implicit none
      integer NMAX
      parameter (NMAX=1000)
      integer nz, i
      real*8 z(NMAX), T(NMAX), Tsurf, Tbottom
!      print *,Tsurf,Tbottom
      do i=1,nz
         T(i) = Tsurf + z(i)/z(nz)*(Tbottom-Tsurf);
      enddo
      end
 




