      program mars_mapt
C***********************************************************************
C   mars_mapt: program to calculate the vertical subsurface temperature
C              profile for a list of input parameters (can be the 
C              entire globe)
C***********************************************************************
      implicit none
      real*8 pi, d2r
      parameter (pi=3.1415926535897932, d2r=pi/180.)

      integer nz
      real*8 dt, zfac, zdepth, icefrac, latitude, thIn, albedo
      real*8 fracIR, fracDust, Fgeotherm, lon, Tfrost
      real*8 junk, Tb, pfrost, rhoc, psv, patm, junk2
      external psv

C-----set input parameters
      dt = 0.01
      nz=80; zfac=1.05;
      fracIR=0.04; fracDust=0.02
      !Fgeotherm = 0.028
      Fgeotherm = 0.
      icefrac = 0.4
      patm = 600.

      print *,'RUNNING MARS_MAP-TEMPERATURE'
      write(*,*) 'Global model parameters'
      write(*,*) 'nz=',nz,' zfac=',zfac
      write(*,*) 'fracIR=',fracIR,' fracDust=',fracDust
      write(*,*) 'ice fraction=',icefrac
      write(*,*) 'Fgeotherm=',Fgeotherm

      open(unit=20,file='mapgrid2.dat',status='old') ! the only input

      open(unit=30,file='z',status='unknown');
      !open(unit=31,file='Tpeak',status='unknown')
      !open(unit=32,file='Tlow',status='unknown')
      open(unit=33,file='mapgrid3.dat',status='unknown') ! will be identical to input file
      open(unit=34,file='means',status='unknown')
      open(unit=35,file='rhosatav',status='unknown')
      !open(unit=38,file='seasonrho4',status='unknown')
      !open(unit=39,file='seasonrho3',status='unknown')

      do
c        read output of mars_mapi or another file that includes ice table depths
         read(20,*,end=80) lon,latitude,albedo,thIn,Tfrost,zdepth
         print *,lon,latitude,albedo,thIn,Tfrost,zdepth
         if (albedo==-9999. .or. thIn==-9999.) cycle
         rhoc = 800.*(150 + 100.*sqrt(34.2+0.714*thIn)) ! empirical relation
!        rhoc = 1004640.

         pfrost = psv(Tfrost)
c        negative or very large zdepth indicates the absence of ice
         Tb = -9999.
         call jsub(zdepth, latitude*d2r, albedo, thIn, pfrost,
     &        nz/2, rhoc, fracIR, fracDust, patm, Fgeotherm, 2*dt,
     &        zfac, icefrac, 1, Tb, junk, junk2)
         call jsub(zdepth, latitude*d2r, albedo, thIn, pfrost,
     &        nz, rhoc, fracIR, fracDust, patm, Fgeotherm, dt, 
     &        zfac, icefrac, 0, Tb, junk, junk2)
         write(33,'(f7.2,1x,f7.3,2x,f5.3,2x,f5.1,1x,f7.1,2x,f10.4)') 
     &        lon,latitude,albedo,thIn,Tfrost,zdepth ! merely mirrors input
         ! additional output from subroutines
      enddo

 80   continue
      close(20)
      close(30); close(31); close(32); close(33);
      close(34); close(35); close(38); close(39)
      end


