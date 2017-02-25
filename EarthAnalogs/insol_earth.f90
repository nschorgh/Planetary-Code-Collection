program insol_earth
!***********************************************************************
!   insol_earth: program to calculate direct insolation with shadowing, 
!                Earth orbit, and Mauna Kea atmosphere
!
!   written by Norbert Schorghofer 2010-2016
!***********************************************************************
  use filemanager
  use allinterfaces
  use dateformat
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.

  integer nsteps, n, i, j, nm
  real(8) tmax, edays, dtmin, latitude
  real(8) R, dZenithAngle, dAzimuth, longitude
  real(8) azSun, sinbeta, smax
  real(8), dimension(NSx,NSy) :: h, surfaceSlope, azFac
  real(8), dimension(NSx,NSy) :: Qn, Qsw   ! incoming
  real(8), dimension(NSx,NSy) :: Qmeans, Qmax, daytime
  real(8) Qref, Qrefm, alltime
  real(8) I0, D0, S0, unsd  ! atmosphere
  type(cTime) udtTime
  real(8), parameter :: zero=0.
  logical, parameter :: atmosphere=.true.


  dtmin=10.
  tmax = 365.+1
  latitude = 19.821; longitude = -155.468  ! Mauna Kea summit
  ! azimuth in degrees east of north, 0=north facing

  ! set some constants
  nsteps=int(tmax*1440./dtmin)       ! calculate total number of timesteps

  ! start time in UTC = HST-10
  udtTime = cTime(2013,1,1,0.,0.,0.) ! 14:00 HST 
  !udtTime = cTime(2012,5,20,0.,0.,0.) ! 14:00 HST
  
  write(*,*) 'Starting time',udtTime
  write(*,*) 'Time step=',dtmin,'(min)  Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Atmosphere',atmosphere
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  Qmax=0.; Qmeans=0.; Qrefm=0.
  nm=0
  daytime(:,:)=0.; alltime=0.
  
  print *,'...reading horizons file...'
  call gethorizon(0,0,azSun,smax,.TRUE.)

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     edays = (n+1)*dtmin/1440.
     print *,edays
     udtTime%dMinutes = udtTime%dMinutes + dtmin
     call addtime(udtTime)

     call sunpos(udtTime, longitude, latitude, dZenithAngle, dAzimuth, R)
     sinbeta = cos(dZenithAngle*d2r)
     azSun = dAzimuth*d2r

     do i=2,NSx-1
        do j=2,NSy-1
           call gethorizon(i,j,azSun,smax,.FALSE.)
           Qn(i,j)=flux_wshad(R,sinbeta,azSun,surfaceSlope(i,j),azFac(i,j),smax)
        enddo
     enddo
     Qref=flux_wshad(R,sinbeta,azSun,zero,zero,zero)

     if (atmosphere) then
        unsd = mk_atmosphere(dZenithAngle*d2r,I0,D0)
        S0=1365./R**2  ! must be the same as in flux_wshad
        Qsw(:,:) = Qn(:,:)*I0 + S0*D0  ! do not use skysize
        Qref = Qref*I0 + S0*D0
     else
        Qsw = Qn
     endif
     
     if (edays>tmax-365.) then
     !if (edays>tmax-1.) then
     !if (udtTime%iMonth==6) then
        Qmeans = Qmeans + Qsw
        where (Qsw>Qmax) Qmax=Qsw
        Qrefm = Qrefm + Qref
        nm=nm+1
        
        where (Qn>0.) daytime = daytime + dtmin
        alltime = alltime + dtmin
     endif

  enddo  ! end of time loop

  Qmeans=Qmeans/nm

  print *,'Qref=',Qrefm/nm
  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,2(1x,f6.1),1x,f6.3)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmeans(i,j),Qmax(i,j),daytime(i,j)/alltime*24
        !write(21,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),1x,f5.1)') &  ! instanteneous values
        !     & i,j,h(i,j),surfaceSlope(i,j),Qn(i,j),Qmax(i,j)
     enddo
  enddo
  close(21)

end program insol_earth
