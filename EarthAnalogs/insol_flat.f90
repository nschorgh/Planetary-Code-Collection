program insol_flat
!*********************************************************************
!   program to calculate short-wavelength solar irradiance 
!           at flat location on Earth
!*********************************************************************
  use dateformat
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: zero=0.

  integer nsteps, n, nm
  real(8) tmax   ! integration period in minutes
  real(8) latitude, R, dtmin
  type(cTime) udtTime
  real(8) azSun, sinbeta, Qn, Qmean(4), Qmax, S0
  real(8) dZenithAngle, dAzimuth, longitude
  real(8) Qatm, I0, D0, Qdirect, Qdiffuse
  real(8), external :: flux_wshad, mk_atmosphere

  ! set some constants
  dtmin = 5.  ! time step in minutes
  dtmin = 10.
  !tmax = 1.*1440
  tmax = 365.*1440

  ! Mauna Kea summit
  latitude = 19.821; longitude = -155.468

  nsteps=int(tmax/dtmin)       ! total number of timesteps

  ! start time in UTC = HST-10
  !udtTime = cTime(2012,11,26,0.,0.,0.) ! = Nov 26, 14 HST 
  udtTime = cTime(2013,1,1,0.,0.,0.) 

  write(*,*) 'Time step=',dtmin,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude

  Qmean(:) = 0.; Qmax = 0.
  nm=0   

  ! loop over time steps 
  do n=0,nsteps-1

     call sunpos(udtTime, longitude, latitude, dZenithAngle, dAzimuth, R)
     sinbeta = cos(dZenithAngle*d2r)
     azSun = dAzimuth*d2r

     Qn = flux_wshad(R,sinbeta,azSun,zero,zero,zero)  
     S0 = 1365./R**2
     Qatm = S0*mk_atmosphere(dZenithAngle*d2r,I0,D0)
     Qdirect =  Qn*I0
     Qdiffuse = S0*D0
     ! Qatm = Qn*I0 + S0*D0

     write(*,'(i4,1x,2(i2,1x),2(f3.0,1x),4(f6.1,1x),2(f6.4,1x))') &
          & udtTime%iYear,udtTime%iMonth,udtTime%iDay,udtTime%dHours,udtTime%dMinutes, &
          & Qn,Qatm,Qdirect,Qdiffuse,I0,D0

     if (n>nsteps-int(365.*1440./dtmin)) then
        Qmean(1) = Qmean(1) + Qn
        Qmean(2) = Qmean(2) + Qatm
        Qmean(3) = Qmean(3) + Qdirect
        Qmean(4) = Qmean(4) + Qdiffuse
        if (Qmax<Qn) Qmax = Qn
        nm=nm+1
     endif

     udtTime%dMinutes = udtTime%dMinutes + dtmin
     call addtime(udtTime)

  enddo  ! end the loop

  Qmean = Qmean/nm

  print *,'Qmean=',Qmean
end program insol_flat



