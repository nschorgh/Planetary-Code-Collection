module miscparameters
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: earthDay=86400., marsDay=88775.244, solsperyear=668.60
end module miscparameters


program insol_driver
!***********************************************************************
! Insolation on Mars over the past 20 Myr
!***********************************************************************
  use miscparameters
  implicit none
  integer, parameter :: NP=1  ! number of input latitudes
  integer, parameter :: earliest=5000   ! start time (kyr)
  integer i, k, ierr
  real(8) icetime, ecc, omega, eps, timestep
  real(8), dimension(NP) :: latitude
  real(8), dimension(earliest+1) :: lasktime, laskecc, laskomega, laskeps
  character(10) ext

  !if (iargc() /= 1) stop 'USAGE: a.out ext'
  !call getarg( 1, ext )
  if (command_argument_count() /= 1) stop 'USAGE: a.out ext'
  call get_command_argument( 1, ext )

  open(unit=21,file='lats.'//ext,action='read',status='old',iostat=ierr)
  if (ierr /= 0) then
     print *,'File lats.'//ext,'not found'
     stop
  endif
  do k=1,NP
     read(21,*) latitude(k)
  enddo
  close(21)

  !ecc = 0.0934;  eps = 25.19*d2r;  omega = 250.87*d2r
  ! Laskar orbital solution http://vo.imcce.fr/insola/earth/online/mars/La2003-04/
  open(20,file='INSOLN.LA2004.MARS.ASC',action='read',status='old')
  do i = 1,earliest+1
     read(20,*) lasktime(i),laskecc(i),laskeps(i),laskomega(i)
  enddo
  close(20)
  timestep = abs(lasktime(1)-lasktime(2))*1000.

  print *,'RUNNING MARS ICE AGES - INSOLATION ONLY'
  print *,'Starting at time',lasktime(earliest)*1000.,'years'
  print *,'Time step=',timestep,'years'
  print *,'Number of sites=',NP
  do k=1,NP
     print *,'  Latitude (deg)',latitude(k)
  enddo

  open(unit=34,file='insol.'//ext,action='write',status='unknown')

  do i = earliest,1,-1
     ecc = laskecc(i)
     omega = laskomega(i) + pi
     omega = mod(omega,2*pi)
     eps = laskeps(i)
     icetime = lasktime(i)*1000.
     call icesheet(icetime,NP,latitude,ecc,omega,eps)
  enddo

  close(34)
end program insol_driver



subroutine icesheet(icetime,NP,latitude,ecc,omega,eps)
  ! pass-through subroutine for similarity to big iceage code
  use miscparameters
  implicit none
  integer, intent(IN) :: NP
  real(8), intent(IN) :: icetime, latitude(NP), ecc, omega, eps
  integer k
  real(8) Qmean, Q4, Qpeak
  do k=1,NP 
     call insolonly_mars(latitude(k)*d2r,ecc,omega,eps,Qmean,Q4,Qpeak)
     ! the following format list has been carefully chosen
     write(34,'(f10.0,1x,f6.2,1x,f6.3,1x,f7.5,1x,f5.1,1x,f6.2,1x,f6.4,2x,f5.1)') &
          & icetime,latitude(k),eps/d2r,ecc,omega/d2r,Qmean,Q4,Qpeak
  enddo  
end subroutine icesheet



subroutine insolonly_mars(latitude, ecc, omega, eps, Qmean, Q4, Qpeak)
  use miscparameters
  implicit none
  real(8), intent(IN) :: latitude  ! [radians]
  real(8), intent(IN) :: ecc
  real(8), intent(IN) :: omega, eps  ! [radians]
  real(8), intent(out) :: Qmean, Q4, Qpeak
  real(8), parameter :: dt = 0.01   ! [sol]
  real(8), parameter :: a = 1.52366 ! Mars semimajor axis [AU]
  integer nsteps, n
  real(8) time, Qn, tdays, marsR, marsLs, marsDec, HA
  real(8), external :: flux_noatm

  nsteps=nint(nint(solsperyear)/dt)   ! total number of timesteps
  Qmean=0.; Q4=0.
  Qpeak = -9.

  !-----loop over time steps 
  do n=0,nsteps-1
     time = n*dt
     tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
     call generalorbit(tdays,a,ecc,omega,eps,marsLs,marsDec,marsR)
     HA = 2.*pi*mod(time,1.d0)  ! hour angle
     Qn = flux_noatm(marsR,marsDec,latitude,HA,0.d0,0.d0)
     if (Qn>0.) then ! in case Qn=-0.
        Qmean = Qmean + Qn
        Q4 = Q4 + Qn**0.25
     endif
     if (Qn>Qpeak) Qpeak = Qn
  enddo  ! end of time loop
  
  Qmean = Qmean/nsteps;  Q4 = Q4/nsteps
end subroutine insolonly_mars

