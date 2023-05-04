PROGRAM testcrankQ_asymp
!***********************************************************************
! test Crank-Nicolson subroutine with nonlinear boundary condition by
! comparing with asymptotic solution for sudden change in incoming flux
! July 2021 -Norbert
!***********************************************************************
  implicit none
  integer, parameter :: nz=50
  real*8, parameter :: pi=3.1415926535897931
  real*8, parameter :: sigSB=5.6704d-8
  real*8, parameter :: emiss=1., Fgeo = 0.

  integer n, i
  real*8 T(nz), time, dt, zmax, zfac
  real*8 thIn, rhocv(nz), Qn, Qnp1, T0, Tanaly
  real*8 ti(nz), Fsurf, Tsurf, z(nz), kappa
  
  dt = 1.
  thIn = 500.
  zmax = 0.2; zfac = 1.05
  rhocv(:) = 1200.*800.
  ti(1:nz) = thIn
  T0 = 200.
  
  print *,'Time step=',dt
  print *,'zmax=',zmax
  print *,'Thermal inertia=',thIn

  open(unit=21,file='Tsurface',status='unknown') ! surface temperature
  
  Tsurf = T0
  T(1:nz) = T0
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z',status='unknown');
  write(30,*) (z(i),i=1,nz)
  close(30)
  kappa = thIn**2/rhocv(1)**2  ! thermal diffusivity
  print *,z(1)**2/(2*kappa)  ! time scale for heat propagation to first grid point
  
  Qn = emiss*sigSB*T0**4
  Qnp1 = Qn

  do n=-10,50-1
     time = (n+1)*dt   !   time at n+1; 
     if (time>0.) Qnp1 = 0. ! turn off the flux
     call conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeo,Fsurf)
     Qn = Qnp1

     ! asymptotic solution from Handelsman & Olmstead (1972)
     if (time>=0.) then
        Tanaly = T0 + 2/sqrt(pi)*emiss*sigSB/thIn*(-T0**4)*sqrt(time)
     else
        Tanaly = T0
     end if
     write(21,'(f12.6,3f9.3)') time,Tsurf,T(nz),Tanaly
  end do    ! end of time loop
  
  close(21)
END PROGRAM testcrankQ_asymp
