program asteroid_fast
  ! Tests the step size increases in asteroid_fast2 
  implicit none
  integer SPINUPN   ! # number of spin-up steps
  real(8) spinupfac 
  parameter(SPINUPN=20, spinupfac=2.)
  integer i, earliest
  real(8) tstart  ! (earth) years
  real(8) icetime, timestep
  real(8) bigstep, bssum, omega

  ! parameters that never change
  tstart = 4.5e9  ! Earth years
  !tstart = 2.0e9
  timestep = 1e5  ! Earth years
 
  omega = 301.

  print *,'RUNNING FAST ASTEROID MODEL'
  print *,'Starting at time',tstart,'years'
  print *,'Time step=',timestep,'years'
  print *,'Spinup',SPINUPN,spinupfac
  print *

  earliest = nint(tstart/timestep)

  icetime = - earliest*timestep
  print *,icetime
  print *,'Spin-up begins here'
  bssum=spinupfac*(spinupfac**SPINUPN-1)/(spinupfac-1.) ! sum_{j=1,n} a^j = a (a^n-1)/(a-1)
  print *,'Spin-up', SPINUPN,'steps over',timestep,'years'
  do i=1,SPINUPN
     bigstep = spinupfac**i/bssum*timestep
     icetime = icetime + bigstep
     print *,i,'of',SPINUPN,'  ',bigstep,omega
     omega = mod(omega + 36.,360.)  ! sweep
  enddo

  
  icetime = -(earliest-1)*timestep
  print *,icetime
  do 
     icetime = icetime + timestep
     print *,icetime
     if (any(-icetime == (/ 4.498d9, 4.450d9, 4d9 /))) then
     !if (any(-icetime == (/ 4.498d9, 4.460d9, 4d9 /))) then
     !if (any(-icetime == 2.0d9 - (/ 0.002d9, 0.05d9, 0.5d9 /))) then
        timestep = 10.*timestep
        print *,'# Increasing time step 10-fold'
     endif
     if (icetime>=0.) exit
     omega = mod(omega + 36.,360.)  ! sweep
  enddo

end program asteroid_fast


