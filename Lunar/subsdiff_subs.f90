! subroutines and functions for random walk model of subsurface diffusion 
!    in lunar regolith


real(8) function residence_time2(T,vvm,alphaT,idum)
  ! molecular residence time for desorption or sublimation
  implicit none
  real(8), intent(IN) :: T  ! temperature [K]
  real(8), intent(IN) :: vvm  ! v/vm, number of monolayers
  real(8), intent(IN) :: alphaT  ! fraction of inelastic collisions
  integer, intent(INOUT) :: idum
  real(8), parameter :: pi=3.1415926535, kB=1.38065e-23, m=18*1.66e-27
  !real(8), parameter :: theta1 = (930./m)**(2./3.) ! ~1e19 molecules/m^2
  real(8), parameter :: theta1 = 1e19  ! [molecules/m^2]
  real(8) p, sublrate, tau, y
  real(8), external :: psv, alpha, ran2, vfv_BET, p_BET

  if (T<=0.) then
     residence_time2 = 1.e32
     return
  end if
  
  p = psv(T)
  !alphaT = alpha(T)  ! fraction of inelastic bounces
  sublrate = alphaT * p / sqrt(2*pi*m*kB*T) ! sublimation rate [molecules/m^2/s]
  !desorprate = sublrate * padsr(vvm)  ! desorption rate
  !desorprate = sublrate * p_BET(vvm)

  ! single residence time
  ! theta1*vvm = theta
  if (vvm<=1.) then
     tau = theta1 / sublrate * vfv_BET(vvm)  ! = theta / desorprate
  else ! if mobile
     tau = theta1 / ( sublrate * p_BET(vvm) )  ! = theta1 / desorprate
  endif
  !residence_time2 = tau  ! without probability distribution

  ! residence time selected from probability distribution
  y = ran2(idum)
  residence_time2 = -tau/log(y)
  ! gives P(t) = tau/t^2 e^(-tau/t), <1/t>=1/tau, P(t)d(1/t) ~ e^(-tau/t)

  !if (mod(idum,271)==0) then
  !   print *,vvm, tau, residence_time2
  !endif

end function residence_time2



subroutine hop_s(zidx,idum)
  ! diffuse up or down, random walk
  implicit none
  integer, intent(INOUT) :: zidx, idum
  real(8) dir
  real(8), external :: ran2
  
  if (zidx<0) return
  dir = ran2(idum)
  if (dir>0.5) then  
     zidx = zidx+1  ! down
  else
     zidx = zidx-1  ! up
  endif
end subroutine hop_s



subroutine randomwalk1(zidx, dt, idum, MAXDEPTHIDX, T, vvm)
  ! subsurface random walk of one molecule
  implicit none
  integer, intent(INOUT) :: zidx, idum
  real(8), intent(IN) :: dt
  integer, intent(IN) :: MAXDEPTHIDX
  real(8), intent(IN), dimension(0:MAXDEPTHIDX)  :: T, vvm
  integer, parameter :: GONE = -9999
  real(8) tau, flighttime, alphaT, telapsed
  real(8), external :: alpha, ran2, residence_time2
  
  if (zidx<0 .or. zidx>MAXDEPTHIDX) return
  telapsed = 0.

  do while (telapsed<dt)  ! repeat until dt time has elapsed
  
     ! elastic or inelastic
     !alphaT = alpha(T(zidx))
     alphaT = 1.  ! always inelastic
     !if (ran2(idum)>alphaT) then  ! elastic collision
     !   tau = 0.
     !else ! inelastic, thermally accommodated
        tau = residence_time2( T(zidx), vvm(zidx), alphaT, idum)
     !endif

     !if (VERBOSE) write(27,*) zidx,T(zidx),tau,vvm(zidx)

     telapsed = telapsed + tau
     
     if (telapsed > dt) return ! not moving during this time interval
     call hop_s(zidx, idum)
     if (zidx<0) then
        zidx = GONE  ! lost to space
        return
     end if
     if (zidx>MAXDEPTHIDX) then
        !zidx = GONE   ! open lower boundary
        zidx = MAXDEPTHIDX   ! impermeable lower boundary
        return
     end if
     
     ! add time-of-flight
     flighttime = 0.
     ! vth = sqrt(8*kB*T/pi/m) ! mean thermal speed, 376 m/s at 120K
     !flighttime = dz/400.
     telapsed = telapsed + flighttime

  end do
  ! Note: transition at dt is not seamless
end subroutine randomwalk1
