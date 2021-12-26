! Subroutines for lunar thermal ice pump
! written by Norbert Schorghofer, 2013
! equations are in Schorghofer & Aharonson, ApJ 788, 169 (2014)


module oscidea_params
  real(8), parameter :: m = 18*1.66e-27
  real(8), parameter :: theta0 = 1e19  ! #molecules/monolayer
  real(8), parameter :: rough = 1.   ! surface area per grain layer / plane area
end module oscidea_params



subroutine surface_integration(Tm,Ta,supply,spaceweather0, &
     & Tmean,avtheta,fract,Esurf,avweather,mintheta,maxtheta)
  ! solve differential equation for surface population of adsorbed H2O molecules
  use oscidea_params, only : theta0
  use, intrinsic :: ieee_arithmetic
  implicit none
  real(8), intent(IN) :: Tm, Ta, supply, spaceweather0
  real(8), intent(OUT) :: Tmean, Esurf, fract, avweather
  real(8), intent(OUT) :: avtheta, mintheta, maxtheta
  real(8), parameter :: pi = 3.1415926535897932
  real(8), parameter :: lunarDay = 86400.*29.53
  integer, parameter :: nr = 100   ! # of lunar days
  real(8) T, Tn, Tnp1, dtsec, time, totalt
  real(8) loss, spaceweather, dthetadt, theta, theta_np1, df, Ebase, nan
  real(8), external :: sublr_amorph !, deriv

  dtsec = 86400./200.
  if (dtsec>lunarDay/50.) stop 'Time step too large'

  Tmean = 0.
  Esurf = 0.; Ebase = 0.; avweather = 0.
  avtheta = 0.; fract = 0.
  theta = min( theta0*supply/spaceweather0, theta0 )  ! initial guess
  time = 0.
  totalt = 0.
  mintheta = huge(8); maxtheta=-huge(8)
  nan = IEEE_VALUE( real(1.,8), IEEE_QUIET_NAN)
  !Tnp1 = Tm + Ta*sin(-2*pi/lunarDay*time)
  Tnp1 = Tm + Ta*max(sin(-2*pi/lunarDay*time),0.d0)
  
  do  ! time loop 
     Tn = Tnp1
     time = time+dtsec
     !Tnp1 = Tm + Ta*sin(-2*pi/lunarDay*time) 
     Tnp1 = Tm + Ta*max(sin(-2*pi/lunarDay*time),0.d0)   ! average is Tm+Ta/pi
     T = (Tn+Tnp1)/2.

     ! simple Euler step
     !dthetadt = deriv(theta,T,supply,spaceweather0,loss,spaceweather)
     !theta_np1 = theta + dthetadt*dt

     ! stiff ODE
     call funcd(theta,dthetadt,df,T,supply,spaceweather0,loss,spaceweather)
     theta_np1 = theta + dthetadt*dtsec/(1.-df*dtsec) 

     if (theta_np1<0.) theta_np1=0.
     !write(*,*) time/86400.,T,theta/theta0
     if (time>(nr-1)*lunarDay) then
        if (theta<=0.) then
           print *,'UNSTABLE'
           print *,'theta=',theta
           print *,'supply=',supply
           print *,'loss=',loss
           print *,'spaceweather=',spaceweather
           stop
        endif
        Tmean = Tmean + Tn*dtsec
        Esurf = Esurf + loss*dtsec
        Ebase = Ebase + sublr_amorph(Tm)*dtsec  ! for validation purposes only
        avtheta = avtheta + theta*dtsec
        avweather = avweather + spaceweather*dtsec
        if (loss>spaceweather) fract = fract+dtsec
        totalt = totalt + dtsec  ! for validation purposes only
        if (theta>maxtheta) maxtheta=theta
        if (theta<mintheta) mintheta=theta
     endif
     if (time >= nr*lunarDay) exit

     theta = theta_np1
  enddo

  Tmean = Tmean/lunarDay
  avtheta = avtheta/lunarDay
  Esurf = Esurf/lunarDay
  Ebase = Ebase/lunarDay
  avweather = avweather/lunarDay
  fract = fract/lunarDay
  if (ieee_is_nan(avtheta)) maxtheta=nan
  !print *,'Total time',time/lunarDay,'months'
end subroutine surface_integration



function deriv(theta,T,supply,spaceweather0,loss,spaceweather)
  ! physics core
  use oscidea_params, only : theta0, rough
  implicit none
  real(8) deriv
  real(8), intent(IN) :: theta,T,supply,spaceweather0
  real(8), intent(OUT) :: loss, spaceweather
  real(8) v
  real(8), external :: sublr_amorph, padsr

  v = theta/theta0
  loss = sublr_amorph(T)*padsr(v)
  if (v<1e-6) then ! roundoff
     spaceweather = spaceweather0*v/rough
  else
     spaceweather = spaceweather0*(1.-exp(-v/rough))
  endif
  deriv = supply/rough - loss - spaceweather  ! = dthetadt
end function deriv



subroutine funcd(theta,f,df,T,supply,spaceweather0,loss,spaceweather)
  implicit none
  real(8), intent(IN) :: theta, T, supply, spaceweather0
  real(8), intent(OUT) :: f, df, loss, spaceweather
  real(8) th1, th2, f1, f2
  real(8), external :: deriv

  th1 = 0.9*theta; th2=1.1*theta
  f1 = deriv(th1,T,supply,spaceweather0,loss,spaceweather)
  f2 = deriv(th2,T,supply,spaceweather0,loss,spaceweather)
  df = (f2-f1)/(th2-th1)
  f = deriv(theta,T,supply,spaceweather0,loss,spaceweather)
  !print *,'s',theta,f,df
end subroutine funcd
