! Monte Carlo model of ballistically hopping molecules
! each particle has 
! x (longitude), y (latitude), status (on surface or in flight), 
! time (until it arrives on surface or until it will leave the surface)

module exo_species
  real(8), parameter :: mmass = 18.015   ! H2O
  !real(8), parameter :: mmass = 19.021   ! HDO
  !real(8), parameter :: mmass = 2.0159   ! H2
  !real(8), parameter :: mmass = 4.0026   ! He-4
  !real(8), parameter :: mmass = 16.0425  ! CH4
  !real(8), parameter :: mmass = 17.007   ! OH
  !real(8), parameter :: mmass = 17.03    ! NH3
  !real(8), parameter :: mmass = 20.18    ! Ne
  !real(8), parameter :: mmass = 22.990   ! Na
  !real(8), parameter :: mmass = 28.01    ! CO
  !real(8), parameter :: mmass = 28.013   ! N2
  !real(8), parameter :: mmass = 35.97    ! Ar-36
  !real(8), parameter :: mmass = 39.962   ! Ar-40
  !real(8), parameter :: mmass = 39.098   ! K
  !real(8), parameter :: mmass = 44.01    ! CO2
  !real(8), parameter :: mmass = 64.06    ! SO2
  !real(8), parameter :: mmass = 131.3    ! Xe
  !real(8), parameter :: mmass = 222.     ! Rn

  ! photodissociation time scale at 1 AU
  !real(8), parameter :: taudissoc = 20.*3600.  ! Potter & delDuca (1964)
  real(8), parameter :: taudissoc = 1/12.6e-6  ! Crovisier (1989)
  !real(8), parameter :: taudissoc = 1/23.0e-6  ! Crovisier (1989), active sun
  !real(8), parameter :: taudissoc = 1.9e7  ! He, Killen & Ip (1999)
  !real(8), parameter :: taudissoc = 3.2e6  ! Ar, Killen & Ip (1999)
  ! Huebner et al. (1992), quite sun
  !real(8), parameter :: taudissoc = 1/1.5e-7 ! H2
  !real(8), parameter :: taudissoc = 1/5.2e-8 ! He
  !real(8), parameter :: taudissoc = 1/7.6e-6 ! CH4
  !real(8), parameter :: taudissoc = 1/7.5e-6 ! OH
  !real(8), parameter :: taudissoc = 1/1.8e-4 ! NH3
  !real(8), parameter :: taudissoc = 1/5.92e-6 ! Na
  !real(8), parameter :: taudissoc = 1/1.0e-6 ! N2
  !real(8), parameter :: taudissoc = 1/3.1e-7 ! Ar
  !real(8), parameter :: taudissoc = 1/2.2e-5 ! K
  !real(8), parameter :: taudissoc = 1/2.0e-6 ! CO2
  !real(8), parameter :: taudissoc = 1/2.1e-4 ! SO2
  !real(8), parameter :: taudissoc = 1/1.5e-6 ! Xe

  ! this module is only used here
end module exo_species




subroutine hop1(p_r, p_s, p_t, idum, Tsurf, Q)
  ! ballistic flight of one particle
  use body, only: g, Rbody, semia, vescape, siderealDay
  use exo_species 
  implicit none
  real(8), intent(INOUT) :: p_r(2) ! (1)=longitude  (2)=latitude
  integer, intent(OUT) :: p_s  ! status
  real(8), intent(INOUT) :: p_t  ! time
  integer, intent(INOUT) :: idum
  real(8), intent(IN) :: Tsurf, Q

  real(8) d,v(3),lat,buf,az,cosaz,sinph2,cosph2,cosdlon,dlon
  real(8) u,flighttime,destr_rate,vspeed,alpha
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer, external :: inbox
  real(8), external :: ran2
  real(8), external :: gasdev  ! gaussian with unit variance

  ! Gaussian/Maxwellian launch velocities
  buf = sqrt(Tsurf*8314.5/mmass)
  ! gasdev has unit variance
  v(1) = gasdev(idum)*buf
  v(2) = gasdev(idum)*buf
  v(3) = abs(gasdev(idum))*buf

  ! v(1) = v(1) - 2*Rbody*pi/siderealDay*cos(p_r(2)*d2r)  ! Coriolis
  vspeed = sqrt(sum(v(:)**2))
  if (vspeed>vescape) then  ! grav. escape
     p_s = -2  ! destroyed by escape
     p_t = 1d100  ! never use again
     return
  endif

  ! molecule moves on plane that goes through center of sphere/body
  ! ground track is thus part of a great circle
  flighttime = 2*v(3)/g   ! time of flight for constant g
  d = 2/g*v(3)*sqrt(v(1)**2+v(2)**2)  ! distance for constant g
  if (vspeed>0.4*vescape) then  ! use non-uniform gravity formula
  !if (d>0.1*Rbody) then
     !gamma = (vspeed/vescape)**2
     ! theta = zenith angle of launch velocity
     !sin2theta = 2*v(3)*sqrt(v(1)**2+v(2)**2)/vspeed**2
     ! derived from an equation in Vogel (1966)
     !d = 2*Rbody*atan(sin2theta/(1/gamma-sin2theta**2))
     !flighttime = flighttime*(1+4./3.*gamma)  ! 1st order correction
     alpha = atan(sqrt(v(1)**2+v(2)**2)/v(3))  ! angle from zenith
     call nonuniformgravity(vspeed,alpha,d,flighttime)
  endif
  !write(70,*) 'flighttime',flighttime,d  ! for flighttime statistics
  az = atan2(v(2),v(1))
  !write(70,*) d,az   ! for statistical tests
  cosaz = v(2)/sqrt(v(1)**2+v(2)**2)
  lat = p_r(2)*d2r
  sinph2 = sin(d/Rbody)*cos(lat)*cosaz + sin(lat)*cos(d/Rbody)
  !write(70,*) asin(sinph2)/d2r-p_r(2)   
  p_r(2) = asin(sinph2)/d2r
  cosph2 = sqrt(1.-sinph2**2)
  if (cosph2/=0) then  ! not on pole
     cosdlon= (cos(d/Rbody)*cos(lat)-sin(lat)*sin(d/Rbody)*cosaz)/cosph2
     if (cosdlon>+1.) cosdlon=+1.  ! roundoff
     if (cosdlon<-1.) cosdlon=-1.  ! roundoff
     dlon = acos(cosdlon)
     if (v(1)<0.) dlon=-dlon 
     !write(70,*) dlon  
     p_r(1) = p_r(1) + dlon/d2r
  else    ! on pole
     ! longitude does not matter 
     p_r(1) = 0.  ! just in case
  endif
  ! p_r(1) = p_r(1) + flighttime/siderealDay*360.  ! Coriolis
    
  if (p_r(2)>90. .or. p_r(2)<-90) then
     print *,'hop1: this cannot happen',p_r(2)
     !stop
  endif
  p_r(1)=modulo(p_r(1),360.)   ! 0 <= p_r(1) < 360.
  
  p_s = 1
  p_t = p_t + flighttime

  ! in-flight destruction
  if (Q>0. .and. taudissoc>0.) then  ! dayside
     u=ran2(idum)  ! between 0 and 1
     !if (u<0.004) then ! 0.4%
     destr_rate = flighttime/(taudissoc*semia**2) ! photodissociation
     if (u < destr_rate) then
        p_s = -1  ! destroyed by photodissociation
        p_t = 1d100  ! never use again
     endif
  endif
end subroutine hop1


function residence_time(T)
  implicit none
  real(8), intent(IN) :: T
  real(8) residence_time
  real(8), parameter :: sigma0 = 1e19  ! H2O
  real(8), external :: sublrate
  residence_time = sigma0/sublrate(T) 
  residence_time = residence_time*400   ! 0.1 monolayers (see S&A, 2014)
  if (T==0.) residence_time = 1e32
  !residence_time = 0. ! Ar, He, (noncondensible species)
end function residence_time


function residence_time2(T,sigma)
  ! residence time that is density dependent
  implicit none
  real(8), intent(IN) :: T, sigma
  real(8) residence_time2
  real(8), parameter :: sigma0 = 1e19
  real(8), external :: sublrate
  real(8) frac
  frac = min(sigma/sigma0,1.)
  residence_time2 = sigma0/sublrate(T)
  residence_time2 = residence_time2/frac
  if (T==0.) residence_time2 = 1e32
end function residence_time2


function residence_timeR(T)
  ! residence time with probability distribution
  implicit none
  real(8), intent(IN) :: T
  real(8) residence_timeR
  real(8), parameter :: sigma0 = 1e19  ! H2O
  real(8), external :: sublrate, ran2
  real(8) tau, y
  integer, save :: idum=-4578
  
  tau = sigma0/sublrate(T)
  y  = ran2(idum)
  residence_timeR = -tau/log(y)  ! gives P(t)=tau/t^2 e^(-tau/t),  <1/t>=1/tau
  if (T==0.) residence_timeR = 1e32
  !residence_time = 0. ! Ar, He, (noncondensible species)
  ! Integrate[tau/t^2 Exp[-tau/t],{t,0,Infinity},Assumptions->tau>0] = 1
  ! Integrate[tau/t^3 Exp[-tau/t],{t,0,Infinity},Assumptions->tau>0] = 1/tau
end function residence_timeR


subroutine montecarlo(Np,idum,p_r,p_s,p_t,p_n,Tsurf,dtsec,ccc,Q) 
  ! called once an hour
  implicit none
  integer, intent(IN) :: np
  real(8), intent(IN) :: Tsurf(*), dtsec, Q(*)
  integer, intent(INOUT) :: idum, p_s(np), p_n(np), ccc(4)
  real(8), intent(INOUT) :: p_r(np,2), p_t(np)
  integer i, k
  real(8) residencetime
  logical, parameter :: VERBOSE = .false.
  logical, external :: incoldtrap
  integer, external :: inbox, insidecoldtrap
  real(8), external :: residence_time

  do i=1,np
     if (VERBOSE) print *,'montecarlo: working on particle',i

     do ! do this for an hour
        if (p_t(i) > dtsec .or. p_s(i)<0) exit

        if (VERBOSE) print *,'next event is',i,p_t(i),p_s(i)
        select case (p_s(i))
        case(-1000:-1)
           exit
        case(0)  ! leaving
           k = inbox(p_r(i,:))
           !write(70,*) i,p_r(i,:), p_s(i), p_t(i), Tsurf(k), Q(k)
           call hop1(p_r(i,:),p_s(i),p_t(i),idum,Tsurf(k),Q(k)) 
           if (p_s(i)==-1) ccc(1)=ccc(1)+1
           if (p_s(i)==-2) ccc(2)=ccc(2)+1
           p_n(i) = p_n(i)+1
        case(1) ! landing
           k = inbox(p_r(i,:))
           if (incoldtrap(p_r(i,:))) then
              !p_s(i)=-100-insidecoldtrap(p_r(i,:))
              if (p_r(i,2)>0.) p_s(i)=-3
              if (p_r(i,2)<0.) p_s(i)=-4
              if (p_s(i)==-3) ccc(3)=ccc(3)+1
              if (p_s(i)==-4) ccc(4)=ccc(4)+1 
              p_t(i)=residence_time(100.d0)
              cycle
           endif
           residencetime = residence_time(Tsurf(k))
           !print *,i,'restime=',k,Tsurf(k),residencetime
           ! reside on surface
           if (VERBOSE) print *,'landed, adding residence time',residencetime,Tsurf(k)
           p_t(i) = p_t(i) + residencetime
           p_s(i)=0
        end select
     enddo  ! end loop over hour

  enddo  ! end loop over particles

  ! subtract dtsec from all times
  where(p_s>=0) p_t = p_t-dtsec
  
end subroutine montecarlo


subroutine production(Np,p_r,p_s,p_n,idum,Tsurf,newcc)
  ! continuous production
  implicit none
  integer, intent(IN) :: Np
  real(8), intent(INOUT) :: p_r(Np,2)
  integer, intent(INOUT) :: p_s(Np), p_n(Np), idum
  real(8), intent(IN) :: Tsurf(*)
  integer, intent(OUT) :: newcc
  integer i, k
  integer, parameter :: NPROD = 2000
  integer, external :: inbox
  real(8), external :: ran2

  newcc=0
  do i=1,Np
     if (p_s(i)<0) then
        p_s(i)=0
        p_n(i)=0

        ! in subsolar region
        do 
           p_r(i,1) = 360.*ran2(idum)
           p_r(i,2) = 20*(2*ran2(idum)-1.)
           k = inbox(p_r(i,:))
           if (Tsurf(k)>360.) exit
        enddo

        ! global
        !p_r(i,1) = 360.*ran2(idum)
        !p_r(i,2) = 90.*(2*ran2(idum)-1.)
        !k = inbox(p_r(i,:))

        ! in equatorial region
        !p_r(i,1) = 360.*ran2(idum)
        !p_r(i,2) = 40*(2*ran2(idum)-1.)

        newcc = newcc+1
     endif
     if (newcc==NPROD) exit
  enddo
  print *,'Created',newcc,'new molecules'
end subroutine production


subroutine destruction(Np,p_r,p_s,p_t,idum,dtsec,veclen,sigma)
  ! destruction on the surface by space weathering
  implicit none
  integer, intent(IN) :: Np, veclen
  real(8), intent(IN) :: p_r(Np,2)
  integer, intent(INOUT) :: p_s(Np)
  real(8), intent(INOUT) :: p_t(Np)
  integer, intent(INOUT) :: idum
  real(8), intent(IN) :: dtsec, sigma(veclen)

  integer k, i
  real(8) u, drate
  real(8), parameter :: sigma0 = 1e19
  real(8) Diso, Ds  ! destruction rate, molecules/m^2/s
  integer, external :: inbox
  real(8), external :: ran2

  Diso = 1e11 ! isotropic destruction rate  1 t/m^2/Ga = 1e12 molec/m^2/s
  Ds = 1e8*1e4*0. ! times cos(incangle)

  do i=1,Np
     if (p_s(i)==0) then
        u = ran2(idum)
        k = inbox(p_r(i,:))
        if (sigma(k)<=sigma0) then
          drate = Diso/sigma0 
        else
          drate = Diso*sigma(k)/sigma0**2
        endif
        if (u<drate*dtsec) then
           p_s(i) = -1  ! destroyed
           p_t(i) = 1d101
        endif
     endif
  enddo
end subroutine destruction


subroutine writeparticles(unit,Np,p_r,p_s,p_t,p_n)
  implicit none
  integer, intent(IN) :: unit, Np
  integer, intent(IN) :: p_s(Np), p_n(Np)
  real(8), intent(IN) :: p_r(Np,2), p_t(Np)
  integer i,k
  integer, external :: inbox
  do i=1,Np
     k = inbox(p_r(i,:))
     write(unit,*) i,p_r(i,:),p_s(i),p_t(i),k,p_n(i)
  enddo
end subroutine writeparticles


subroutine totalnrs(Np,p_s,cc)
  implicit none
  integer, intent(IN) :: Np, p_s(Np)
  integer, intent(OUT) :: cc(6)

  cc(1) = count(p_s==0)   ! on surface
  cc(2) = count(p_s==1)   ! inflight
  cc(3) = count(p_s==-1)  ! destroyed, photo
  cc(4) = count(p_s==-2)  ! destroyed, escape
  cc(5) = count(p_s==-3)  ! coldtrapped, north
  cc(6) = count(p_s==-4)  ! coldtrapped, south
end subroutine totalnrs


logical function incoldtrap(p_r)
  implicit none
  real(8), intent(IN) :: p_r(2)
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  !real(8) u,f
  !integer, external :: insidecoldtrap, ran2
  !integer, save :: idum = -999

  incoldtrap = .FALSE.

  ! approx. relative area of spherical cap (a*pi/180)**2/2, a=sqrt(2*F)*180/pi
  ! approx. relative area of spherical cap 1-cos(a*pi/180), a=acos(1-F)*180/pi

  ! MOON
  ! Mazarico et al. (2011)
  if (p_r(2)> +90-2.11) incoldtrap = .TRUE.  ! 12866 km^2, 0.068%
  if (p_r(2)< -90+2.36) incoldtrap = .TRUE.  ! 16055 km^2, 0.085%
  !if (abs(p_r(2))> 90-2.24) incoldtrap = .TRUE.  ! 14460 km^2, 0.076%
  ! dlat = 0.076e-2/cos(85.*d2r)/d2r
  !if (abs(p_r(2))>85.-dlat/2 .and. abs(p_r(2))<85.+dlat/2.) incoldtrap = .TRUE.

  ! MERCURY
  ! 28,000 km^2 PSR south of 85S, Chabot et al. (2012) = 0.075% of the hemisphere
  !if (abs(p_r(2))> 90-2.2) incoldtrap = .TRUE. 

  ! CERES
  !if (p_r(2)> +90-2.92) incoldtrap = .TRUE.  ! 0.13% of hemisphere
  !if (p_r(2)< -90+2.92) incoldtrap = .TRUE.  ! 0.13% of hemisphere
  ! dlat = 0.13e-2/cos(80.*d2r)/d2r
  !if (abs(p_r(2))>80.-dlat/2 .and. abs(p_r(2))<80.+dlat/2.) incoldtrap = .TRUE.

  !if (abs(p_r(2)) > 85.) then
  !   ! this only works for f<0.38%
  !   u=ran2(idum)  ! between 0 and 1
  !   f = 0.1/100/0.003805302  ! (1-cos(5*d2r))
  !   if (u<f) incoldtrap = .TRUE.
  !endif

  !if (insidecoldtrap(p_r)>0) incoldtrap = .TRUE.
end function incoldtrap


subroutine nonuniformgravity(vspeed,alpha,d,t)
  ! ballistic travel distance and flighttime for non-uniform gravity
  ! not suitable for small velocities due to roundoff
  use body
  implicit none
  real(8), intent(IN) :: vspeed
  real(8), intent(IN) :: alpha  ! zenith angle of launch velocity
  real(8), intent(OUT) :: d, t  ! distance and flighttime
  real(8), parameter :: pi=3.1415926535897932
  real(8) gamma, a, ecc, Ep

  gamma = (vspeed/vescape)**2
  a = Rbody/2./(1-gamma)
  ecc = sqrt(1-4*(1-gamma)*gamma*sin(alpha)**2)
  d = 2*Rbody*acos(1/ecc*(1-2*gamma*sin(alpha)**2))
  Ep = pi - 2*atan(sqrt((1-ecc)/(1+ecc))/tan(d/(4*Rbody)))
  if (ecc>1.-1d-5) then
     d = Rbody*4*gamma*sin(alpha)
     Ep = pi - 2*atan(sqrt((1-gamma)/gamma))
  endif
  t = 2*sqrt(2*a**3/Rbody/vescape**2)*(Ep+ecc*sin(Ep))
  if (1-2*gamma*sin(alpha)**2 > ecc) then ! otherwise d=NaN
     d = 0.
     t = 0.
  endif
end subroutine nonuniformgravity
