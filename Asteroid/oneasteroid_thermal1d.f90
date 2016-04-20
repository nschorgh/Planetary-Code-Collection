subroutine oneasteroid(latitude, omega, eps, albedo, thIn, Qmean, Tmean, Tmin, Tmax)
  ! calculate surface temperatures of body with rotation and conduction
  ! modules are defined in main program
  use constants, only : pi, sigSB, So, d2r
  use body, only : a, ecc, Trot, emiss
  implicit none
  integer n, nr, cc
  real(8), intent(IN) :: latitude, omega, eps, albedo, thIn
  real(8), intent(OUT) :: Qmean, Tmean, Tmin, Tmax
  integer, parameter :: nz=110
  integer, parameter :: EQUILTIME=20  ! [orbits]
  real(8) edays, Rau, Ls, decl, HA, Torb, Q, Qp1, dt
  real(8) z(nz), zmax, delta1, delta2
  real(8) Temp(nz), Tsurf, Fsurf, ti(nz), rhoc(nz)
  !real(8) c, k
  !real(8) Emean, E
  !real(8) rhocmean, Imean
  real(8) slope,az
  real(8), external :: flux_noatm, flux2T, a2Torb, sublrate

  print *,'Latitude=',latitude/d2r
  Torb = a2Torb(a)
  nr = 100*nint(Torb/Trot)   ! usually 50, more for low thermal inertia <<10
  dt = Torb/nr*86400
  print *,'dt=',dt/3600.,'hours','  Trot=',Trot*24.,'hours'
  ti(:) = thIn
  ! for one-layer model the value of rho*c does not matter
  !rhoc(:) = 500.*500.
  rhoc(:) = 930*1300. ! solid ice
  delta1 = thIn/rhoc(1)*sqrt(Trot*86400/pi)  ! Trot has meaning of solar day
  delta2 = thIn/rhoc(1)*sqrt(Torb*86400/pi)
  zmax = 5.*delta2
  print *,'skin depth, orbital =',delta2,'rotat=',delta1
  print *,'nz =',nz
  call setgrid(nz,z,zmax,1.05d0)
  cc=0
  do n=1,nz
     if (z(n)<delta1) cc=cc+1
  enddo
  print *,cc,'points in diurnal skin depth'
  if (cc<5) then
     stop 'Exiting because grid is too coarse'
  endif
  !open(unit=30,file='z',action='write',status='unknown')
  !write(30,'(999(f8.4,1x))') z(:)
  !close(30)

  slope=0.  ! radians
  az = 0.  ! radians east of north
  !slope = 30.*d2r; az = 270.*d2r

  Tsurf=flux2T(So/a**2*(1.-albedo)*cos(latitude)/pi)
  if (abs(cos(latitude))<0.01) Tsurf=50;  ! avoids NaNs
  print *,'initialization temperature=',Tsurf
  Temp(:)=Tsurf
  Tmean=0.; Qmean=0.; 
  Tmin=1e32; Tmax=-1e32
  !Emean=0.
  !rhocmean =0.; Imean=0.; 
  call generalorbit(0.d0,a,ecc,omega,eps,Ls,decl,Rau) 
  HA = 0.
  Q= (1-albedo)*flux_noatm(Rau,decl,latitude,HA,slope,az)

  do n=1,EQUILTIME*nr
     !if (mod(n,10000)==0) print *,real(n)/real(20*nr)*100,'% done'

     edays = n*Torb/real(nr)
     call generalorbit(edays,a,ecc,omega,eps,Ls,decl,Rau) 
     HA = mod(edays/Trot,1.d0)*2.*pi  ! Trot in solar days
     Qp1= (1-albedo)*flux_noatm(Rau,decl,latitude,HA,slope,az)

     !do j=1,nz
     !   call props(Temp(j),thIn,c,k)
     !   ti(j) = sqrt(k*c*500.)
     !   rhoc(j) = 500.*c
     !enddo

     call conductionQ(nz,z,dt,Q,Qp1,Temp,ti,rhoc,emiss,Tsurf,0.d0,Fsurf)
     Q = Qp1
     if (n>(EQUILTIME-1)*nr) then
        !if (n>=EQUILTIME*nr-100) write(71,*) edays,HA,Tsurf,Temp(nz),flux2T(Q)
        Qmean = Qmean+Q
        Tmean = Tmean+Tsurf
        if (Tsurf>Tmax) Tmax=Tsurf
        if (Tsurf<Tmin) Tmin=Tsurf
        !rhocmean = rhocmean + rhoc(1)
        !Imean = Imean + ti(1)
        !E=sublrate(Tsurf)*18.015*1.66054e-27
        !if (mod(n,1007)==0) then  ! not a multiple of 50
        !   write(60,*) edays,Tsurf,E
        !endif
        !Emean = Emean+E
     endif
  enddo
  Qmean = Qmean/nr; Tmean = Tmean/nr
  !Imean = Imean/nr; rhocmean = rhocmean/nr
  !Emean = Emean/nr
  !print *,'Thermal inertia=',Imean !,rhocmean/500. 
  !print *,'Sublimation loss=',Emean*86400*365.24,'kg/m^2/yr = mm/yr'
end subroutine oneasteroid


pure function flux2T(Q)
  ! flux Q must already contain albedo effect
  use constants, only : sigSB
  use body, only : emiss
  implicit none
  real(8), intent(IN) :: Q
  real(8) flux2T

  flux2T = (Q/sigSB/emiss)**0.25
end function flux2T


pure function a2Torb(a)
  ! returns orbital period in Earth days
  use constants, only : pi
  implicit none
  real(8), intent(IN) :: a  ! semimajor axis (AU)
  real(8) a2Torb  

  a2Torb = sqrt(4*pi**2/(6.674e-11*1.989e30)*(a*149.598e9)**3)/86400.
end function a2Torb


subroutine props(T,c,k)
  implicit none
  real(8), intent(IN) :: T
  real(8), intent(OUT) :: c, k
  !real(8) kc
  real(8) A, B

  ! heat capacity from Ledlow et al. (1992), <350K
  !c = 0.1812 + 0.1191*(T/300.-1) + 0.0176*(T/300.-1)**2 + &
  !     0.2721*(T/300.-1)**3 + 0.1869*(T/300.-1)**4
  !c = c*1000*4.184  ! cal/(g K) -> J/(kg K)

  ! heat capacity from Winter & Saari (1969),  20K<T<500K
  c = -0.034*T**0.5 + 0.008*T - 0.0002*T**1.5
  c = c*1000   ! J/(g K) -> J/(kg K)

  ! c150 = 433.96
  ! rho = 500

  ! thermal conductivity
  !kc = I150**2/(500.*433.96)/(1+1.48*(150./350.)**3)
  !k = kc*(1+1.48*(T/350.)**3)

  ! based on Sakatani et al. (2012) 
  A = 0.001; B = 3e-11     ! D= 100um
  !A = 0.003; B = 1.5e-10   ! D= 1mm
  k = A + B*T**3
end subroutine props
