subroutine oneasteroid(latitude, omega, eps, thIn, Qmean, Tmean, Tmin, Tmax)
  ! calculate surface temperatures of airless body with rotation and conduction
  use constants, only : pi, So, d2r
  use body, only : semia, ecc, Trot, emiss, albedo
  implicit none
  integer n, nr, cc
  real(8), intent(IN) :: latitude, omega, eps, thIn
  real(8), intent(OUT) :: Qmean, Tmean, Tmin, Tmax
  integer, parameter :: nz=110
  integer, parameter :: EQUILTIME=20  ! [orbits]
  real(8) edays, Rau, Ls, decl, HA, Torb, Q, Qp1, dt
  real(8) z(nz), zmax, delta1, delta2
  real(8) coslat, Temp(nz), Tsurf, Fsurf, ti(nz), rhoc(nz)
  !real(8) c, k
  !real(8) Emean, E
  !real(8) rhocmean, Imean
  real(8), parameter :: slope = 0.*d2r ! [radians]
  real(8), parameter :: az = 0.*d2r ! radians east of north
  real(8), external :: flux_noatm, flux2T, a2Torb, sublrate, heatcapacity

  print *,'Latitude=',latitude/d2r
  Torb = a2Torb(semia)
  nr = 50*nint(Torb/Trot)   ! usually 50, more for low thermal inertia <10
  dt = Torb/nr*86400
  print *,'dt=',dt/3600.,'hours','  Trot=',Trot*24.,'hours'
  
  ti(:) = thIn
  rhoc(:) = 2500.*(1-0.4)*400.
  !rhoc(:) = 2500.*(1-0.4)*heatcapacity(120.d0)
  !rhoc(:) = 930*1300. ! solid ice
  
  delta1 = thIn/rhoc(1)*sqrt(Trot*86400/pi)
  delta2 = thIn/rhoc(1)*sqrt(Torb*86400/pi)
  zmax = 5.*delta2
  print *,'skin depth, orbital =',delta2,'diurnal=',delta1
  print *,'nz =',nz
  call setgrid(nz,z,zmax,1.05d0)
  cc=0
  do n=1,nz
     if (z(n)<delta1) cc=cc+1
  enddo
  print *,cc,'points in diurnal skin depth'
  if (cc<5) stop 'Exiting because grid is too coarse'
  !open(unit=30,file='z',action='write',status='unknown')
  !write(30,'(999(f8.4,1x))') z(:)
  !close(30)

  coslat = max(cos(latitude),cos(latitude+eps/2),cos(latitude-eps/2))
  Tsurf=flux2T(So/semia**2*coslat/pi,albedo,emiss) ! effective temperature
  if (abs(coslat)<0.01) Tsurf=50;  ! avoids NaNs
  print *,'initialization temperature=',Tsurf
  Temp(:)=Tsurf
  Tmean=0.; Qmean=0.; 
  Tmin=1e32; Tmax=-1e32
  !Emean=0.
  !rhocmean =0.; Imean=0.; 
  call generalorbit(0.d0,semia,ecc,omega,eps,Ls,decl,Rau) 
  HA = 0.
  Q = (1-albedo)*flux_noatm(Rau,decl,latitude,HA,slope,az)

  do n=1,EQUILTIME*nr
     !if (mod(n,10000)==0) print *,real(n)/real(EQUILTIME*nr)*100,'% done'

     edays = n*Torb/real(nr)
     call generalorbit(edays,semia,ecc,omega,eps,Ls,decl,Rau) 
     HA = mod(edays/Trot,1.d0)*2.*pi  ! Trot in solar days
     Qp1 = (1-albedo)*flux_noatm(Rau,decl,latitude,HA,slope,az)

     !do j=1,nz
     !   c = heatcapacity(Temp(j))
     !   k = conductivity(Temp(j))
     !   ti(j) = sqrt(k*c*500.)
     !   rhoc(j) = 500.*c
     !enddo

     call conductionQ(nz,z,dt,Q,Qp1,Temp,ti,rhoc,emiss,Tsurf,0.d0,Fsurf)
     Q = Qp1
     if (n>(EQUILTIME-1)*nr) then  ! annual average
        !if (n>=EQUILTIME*nr-100) write(71,*) edays,HA,Tsurf,Temp(nz)
        Qmean = Qmean+Q
        Tmean = Tmean+Tsurf
        if (Tsurf>Tmax) Tmax=Tsurf
        if (Tsurf<Tmin) Tmin=Tsurf
        !rhocmean = rhocmean + rhoc(1)
        !Imean = Imean + ti(1)
        !E=sublrate(Tsurf)*18.015*1.66054e-27
        !Emean = Emean+E
     endif
  enddo
  Qmean = Qmean/nr; Tmean = Tmean/nr
  !Imean = Imean/nr; rhocmean = rhocmean/nr
  !Emean = Emean/nr
  !print *,'Thermal inertia=',Imean !,rhocmean/500. 
  !print *,'Sublimation loss=',Emean*86400*365.24,'kg/m^2/yr = mm/yr'
end subroutine oneasteroid



function conductivity(T)
  implicit none
  real(8), intent(IN) :: T
  real(8) k, conductivity
  !real(8) kc, rhoc150
  real(8) A, B

  ! rhoc150 = 500.*433.96

  ! thermal conductivity
  !kc = I150**2/rhoc150/(1+1.48*(150./350.)**3)
  !k = kc*(1+1.48*(T/350.)**3)

  ! based on Sakatani et al. (2012) 
  A = 0.001; B = 3e-11     ! D= 100um
  !A = 0.003; B = 1.5e-10   ! D= 1mm
  k = A + B*T**3

  conductivity = k
end function conductivity
