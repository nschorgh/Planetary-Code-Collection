program sphere1d_implicit
  ! 1D spherically symmetric heat equation and ice retreat
  ! solved with semi-implicit method
  implicit none
  integer, parameter :: nz=100
  integer, parameter :: NB=58 ! max number of bodies
  integer i, ierr, id, unit, narg

  real*8 semia, ecc, Q, time0
  real*8, parameter :: So=1365
  real*8, parameter :: A=0.05, emiss=0.96
  character(2) :: nn
  
  logical :: init(NB) = .true.
  real*8 z(nz), dr, kappa, Radius, dtsec, time(NB), oldtime(NB), Deltat
  real*8 T(nz,NB), Tsurf(NB), Tsurfm1(NB), Teff

  real*8 zT(NB)  ! depth of ice table
  real*8 Tatz, Elatent
  real*8, external :: interp1, flux2T

  zT(:) = 0.
  time(:) = 0.  ! earliest time in input file
  
  Radius = 500.  ! body radius
  dr = Radius/nz
  forall(i=1:nz) z(i) = i*dr

  kappa = 0.25/(0.5*2600*500)  ! kappa=k/rhoc

  narg = COMMAND_ARGUMENT_COUNT()
  if (narg==0) then
     nn='21'
  else
     call getarg(1,nn)
  endif
  print *,'processing extension ',nn
  if (nn(2:2)==' ') then ! single digit input
     nn(2:2)=nn(1:1)
     nn(1:1)='0'
  endif
  
  open(unit=29,file='info_sphere',action='write')
  write(29,*) 'albedo=',A,'emissivity=',emiss
  write(29,*) 'Radius=',Radius,'spatial resolution=',dr
  write(29,*) 'thermal diffusivity=',kappa
  write(29,*) 'initial ice depth=',zT(1)
  close(29)
  
  !open(unit=20,file='/arsia/Orbits/OMBr10201_follow.out',action='read',iostat=ierr,status='old') ! 11 columns
  !open(unit=20,file='/arsia/Orbits/OMBr10201_clean.out',action='read',iostat=ierr,status='old') ! 4 columns
  open(unit=20,file='/arsia/Orbits/OMBr102'//nn//'_clean.out',action='read',iostat=ierr,status='old') ! 4 columns
  if (ierr>0) stop 'input file not found'

  if (nz>1000) stop 'tridag is only set up for N<=1000'

  do  ! time loop
     !read(20,*,iostat=ierr) time0,semia,ecc,u,u,u,u,u,u,u,id
     read(20,*,iostat=ierr) time0,semia,ecc,id
     if (ierr<0) then
        print *,'reached end of file'
        exit
     endif
     if (ierr>0) then
        print *,'read error'
        exit
     endif

     !if (id/=26) cycle  ! follow only one body
     
     oldtime(id) = time(id)
     time(id) = time0 ! years
     Deltat = time(id)-oldtime(id)  ! usually 300 yrs
     !print *,Deltat
     if (Deltat<0. .or. Deltat>1000.) error stop 'Deltat'
     dtsec=Deltat*86400*365.24
     
     !if (time(id)<0.) cycle
     if (id>NB .or. id<1) error stop 'particle id out of range'
     
     !--orbital elements -> surface temperature
     Q = So/semia**2/sqrt(1-ecc**2) ! annual mean insolation
     !Q = Q*faintsun(1e8-time)
     if (.not.init(id)) Tsurfm1(id) = Tsurf(id)
     Teff = flux2T(Q/4.,A,emiss)

     ! options for surface temperature
     Tsurf(id) = Teff ! upper bound
     !call insolonly1(latitude,semia,omega,ecc,obliq,Q0mean,Qmean,Q4mean)
     !print *,Tsurf(id),flux2T(Qmean,A,emiss),flux2T(Q4mean,A,emiss)
     !Tsurf(id) = flux2T(Q4mean,A,emiss)  ! cold end-member
     
     if (init(id)) then ! initialization
        print *,'initializing particle',id,Tsurf(id)
        T(:,id) = Tsurf(id)
        Tsurfm1(id) = Tsurf(id)
        init(id) = .false.
     end if
     
     !--thermal evolution
     call conductionT_sphere(nz,dr,dtsec,T(:,id),Tsurfm1(id),Tsurf(id),kappa)
     
     !--ice evolution
     if (zT(id)>=0. .and. zT(id)<Radius) then
        Tatz = interp1(real(0.,8),z,Tsurf(id),T(:,id),zT(id),nz)
        call retreat_s(Tatz,zT(id),Radius,dtsec,Elatent)
        if (zT(id)>Radius) zT(id)=-9999.
     else
        Tatz = -9999.
     end if

     ! lots of output
     unit = 100+id
     if (id==26) then
        !write(unit,'(f11.0,1x,f7.4,1x,f6.4,3(1x,f5.1),1x,f8.3,1x,i2)') time(id),semia,ecc,Tsurf(id),T(nz,id),Tatz,zT(id),id
     endif
     
     !if (zT /= zT) stop  ! NaN
  enddo
  close(20)
  
  open(unit=30,file='depths_sphere.'//nn,action='write')
  do id=1,NB
     if (init(id)) cycle
     write(30,*) nn,id,time(id),Tsurf(id),T(nz,id),zT(id)
     do i=1,nz
        write(40,*) id,z(i),T(i,id)
     end do
  end do
  close(30)
end program sphere1d_implicit



subroutine conductionT_sphere(nz,dr,dt,T,Tsurf,Tsurfp1,kappa)
  !***********************************************************************
  !   conductionT_sphere: program to calculate the diffusion of
  !                 temperature in a spherically symmetric geometry
  !   Crank-Nicolson scheme, flux conservative
  !
  !   Eqn: T_t = (kappa/r^2)*d/dr(r^2*dT_r)
  !   BC (z=0, r=R): T=T(t)
  !   BC (z=R, r=0): center of sphere, dT_r=0
  !
  !   kappa = thermal diffusivity [m^2/s]
  !   T = radial temperature profile [K]
  !   Tsurf, Tsurfp1 = surface temperatures at times n and n+1
  !
  !   Grid: surface is at z=0
  !         T(nz) is at center of sphere
  !***********************************************************************
  implicit none
  integer, intent(IN) :: nz
  real*8, intent(IN) :: dr, dt, Tsurf, Tsurfp1, kappa
  real*8, intent(INOUT) :: T(nz)
  integer i
  real*8 cp(nz), cm(nz), cn(nz)
  real*8 alpha, a(nz), b(nz), c(nz), r(nz)
  
  alpha = kappa*dt/dr**2

  do i=1,nz-1
     r(i) = (nz-i)*dr
     cp(i) = (r(i)-dr/2)**2/r(i)**2  ! inward of r(i)
     cm(i) = (r(i)+dr/2)**2/r(i)**2  ! outward of r(i)
     cn(i) = (cp(i)+cm(i))/2.
  enddo
  
  ! elements of tridiagonal matrix
  do i=1,nz-1
     a(i) = -alpha/2*cm(i)  !  a(1) is not used
     b(i) = 1. + alpha*cn(i)
     c(i) = -alpha/2*cp(i)  !  c(nz) is not used
  enddo 
  a(nz) = -3*alpha 
  b(nz) = 1. + 3*alpha
  
  ! Set RHS         
  r(1) = alpha/2*cp(1)*T(2) + (1.-alpha*cn(1))*T(1) + alpha/2*cm(1)*(Tsurf+Tsurfp1)
  do i=2,nz-1
     r(i) = alpha/2*cp(i)*T(i+1) + (1.-alpha*cn(i))*T(i) + alpha/2*cm(i)*T(i-1)
  enddo
  r(nz) = 3*alpha*T(nz-1) + (1.-3*alpha)*T(nz)

  ! Solve for T at n+1
  call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion
  
end subroutine conductionT_sphere



subroutine retreat_s(T,zT,Radius,dt,Elatent)
  ! retreat of ice table, 1D spherically symmetric
  implicit none
  real*8, intent(IN) :: T  ! scalar temperature [K]
  real*8, intent(IN) :: Radius  ! radius of body [m]
  real*8, intent(IN) :: dt  ! time step [sec]
  real*8, intent(INOUT) :: zT  ! depth of ice table below surface
  real*8, intent(OUT) :: Elatent  ! latent heat, for diagnostics
  real*8, parameter :: pi=3.1415926535897932
  real*8, parameter :: Lh2o=2.834e6 ! latent heat of sublimation [J/kg]
  real*8, parameter :: R = 8314.5 ! universal gas constant
  real*8, parameter :: rhoice = 930.
  real*8 D, rhos, buf, diam, porosity
  real*8, external :: psv, vapordiffusivity
  real*8 zTold, dV
  
  diam=100d-3; porosity=0.5
  D=vapordiffusivity(diam,porosity,T)
  
  rhos = psv(T)*18/(R*T) ! p = n k T
  buf = D*rhos/(rhoice*porosity)/(1-zT/Radius)
  zTold = zT
  zT = sqrt(zT**2 + 2*buf*dt) ! fixes divergence at surface
  if (zT>=0. .and. zTold>=0.) then
     dV = (zT-zTold)*4*pi*(Radius-(zT+zTold)/2)**2
     !dV = 4*pi/3*((Radius-zTold)**3 - (Radius-zT)**3) ! more accurate but roundoff sensitive
     Elatent = dV*porosity*rhoice*Lh2o  ! [J]
  else
     Elatent = -999.
  end if
end subroutine retreat_s


