PROGRAM testcrankQ
!***********************************************************************
! test Crank-Nicolson subroutine
!***********************************************************************
  !use conductionQ   ! for test of conductionQ2.f90 
  implicit none
  integer, parameter :: nz=60
  real*8, parameter :: pi=3.1415926535897931, d2r=pi/180.
  real*8, parameter :: sigSB=5.6704d-8, zero=0.
  real*8, parameter :: Period = 88775.244*670, emiss=1., Fgeo = 0.2

  integer n, i, STEPSPERSOL
  real*8 T(nz), time, dt, zmax, zfac
  real*8 latitude, albedo
  real*8 thIn  ! thermal inertia
  real*8 rhocv(nz), Qn, Qnp1
  real*8 Rau, Decl, HA
  real*8 flux_noatm, delta, ti(nz), Fsurf, Tsurf, z(nz), Fmean
  real*8 Tmean(nz), Hflux(nz), kcond(nz)
  external flux_noatm
  
  STEPSPERSOL = 120
  dt = Period / STEPSPERSOL
  thIn = 120.
  albedo = 0.2
  latitude = 5.  ! [degree]
  zmax = 2.5; zfac = 1.05

  ! rhoc=thIn*sqrt(Period/pi)  ! skin depth = 1
  rhocv(:) = 1200.*800.
  delta=thIn/rhocv(1)*sqrt(Period/pi)
  print *,'Skin depth= ',delta
  ti(1:nz) = thIn
  kcond(:) = ti(1:nz)**2/rhocv(1:nz)
  Tmean(:) = 0.
  
  Rau = 1.52; Decl = 0.
  
  print *,'Time step=',dt
  !print *,'zmax=',(nz-0.5)*dz
  print *,'zmax=',zmax
  print *,'Thermal inertia=',thIn,' Period=',Period
  print *,'Heat Flux at bottom boundary=',-Fgeo

  ! Initialize
  open(unit=21,file='Tsurface',status='unknown') ! surface temperature
  open(unit=22,file='Tprofile',status='unknown') ! temperature profile
  
  Tsurf = 210.
  do i=1,nz
     T(i) = 210.
     !z(i) = (i-0.5)*dz
  enddo
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z',status='unknown');
  write(30,*) (z(i),i=1,nz)
  close(30)
  
  latitude = latitude*d2r
  
  !write(22,'(f8.4,1000(x,f7.2))') time,(T(i),i=1,nz)
  
  do i=1,nz
     !if (i>nz/4) ti(i)=2.*ti(i)
  enddo
  Fmean = 0.
  
  HA = 0.
  
! solar insolation 
  Qn = (1-albedo)*flux_noatm(Rau,Decl,latitude,HA,zero,zero)
  ! Qn=sigSB*200.**4

  ! call conductionQ2_init(nz,z,dt,ti,rhocv,Fgeo)
      
  do n=0,50000
     !print *,time/Period,Qn
     time = (n+1)*dt   !   time at n+1; 
     HA = 2*pi*mod(time/Period,1.d0) !  hour angle
     Qnp1 = (1-albedo)*flux_noatm(Rau,Decl,latitude,HA,zero,zero)
!    Qnp1=sigSB*200.**4
!    call crankn1(nz,Qn,Qnp1,T,1.,a,b,c,alpha,k1)
!    call crankn1(nz,dz,dt,Qn,Qnp1,T,thIn,rhoc,1.d0)
!    call cranknv(nz,dz,dt,Qn,Qnp1,T,ti,rhoc,1.d0)
     call conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeo,Fsurf)
!    call conductionQ2(nz,Qn,Qnp1,T,emiss,Tsurf,Fsurf)
     Qn = Qnp1

     if (mod(n,3)==0) then
        write(21,'(f12.6,2f9.3)') time/Period,Tsurf,T(nz)
     endif
     if (n > 50000-STEPSPERSOL) then
        if (mod(n,10)==0) then
           write(22,'(1000(x,f7.2))') Tsurf,(T(i),i=1,nz)
        endif
        Fmean = Fsurf + Fmean
        Tmean(:) = Tmean(:) + T(:)
     endif
     
  end do                    ! end of the loop
  
  Fmean = Fmean/STEPSPERSOL
  Tmean(:) = Tmean(:)/STEPSPERSOL
  call heatflux_from_temperature(nz,z,Tmean,kcond,Hflux)

  print *,'# Depth           Temperature         Heat_Flux'
  !print *,0,Tsurf,Fmean
  do i=1,nz
     print *,z(i),T(i),Hflux(i)
  enddo
  
  close(21)
  close(22)
END PROGRAM testcrankQ
