PROGRAM testcrankQ_conv
!***********************************************************************
! test Crank-Nicolson subroutine convergence with time step
!***********************************************************************
  implicit none
  integer, parameter :: nz=30
  real*8, parameter :: pi=3.1415926535897931, d2r=pi/180.
  real*8, parameter :: sigSB=5.6704d-8, zero=0.
  real*8, parameter :: Period=88775.244*670, emiss=1., Fgeo=0.0

  integer n, i, STEPSPERSOL, kk
  real*8 T(nz), time, dt, zmax, zfac
  real*8 latitude, albedo
  real*8 thIn  ! thermal inertia
  real*8 rhocv(nz), Qn, Qnp1
  real*8 Rau, Decl, HA
  real*8 flux_noatm, delta, ti(nz), Fsurf, Tsurf, z(nz)
  external flux_noatm
  
  thIn = 120.
  albedo = 0.2
  latitude = 5.  ! [degree]
  zmax = 4.; zfac = 1.05

  rhocv = thIn*sqrt(Period/pi)  ! skin depth = 1
  !rhocv(:) = 1200.*800.
  delta = thIn/rhocv(1)*sqrt(Period/pi)
  print *,'Skin depth= ',delta
  ti(1:nz) = thIn
  
  Rau = 1.52; Decl = 0.
  
  print *,'Time step=',dt
  !print *,'zmax=',(nz-0.5)*dz
  print *,'zmax=',zmax
  print *,'Thermal inertia=',thIn,' Period (days)',Period/86400.
  print *,'Heat Flux at bottom boundary=',-Fgeo

  ! Initialize
  open(unit=21,file='Tsurface',status='unknown') ! surface temperature
  open(unit=22,file='Tprofiles',status='unknown') ! temperature profile

  call setgrid(nz,z,zmax,zfac)
  
  latitude = latitude*d2r
  
  do kk = 0,7
     STEPSPERSOL = 2**(12-kk)
     dt = Period / STEPSPERSOL
     print *,'Steps per sol',STEPSPERSOL,'dt=',dt/86400.
  
     Tsurf = 210.
     do i=1,nz
        T(i) = 210.
        !z(i) = (i-0.5)*dz
     enddo
  
     HA = 0.
     
     ! solar insolation 
     Qn = (1-albedo)*flux_noatm(Rau,Decl,latitude,HA,zero,zero)
      
     do n=0,512*STEPSPERSOL-1
        !print *,time/Period,Qn
        time = (n+1)*dt   !   time at n+1; 
        HA = 2*pi*mod(time/Period,1.d0) !  hour angle
        Qnp1 = (1-albedo)*flux_noatm(Rau,Decl,latitude,HA,zero,zero)
        call conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeo,Fsurf)
        Qn = Qnp1

        if (mod(n,3)==0) then
           write(21,'(f12.6,2f9.3)') time/Period,Tsurf,T(nz)
        endif

     end do                    ! end of the loop

     print *,'Time (days)',time/86400.
     write(22,*) 0.,Tsurf
     do i=1,nz
        write(22,*) z(i),T(i)
     enddo

  end do
  close(21)
  close(22)
END PROGRAM testcrankQ_conv
   
