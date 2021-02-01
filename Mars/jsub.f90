subroutine jsub(zdepth, latitude, albedo0, thIn, pfrost, nz, &
     &     rhoc, fracIR, fracDust, patm, Fgeotherm, dt, zfac, icefrac, & 
     &     mode, Tb, avdrho, zequil)
!***********************************************************************
!  jsub: runs thermal model for Mars and returns difference in mean 
!        annual vapor density between surface and ice at depth zdepth
!
!  OUTPUTS: 
!           avdrho = difference in annual mean vapor density between
!                    ice table (or bottom) and atmosphere
!           zequil = instantaneous depth of equilibrium ice table
!
!  IN- or OUTPUT: 
!                 Tb = temperature at the bottom of the domain,
!                      negative value initializes
!
!  mode = 0  ends after fixed amount of time and then outputs all data
!  mode = 1  ends when certain accuracy is reached and outputs no data
!
!  latitude must be in RADIANS
!  Grid: surface is at z=0; T(i) is at z(i)
!  variable icefrac can also have meaning of porosity
!***********************************************************************
  implicit none
  real*8, parameter :: pi=3.1415926535897932, NULL=0.
  real*8, parameter :: earthDay=86400., marsDay=88775.244, solsperyear=668.60
  real*8, parameter :: sigSB=5.6704d-8, Lco2frost=6.0e5
  real*8, parameter :: co2albedo=0.60, co2emiss=1.

  integer, intent(IN) :: nz, mode
  real*8, intent(IN) :: zdepth, latitude, albedo0, thIn, pfrost, rhoc
  real*8, intent(IN) :: fracIR, fracDust, patm, Fgeotherm, dt, zfac, icefrac
  real*8, intent(INOUT) :: Tb
  real*8, intent(OUT) :: avdrho, zequil
  
  integer nsteps, n, i, nm
  integer iyr, imm, iday
  real*8 T(nz),tmax, time, zmax, Tsurf
  real*8 albedo, emiss0, emiss, Qn, Qnp1, tdays
  real*8 marsR, marsLs, marsDec, HA
  real*8 jd, temp1, dcor, dt0_j2000
  real*8 ti(nz), rhocv(nz), z(nz), Fsurf, m, dE
  real*8 Told(nz), Fsurfold, Tsurfold, Tmean1, Tmean2
  real*8 Tpeak(nz), Tlow(nz), Tco2frost
  real*8 rhosatav(nz), rhosatav0, rhoavs, Tbold, marsLsold, oldtime
  integer, external :: julday
  real*8, external :: flux_mars77, psv, tfrostco2, equildepth

  select case (mode)
  case (0) ! full mode
     tmax = 7020.
  case (1) ! equilibrate fast
     tmax = 40.*solsperyear  ! use integer multiple of solsperyear
  case default
     stop 'no valid mode'
  end select
  
  iyr=1996; imm=10; iday=5   ! starting year and date
  nsteps = nint(tmax/dt)      ! calculate total number of timesteps
  emiss0 = 1.   ! frost-free emissivity
  zmax = 5.*26.*thIn/rhoc*sqrt(marsDay/pi)
  !print *,'skindepth,diurnal=',thIn/rhoc*sqrt(marsDay/pi)
  !print *,'skindepth,annual=',thIn/rhoc*sqrt(solsperyear*marsDay/pi)
  
  !if (latitude>=0.) then    ! northern hemisphere
  !   Tco2frost=147.
  !else                      ! southern hemisphere
  !   Tco2frost=143.
  !endif
  Tco2frost = tfrostco2(patm)
  
  if (Tb<=0.) then
     !Tmean2=210.15          ! black-body temperature of planet
     Tmean2=(589.*(1.-albedo0)*cos(latitude)/pi/5.67e-8)**0.25 ! better estimate
  else
     Tmean2=Tb
  endif
  
  jd = dble(julday(imm,iday,iyr)) !  JD for noon UTC on iyear/imm/iday
  temp1 = (jd-2451545.d0)/36525.d0
  dcor = (64.184d0 + 95.*temp1 + 35.*temp1**2) ! correction in sec
  ! All time is referenced to dt0_j2000
  dt0_j2000 = jd + dcor/earthDay - 2451545.d0 
  !call marsorbit(dt0_j2000,0.d0,marsLs,marsDec,marsR)

  albedo = albedo0
  emiss = emiss0
  T(1:nz) = Tmean2
  do i=1,nz
     if (T(i)<Tco2frost) T(i)=Tco2frost
  enddo
  ti(1:nz) = thIn
  rhocv(1:nz) = rhoc
  Tpeak(:) = -1.e32
  Tlow(:) = +1.e32
  rhosatav(:) = 0.; rhosatav0 = 0.

  Tsurf = T(1)
  m=0.; Fsurf=0.
  nm=0; Tmean1=0.
  rhoavs=0.
  marsLsold=-1.e32; Tmean2=0.
  Tbold=-1.e32
  oldtime=1.e32
  zequil=-9999.
      
  call setgrid(nz,z,zmax,zfac)
  call smartgrid(nz,z,zdepth,thIn,rhoc,icefrac,ti,rhocv,1,NULL)
  !call smartgrid(nz,z,zdepth,thIn,rhoc,NULL,ti,rhocv,2,icefrac)
  !call smartgrid(nz,z,zdepth,thIn,rhoc,NULL,ti,rhocv,3,NULL)
  if (mode==0) write(30,'(999(f8.5,1x))') (z(i),i=1,nz)
  
  time = 0.
  tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
  call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR)
  HA = 2.*pi*time           ! hour angle
  Qn = flux_mars77(marsR,marsDec,latitude,HA,albedo,fracir,fracdust)
  !-loop over time steps 
  do n=0,nsteps-1
     time = (n+1)*dt         !   time at n+1 
     tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
     call marsorbit(dt0_j2000,tdays,marsLs,marsDec,marsR) 
     HA = 2.*pi*mod(time,1.d0) ! hour angle
         
     Qnp1 = flux_mars77(marsR,marsDec,latitude,HA,albedo,fracir,fracdust)

     Tsurfold = Tsurf
     Fsurfold = Fsurf
     Told(:) = T(1:nz)
     if (Tsurf>Tco2frost .or. m<=0.) then
        call conductionQ(nz,z,dt*marsDay,Qn,Qnp1,T,ti,rhocv,emiss, &
             &           Tsurf,Fgeotherm,Fsurf)
     endif
     if (Tsurf<Tco2frost .or. m>0.) then ! CO2 condensation
        T(1:nz) = Told(:)
        call conductionT(nz,z,dt*marsDay,T,Tsurfold,Tco2frost,ti, &
             &              rhocv,Fgeotherm,Fsurf) 
        Tsurf = Tco2frost
        dE = (- Qn - Qnp1 + Fsurfold + Fsurf + &
             &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
        m = m + dt*marsDay*dE/Lco2frost
     endif
     if (Tsurf>Tco2frost .or. m<=0.) then
        albedo = albedo0
        emiss = emiss0
     else
        albedo = co2albedo
        emiss = co2emiss
     endif
     Qn = Qnp1
     
     if (time>=tmax-solsperyear .and. mode==0) then
        do i=1,nz
           if (T(i)>Tpeak(i)) Tpeak(i)=T(i)
           if (T(i)<Tlow(i)) Tlow(i)=T(i)
           rhosatav(i) = rhosatav(i) + psv(T(i))/T(i)
        enddo
        Tmean1 = Tmean1 + Tsurf
        rhosatav0 = rhosatav0 + psv(Tsurf)/Tsurf
        rhoavs = rhoavs + min(psv(Tsurf),pfrost)/Tsurf
        nm = nm+1
     endif
     
     if (mode==1) then
        if (marsLs<marsLsold) then
           Tmean2 = Tmean2/nm
           if (abs(Tmean2-Tbold)<0.05) exit
           Tbold = Tmean2
           Tmean2=0.; nm=0
        else
           Tmean2 = Tmean2+T(nz)
           nm=nm+1
        endif
     endif
     marsLsold = marsLs
     
  enddo  ! end of time loop
  
  if (mode==0) then
     rhoavs = rhoavs/nm
     rhosatav(1:nz) = rhosatav(1:nz)/nm
     rhosatav0 = rhosatav0/nm
     !if (avrho1pre>=0.) rhoavs=avrho1pre
     Tmean1 = Tmean1/nm
     !write(31,'(999(f6.2,1x))') (Tpeak(i),i=1,nz)
     !write(32,'(999(f6.2,1x))') (Tlow(i),i=1,nz)
     !write(35,'(999(g10.4,1x))') (rhosatav(i),i=1,nz)
     if (zdepth<=0 .or. zdepth>=z(nz)) then
        avdrho = rhosatav(nz)-rhoavs ! no ice
     else
        do i=1,nz
           if (z(i)>zdepth) then ! ice
              avdrho = rhosatav(i)-rhoavs
              !if (i==1) avdrho = zint(0.d0,z(1),rhosatav0,rhosatav(1)) - rhoavs
              !if (i>1) then
              !   avdrho = zint(z(i-1),z(i),rhosatav(i-1),rhosatav(i)) - rhoavs
              !endif
              exit
           endif
        enddo
     endif

     zequil = equildepth(nz, z, rhosatav, rhosatav0, rhoavs)
     
     write(34,'(f7.3,1x,2(f6.2,1x),2(g10.4,1x))') &
          &        Tmean1,Tlow(nz),Tpeak(nz),avdrho+rhoavs,rhoavs
  endif
  Tb = T(nz)
end subroutine jsub



function zint(y1,y2,z1,z2)
  ! linearly interpolate between two values
  implicit none
  real*8 zint
  real*8, intent(IN) :: y1, y2, z1, z2
  zint = (y1*z2 - y2*z1)/(y1-y2)  ! yint = 0
end function zint



function equildepth(nz, z, rhosatav, rhosatav0, avrhoatm)
!***********************************************************************
! returns equilibrium depth for given thermal properties
! find depth where rhosatav = avrhoatm	     
! this is not the true (final) equilibrium depth because of tifeedback
!***********************************************************************
  implicit none
  real*8 equildepth
  integer, intent(IN) :: nz
  real*8, intent(IN) :: z(nz), rhosatav(nz), rhosatav0, avrhoatm
  integer i, kE
  real*8 zdepthE
  real*8, external :: zint
  
  zdepthE = -9999.  ! unstable at all depths
  kE = -9

  ! find shallowest point where ice is stable
  do i=1,nz
     if (rhosatav(i) <= avrhoatm) then
        kE = i
        exit
     endif
  enddo

  !if (kE>=1) zdepthE = z(kE)
  if (kE>1) then  ! interpolate
     zdepthE = zint( avrhoatm-rhosatav(kE-1), avrhoatm-rhosatav(kE), &
          & z(kE-1), z(kE) )
  endif
  if (kE==1) zdepthE = zint(avrhoatm-rhosatav0,avrhoatm-rhosatav(1),0.d0,z(1))
  if (zdepthE>z(nz)) zdepthE = -9999. ! just in case

  equildepth = zdepthE
end function equildepth
