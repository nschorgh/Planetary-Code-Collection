!***************************************************************
! Growth of pore ice with fast method and nothing else
! no massive ice, no retreat
!
! Written by Norbert Schorghofer 2007-2009
!***************************************************************


program stabgrow
  use miscparameters
  implicit none
  integer nz, i, typeG, typeF
  real(8) zfac, zmax, z(NMAX), icetime, porosity
  real(8) pfrost   ! partial pressure of H2O
  real(8) thIn, rhoc, timestep, porefill(NMAX), zdepthF
  real(8) Tmean, Tampl, Diff, tmax, B, zdepthG
  real(8) junk(NMAX)
  namelist /stuff/ nz,zmax,zfac,thIn,rhoc,pfrost,Tmean,Tampl,Diff, &
       & porosity,tmax,timestep
  real(8), external :: colint  ! defined in allinterfaces.mod

  ! parameters that do not change
  open(unit=20,file='input_fast.par')
  read(20, nml=stuff)
  write(6,nml=stuff)

  ! set eternal grid
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z',action='write',status='unknown')
  write(30,'(999(f8.5,1x))') z(1:nz)
  close(30)

  icetime=0.
  porefill(1:nz) = 0.
  write(36,*) icetime,porefill(1:nz)
  zdepthF = zmax  ! for first evaluation of d z_F/dt

  ! equilibrate initial temperature
  ! skips use of soilthprop
  print *,'equilibrating initial temperature'
  B = Diff*timestep*86400.*365.24/porosity/icedensity
  do i=1,10
     zdepthF = -9999.
     call ajsub2(nz, z, spread(thIn,1,NMAX), spread(rhoc,1,NMAX), pfrost, Tmean, Tampl, &
          & typeF, zdepthF, junk(1:nz), porefill(1:nz), .TRUE., B, typeG, zdepthG)
  enddo

  do i=1,10000000
     icetime = i*timestep
     if (icetime>tmax) exit
     call icelayer2(timestep,nz,thIn,rhoc,z,porosity,pfrost,porefill(1:nz), &
          & zdepthF,Tmean,Tampl,Diff,zdepthG)
     if (abs(mod(icetime/2000.,1.d0))<1.e-3) then ! output every 2000 years
        write(36,*) icetime,porefill(1:nz)
     endif
     write(35,*) icetime,zdepthF,colint(porefill(1:nz),z,nz,1,nz) 
     print *,icetime,zdepthF
  enddo
end program stabgrow



subroutine diagnoselatent(nz,rhof,rhocv)
  use miscparameters, only : NMAX
  implicit none
  real(8), parameter :: Lh2o=2.83e6
  integer, intent(IN) :: nz
  real(8), intent(IN) :: rhof(nz), rhocv(nz)
  real(8), SAVE :: rhofold(NMAX)
  integer j
  do j=1,nz
     if (j>1) print *,j,(rhof(j)-rhofold(j))*Lh2o/(rhocv(j)+rhocv(j-1))*2.
     rhofold(j) = rhof(j)
  enddo
end subroutine diagnoselatent



subroutine icelayer2(bigstep, nz, thIn, rhoc, z, porosity, pfrost, &
     & porefill, zdepthF, Tmean, Tampl, Diff, zdepthG)
!***********************************************************************
! a.k.a. icelayer_exper_poreonly
! doesn't do retreat
!***********************************************************************
  use miscparameters, only : NMAX, icedensity
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: bigstep, thIn, rhoc, z(NMAX), porosity, pfrost
  real(8), intent(IN) :: Tmean, Tampl, Diff
  real(8), intent(INOUT) :: porefill(nz), zdepthF
  real(8), intent(OUT) :: zdepthG
  integer j, typeF, typeG
  real(8) ti(NMAX), rhocv(NMAX), ypp(nz), fill, B
  real(8), parameter :: NULL = 0.

  B = Diff*bigstep*86400.*365.24/porosity/icedensity

  ! assign thermal properties of soil
  ti(1:nz) = thIn
  rhocv(1:nz) = rhoc
  do j=1,nz
     if (j==1) then
        fill = porefill(1)/2.  ! porefill(0)=0
     else
        !fill = (porefill(j)+porefill(j-1))/2.  ! tends to be unstable
        fill = porefill(j)  ! tends to be stable
     endif
     call soilthprop(porosity,fill,rhoc,thIn,1,rhocv(j),ti(j),NULL)
  enddo
  
  ! run thermal model
  call ajsub2(nz, z, ti, rhocv, pfrost, Tmean, Tampl, typeF, zdepthF, &
       & ypp, porefill(:), .FALSE., B, typeG, zdepthG)

  ! diffusive filling
  if (typeF>0) then
     do j=1,typeF-1  ! works also for typeF=1
        porefill(j) = 0.
     enddo
     do j=typeF,nz
        porefill(j) = porefill(j) + B*ypp(j)
        if (porefill(j)<0.) porefill(j)=0.
        if (porefill(j)>1.) porefill(j)=1.
     enddo
  endif
  !rhof = porefill(:)*porosity*icedensity
  !call diagnoselatent(nz,rhof,rhocv(:))

  ! enact bottom boundary 
  if (typeG>0) porefill(typeG:nz)=0.

end subroutine icelayer2



subroutine ajsub2(nz, z, ti, rhocv, pfrost, Tmean, Tampl, &
     &     typeF, zdepthF, ypp, porefill, init, B, typeG, zdepthG)
  use miscparameters, only : NMAX, solsperyear, marsDay
  use thermalmodelparam_exper  ! defined in fast_subs_exper
  !use allinterfaces
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(NMAX), ti(NMAX), rhocv(NMAX), pfrost
  real(8), intent(IN) :: Tmean, Tampl, B
  logical, intent(IN) :: init
  integer, intent(OUT) :: typeF  ! index of depth below which filling occurs
  real(8), intent(INOUT) :: zdepthF 
  real(8), intent(OUT) :: ypp(nz) ! (d rho/dt)/Diff
  real(8), intent(IN) :: porefill(nz)
  integer, intent(OUT) :: typeG
  real(8), intent(OUT) :: zdepthG

  integer nsteps, n, i, nm
  real(8) tmax, time, Tsurfp1, Fsurf
  real(8), SAVE :: T(NMAX), Tsurf
  real(8) avrho1, rhosatav0, rhosatav(nz), rlow 
  ! real(8) yp(nz)
  real(8), external :: Tsurface, psv
  
  tmax = EQUILTIME*solsperyear
  nsteps=int(tmax/dt)       ! calculate total number of timesteps

  if (init) T(1:nz) = Tmean
  Tsurf = T(1)

  nm=0
  avrho1=0.
  rhosatav0 = 0.
  rhosatav(:) = 0.

  time=0.
  !-----loop over time steps 
  do n=0,nsteps-1
     time =(n+1)*dt         !   time at n+1 
     Tsurfp1 = Tsurface(time,Tmean,Tampl)
     call conductionT(nz,z,dt*marsDay,T,Tsurf,Tsurfp1,ti,rhocv,Fgeotherm,Fsurf) 
     Tsurf=Tsurfp1
     
     if (time>=tmax-solsperyear) then
        avrho1 = avrho1+min(psv(Tsurf),pfrost)/Tsurf
        rhosatav0 = rhosatav0+psv(Tsurf)/Tsurf
        do i=1,nz
           rhosatav(i) = rhosatav(i)+psv(T(i))/T(i)
        enddo
        nm=nm+1
     endif

  enddo  ! end of time loop
  
  avrho1 = avrho1/nm
  rhosatav0 = rhosatav0/nm; rhosatav(:)=rhosatav(:)/nm

  if (Fgeotherm>0.) then
     typeG=1   ! will be overwritten by avmeth
     rlow=2*rhosatav(nz)-rhosatav(nz-1)
  else
     typeG=-9
     rlow=rhosatav(nz-1)
  end if
  !call deriv1(z,nz,rhosatav,rhosatav0,rlow,yp) 

  call depths_avmeth(-9, nz, z(1:nz), rhosatav(:), rhosatav0, rlow, avrho1, &
       & porefill(:), typeF, zdepthF, B, ypp(:), typeG, zdepthG)

end subroutine ajsub2


