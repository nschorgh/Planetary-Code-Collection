subroutine gethorizon(i0,j0,azSun,smax,first)
  use filemanager
  implicit none
  integer, intent(IN) :: i0,j0
  real(8), intent(IN) :: azSun
  real(8), intent(OUT) :: smax
  logical, intent(IN) :: first 
  integer i,j,ia,ja,k
  integer, parameter :: nres=360  ! # azimuth values
  real(8), parameter :: pi=3.1415926535897932
  real(8), dimension(NSx,NSy,nres+1), save :: s
  real(8) a, daz

  if (first) then
     open(unit=20,file='horizons.'//fileext,status='old',action='read')
     do i=2,NSx-1
        do j=2,NSy-1
           read(20,*) ia,ja,s(i,j,1:nres)
           if (ia /= i .or. ja /= j) then
              print *,i,j,ia,ja
              stop 'gethorizon: inconsistent input values'
           endif
           s(i,j,nres+1) = s(i,j,1) !360 deg
        enddo
     enddo
     s = atan(s)  ! slope -> angle
  else ! interpolate
     daz = 2*pi/real(nres)
     k = floor(modulo(azSun,2*pi)/daz)
     a = modulo(azSun,2*pi)/daz-k
     if (k<0 .or. k>=nres) then
        print *,'Error in gethorizon: impossible k value'
     !   print *,'azSun=',azSun,'k=',k
     !   stop 'gethorizon: impossible k value'
     endif
     smax = s(i0,j0,k+1)*(1.-a) + s(i0,j0,k+2)*a
     !smax = atan(smax)   ! slope -> angle
  endif
end subroutine gethorizon



pure subroutine difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)
  ! calculate slopes and azimuths of surface elements
  implicit none
  integer, intent(IN) :: NSx,NSy
  real(8), intent(IN) :: h(NSx,NSy),dx,dy
  real(8), intent(OUT), dimension(NSx,NSy) :: surfaceSlope,azFac
  integer i,j
  real(8) sx,sy

  do i=2,NSx-1
     do j=2,NSy-1
        sx=(h(i+1,j)-h(i-1,j))/(2.*dx)
        sy=(h(i,j+1)-h(i,j-1))/(2.*dy)
        surfaceSlope(i,j)=atan(sqrt(sx**2+sy**2))
        azFac(i,j)=atan2(sx,-sy)  ! north is up, clockwise (y increases equatorward)
     enddo
  enddo
end subroutine difftopo



elemental subroutine equatorial2horizontal(decl,latitude,HA,sinbeta,azimuth)
!***********************************************************************
!   converts declination and hour angle to azimuth and sin(altitude)
!     decl: planetocentric solar declination (radians)
!     latitude: (radians)
!     HA: hour angle (radians from noon, clockwise)
!     sinbeta: sin(altitude)
!     azimuth: (radians east of north)
!***********************************************************************
  implicit none
  real(8), parameter :: pi=3.1415926535897931
  real(8), intent(IN) :: decl,latitude,HA
  real(8), intent(OUT) :: sinbeta,azimuth
  real(8) c1,s1,cosbeta,buf
  
  c1=cos(latitude)*cos(decl)

  s1=sin(latitude)*sin(decl)
  ! beta = 90 minus incidence angle for horizontal surface
  ! beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1
  
  cosbeta = sqrt(1-sinbeta**2)
  buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta)
  ! buf can be NaN if cosbeta = 0
  if (buf>+1.) buf=+1.;  ! damn roundoff
  if (buf<-1.) buf=-1.;  ! damn roundoff
  azimuth = acos(buf)
  if (sin(HA)>=0) azimuth=2*pi-azimuth  
end subroutine equatorial2horizontal



elemental real(8) function flux_wgeom(R,sinbeta,azSun,surfaceSlope,azFac,smax)
!***********************************************************************
!   flux:  program to calculate incoming solar flux without atmosphere
!     R: distance from sun (AU)
!     sinbeta: sin(altitude) 
!     azSun: azimuth of Sun (radians east of north)
!     surfaceSlope: >=0, (radians) 
!     azFac: azimuth of topographic gradient (radians east of north)
!     smax: maximum slope in direction of azimuth
!***********************************************************************
  implicit none
  real(8), parameter :: So=1365.  ! solar constant
  real(8), intent(IN) :: R,azSun,sinbeta,surfaceSlope,azFac,smax
  real(8) cosbeta,sintheta
  
  cosbeta = sqrt(1.-sinbeta**2)

!-incidence angle
  ! theta = 90 minus incidence angle for sloped surface
  sintheta = cos(surfaceSlope)*sinbeta - &
       &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)
  if (cosbeta==0.) sintheta = cos(surfaceSlope)*sinbeta ! doesn't use azimuths

!-shadowing
  sintheta = max(sintheta,0.d0)  ! self-shadowing
  if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity
  !if (sinbeta<smax/sqrt(1.+smax**2)) sintheta=0.  ! shadowing, tan -> sin 
  if (sinbeta<sin(smax)) sintheta=0. 

!-intensity
  flux_wgeom = sintheta*So/(R**2)
end function flux_wgeom



subroutine getfieldofview(NSx,NSy,fileext,cc,ia,ja,dOh,skysize,CCMAX)
  implicit none
  integer, intent(IN) :: NSx, NSy
  character(len=*), intent(IN) :: fileext
  integer, intent(IN) :: CCMAX
  integer, intent(OUT) :: cc(NSx,NSy) ! number of cells in field of view
  integer(2), intent(OUT), dimension(NSx,NSy,CCMAX) :: ia, ja
  real(4), intent(OUT), dimension(NSx,NSy,CCMAX) :: dOh
  real(8), intent(OUT) :: skysize(NSx,NSy)
  integer i, j, k, i0_2, j0_2, ierr

  open(unit=20,file='fieldofviews.'//fileext,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getfieldofview: input file not found'
  do i=2,NSx-1
     do j=2,NSy-1
        read(20,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') & ! format must match
             & i0_2,j0_2,cc(i,j),skysize(i,j) 
        if (i/=i0_2 .or. j/=j0_2) stop 'getfieldofview: wrong data order'
        if (cc(i,j)>CCMAX) stop 'getfieldofview: not enough memory allocated'
        do k=1,cc(i,j)
           read(20,'(2(i5,1x),g10.4,1x)',advance='no') & ! format must match
                & ia(i,j,k),ja(i,j,k),dOh(i,j,k)
        enddo
        read(20,'()')
     enddo
  enddo
  close(20)
end subroutine getfieldofview



subroutine getmaxfieldsize(NSx,NSy,fileext,maxsize,type)
  implicit none
  integer, intent(IN) :: NSx,NSy,type
  character(len=*), intent(IN) :: fileext
  integer, intent(OUT) :: maxsize
  integer cc, i, j, i0_2, j0_2, ierr

  if (type==1) then
     open(unit=20,file='fieldofviews.'//fileext,status='old',action='read',iostat=ierr)
     if (ierr>0) stop 'getmaxfieldsize: input file not found'
  elseif (type==2) then
     open(20,file='inversefoviews.'//fileext,status='old',action='read',iostat=ierr)
     if (ierr>0) stop 'getmaxfieldsize: input file not found'
  else
     stop 'no valid type'
  end if
  maxsize=0
  do i=2,NSx-1
     do j=2,NSy-1
        read(20,'(2(i5,1x),i6,1x)') i0_2,j0_2,cc
        if (i/=i0_2 .or. j/=j0_2) stop 'getmaxfieldsize: wrong data order'
        if (cc>maxsize) maxsize=cc
     enddo
  enddo
  close(20)
end subroutine getmaxfieldsize



subroutine getinversefieldofview(NSx,NSy,fileext,cc,ia,ja,dO2,CCMAX)
  implicit none
  integer, intent(IN) :: NSx, NSy
  character(len=*), intent(IN) :: fileext
  integer, intent(IN) :: CCMAX
  integer, intent(OUT) :: cc(NSx,NSy)
  integer(2), intent(OUT), dimension(NSx,NSy,CCMAX) :: ia, ja
  real(4), intent(OUT), dimension(NSx,NSy,CCMAX) :: dO2
  real(8) mysize
  integer i, j, k, i0_2, j0_2, ierr

  open(20,file='inversefoviews.'//fileext,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getinversefiledofview: input file not found'
  do i=2,NSx-1
     do j=2,NSy-1
        read(20,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') & ! format must match
             & i0_2,j0_2,cc(i,j),mysize
        if (i/=i0_2 .or. j/=j0_2) stop 'getinversefoviews: wrong data order'
        if (cc(i,j)>CCMAX) stop 'getfieldofview: not enough memory allocated'
        do k=1,cc(i,j)
           read(20,'(2(i5,1x),g10.4,1x)',advance='no') & ! format must match
                & ia(i,j,k),ja(i,j,k),dO2(i,j,k)
        enddo
        read(20,'()')
     enddo
  enddo
  close(20)
end subroutine getinversefieldofview



subroutine subsurfaceconduction(T,Tsurf,dtsec,Qn,Qnp1,emiss)
  ! 1d subsurface conduction
  use allinterfaces, only : conductionQ
  use filemanager, only : solarDay, Fgeotherm, nz
  implicit none
  integer, parameter :: NMAX=1000, Ni=5
  real(8), parameter :: pi=3.1415926535897932
  real(8), intent(INOUT) :: T(NMAX), Tsurf
  real(8), intent(IN) :: dtsec,Qn,Qnp1,emiss
  integer i, k
  real(8) zmax, zfac, Fsurf, Tinit, delta
  real(8) Tsurfold, Told(1:nz), Qarti, Qartiold 
  real(8), save :: ti(NMAX), rhocv(NMAX), z(NMAX)
  logical, save :: first = .true.

  if (first) then ! initialize grid
     !Earth
     !ti(:) = 1000.;  rhocv(:) = 1200.*800.  ! adjust
     !zmax=2.; zfac = 1.05  ! adjust

     !Moon
     ti(:) = 100.;  rhocv(:) = 1200.*800.  ! adjust
     zmax=0.5; zfac = 1.05  ! adjust

     delta = ti(1)/rhocv(1)*sqrt(solarDay/pi)  ! skin depth

     call setgrid(nz,z,zmax,zfac)
     if (z(6)>delta) then
        print *,'WARNING: less than 6 points within diurnal skin depth'
     endif
     do i=1,nz
        if (z(i)<delta) cycle
        print *,i-1,' grid points within diurnal skin depth'
        exit
     enddo
     if (z(1)<1.e-5) print *,'WARNING: first grid point is too shallow'
     open(unit=30,file='z',status='unknown');
     write(30,*) (z(i),i=1,nz)
     close(30)

     write(*,*) 'Subsurface model parameters'
     write(*,*) '   nz=',nz,' zmax=',zmax,' zfac=',zfac
     write(*,*) '   Thermal inertia=',ti(1),' rho*c=',rhocv(1)
     print *,'   Diurnal skin depth=',delta,' Geothermal flux=',Fgeotherm

     first = .false.
  endif
  
  if (Tsurf<=0.) then  ! initialize temperature profile
     if (Tsurf==0.) Tinit=273.
     if (Tsurf<0.) Tinit=-Tsurf
     T(1:nz) = Tinit
     Tsurf = Tinit
  endif

  Tsurfold=Tsurf; Told(1:nz)=T(1:nz)
  call conductionQ(nz,z,dtsec,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeotherm,Fsurf)
  
  ! fixes rapid sunrise on sloped surface
  ! artificial flux smoothing
  if (Tsurf>2*Tsurfold .or. Tsurf<Tsurfold/2) then
     Tsurf = Tsurfold; T(1:nz) = Told(1:nz)
     do k=1,Ni
        Qartiold = ((Ni-k+1)*Qn + (k-1)*Qnp1)/real(Ni)
        Qarti = ((Ni-k)*Qn + k*Qnp1)/real(Ni)
        call conductionQ(nz,z,dtsec/Ni,Qartiold,Qarti,T,ti,rhocv,emiss,Tsurf,Fgeotherm,Fsurf)
     enddo
  endif
  
end subroutine subsurfaceconduction
