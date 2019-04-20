module newhorizons
  ! privately stores horizon data
  use filemanager, only : NSx, NSy, sfn
  implicit none
  integer, parameter :: naz=180  ! # azimuth values
  real(8), dimension(NSx,NSy,naz+1), private :: s ! angles

contains
  subroutine readhorizons
    ! reads and stores all horizon elevations
    ! also converts slopes to angles
    implicit none
    integer i,j,ia,ja,ierr
      
    open(unit=20,file=sfn,status='old',action='read',iostat=ierr)
    if (ierr>0) then
       print *,sfn
       stop 'readhorizons: Input file not found'
    endif
    do i=2,NSx-1
       do j=2,NSy-1
          read(20,*) ia,ja,s(i,j,1:naz)
          if (ia /= i .or. ja /= j) then
             print *,i,j,ia,ja
             stop 'readhorizons: index mismatch'
          endif
          s(i,j,naz+1) = s(i,j,1) !360 deg
       enddo
    enddo
    s = atan(s)  ! slope -> angle
    close(20)
  end subroutine readhorizons
    
  elemental function getonehorizon(i0,j0,azSun)
    ! returns elevation of horizon along azimuth azSun
    implicit none
    integer, intent(IN) :: i0,j0
    real(8), intent(IN) :: azSun
    real(8) getonehorizon  ! angle
    real(8), parameter :: pi=3.1415926535897932
    integer k
    real(8) a, daz, smax
    
    daz = 2*pi/real(naz)
    k = floor(modulo(azSun,2*pi)/daz)
    a = modulo(azSun,2*pi)/daz-k
    !if (k<0 .or. k>=naz) then
       !print *,'azSun=',azSun,'k=',k
       !stop 'getonehorizon: impossible k value'  ! impure
    !endif
    smax = s(i0,j0,k+1)*(1.-a) + s(i0,j0,k+2)*a
    getonehorizon = smax
  end function getonehorizon
  
  elemental function getoneskysize(i,j)
    !***********************************************************************
    !   calculates sky size (steradian) from horizons
    !***********************************************************************
    use allinterfaces, only: area_spherical_triangle
    implicit none
    real(8) getoneskysize
    integer, intent(IN) :: i,j
    real(8), parameter :: pi=3.1415926535897932
    integer k
    real(8) phi(3), theta(3), dOmega, skysize

    skysize=0
    phi = (/ 0, 0, 1 /) *2*pi/naz
    do k=1,naz
       theta = (/ 0.d0, atan(1/s(i,j,k)), atan(1/s(i,j,k+1)) /) ! from zenith
       dOmega = area_spherical_triangle(phi,theta)
       skysize = skysize + dOmega
    enddo
    getoneskysize = skysize
    
  end function getoneskysize
  
  elemental function getoneskysize_v2(i,j)
    !***********************************************************************
    !   calculates sky size (steradian) from horizons
    !***********************************************************************
    implicit none
    real(8) getoneskysize_v2
    integer, intent(IN) :: i,j
    real(8), parameter :: pi=3.1415926535897932
    real(8), parameter :: dphi=2*pi/naz
    real(8) landsize
    
    landsize = sum(sin(s(i,j,1:naz)))*dphi
    getoneskysize_v2 = 2*pi-landsize
    
  end function getoneskysize_v2

  elemental function getoneGterm(i,j,alpha,azFac)
    ! quantity used for specific type of approximation
    implicit none
    real(8) getoneGterm
    integer, intent(IN) :: i,j
    real(8), intent(IN) :: alpha, azFac
    real(8), parameter :: pi=3.1415926535897932
    integer k
    real(8), dimension(1:naz) :: ssub, az, daz, emin
    real(8) G1, G2

    ssub(1:naz) = s(i,j,1:naz)
    do concurrent(k=1:naz)
       az(k) = real(2*pi*(k-1),8)/naz
       daz(k) = azFac-az(k)+pi  ! the +pi was two days of work
       ! daz should be zero when line of sight coincides with direction of steepest descent
       emin(k) = atan( -tan(alpha)*cos(daz(k)) ) 
       if (emin(k)>ssub(k)) ssub(k)=emin(k) ! self-shadowing higher than distant horizon
    end do
    G1 = sum(sin(ssub)**2) - sum(sin(emin)**2)
    G2 = sum(cos(daz)*( ssub+sin(ssub)*cos(ssub) - emin-sin(emin)*cos(emin) ))
    G1 = G1/naz 
    G2 = G2/naz
    
    getoneGterm = cos(alpha)*G1 + sin(alpha)*G2
    ! Gterm should never be negative
  end function getoneGterm
  
end module newhorizons



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
!     decl: planetocentric solar declination [radians]
!     latitude [radians]
!     HA: hour angle [radians from noon, clockwise]
!     sinbeta: sin(altitude)
!     azimuth [radians east of north]
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
  if (buf>+1.) buf=+1.  ! roundoff
  if (buf<-1.) buf=-1.  ! roundoff
  azimuth = acos(buf)
  if (sin(HA)>=0) azimuth=2*pi-azimuth  
end subroutine equatorial2horizontal



subroutine getfieldofview(NSx,NSy,ffn,cc,ia,ja,dOh,skysize,CCMAX)
  implicit none
  integer, intent(IN) :: NSx, NSy
  character(len=*), intent(IN) :: ffn
  integer, intent(IN) :: CCMAX
  integer, intent(OUT) :: cc(NSx,NSy) ! number of cells in field of view
  integer(2), intent(OUT), dimension(NSx,NSy,CCMAX) :: ia, ja
  real(4), intent(OUT), dimension(NSx,NSy,CCMAX) :: dOh
  real(8), intent(OUT) :: skysize(NSx,NSy)
  integer i, j, k, i0_2, j0_2, ierr

  open(unit=20,file=ffn,status='old',action='read',iostat=ierr)
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



integer function getmaxfieldsize(NSx,NSy,ffn)
  implicit none
  integer, intent(IN) :: NSx,NSy
  character(len=*), intent(IN) :: ffn
  integer maxsize
  integer cc, i, j, i0_2, j0_2, ierr

  open(unit=20,file=ffn,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getmaxfieldsize: input file not found'

  maxsize=0
  do i=2,NSx-1
     do j=2,NSy-1
        read(20,'(2(i5,1x),i6,1x)') i0_2,j0_2,cc
        if (i/=i0_2 .or. j/=j0_2) stop 'getmaxfieldsize: wrong data order'
        if (cc>maxsize) maxsize=cc
     enddo
  enddo
  close(20)

  getmaxfieldsize = maxsize
end function getmaxfieldsize


