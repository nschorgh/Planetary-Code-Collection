MODULE newhorizons
  ! reads and privately stores horizon data
  ! and contains functions that calculate values from horizons
  use filemanager, only : NSx, NSy, sfn
  implicit none
  integer, parameter :: naz=180  ! # azimuth values
  real(8), dimension(:,:,:), private, allocatable :: s 

contains
  subroutine readhorizons(Mx1,Mx2,My1,My2)
    ! reads and stores all horizon elevations [radians]
    ! also converts slopes to angles
    implicit none
    integer, intent(IN), optional :: Mx1,Mx2,My1,My2
    logical crop
    integer i,j,ia,ja,ierr,naz2
    real(8) sread(naz)
    integer, external :: countcolumns

    naz2=countcolumns()-2
    if (naz/=naz2) then
      print *,'inconsistent number of azimuth rays',naz,naz2
      stop
    endif
    print *,'number of azimuth rays= ',naz

    if (present(Mx1).and.present(Mx2).and.present(My1).and.present(My2)) then
       crop = .true.
    else
       crop = .false.
    endif
    
    if (.not.crop) then ! the whole domain
       allocate(s(NSx,NSy,naz+1))
    else  ! cropped domain
       if (Mx1<=1 .or. Mx2>=NSx .or. My1<=1 .or. My2>=NSy) then
          stop 'readhorizons: Inappropriate size for region of interest'
       endif
       allocate(s(Mx1:Mx2,My1:My2,naz+1))
    endif
    
    open(unit=20,file=sfn,status='old',action='read',iostat=ierr)
    if (ierr>0) then
       print *,sfn
       stop 'readhorizons: Input file not found'
    endif
    do i=2,NSx-1
       do j=2,NSy-1
          !read(20,*) ia,ja,s(i,j,1:naz)
          read(20,*) ia,ja,sread(:)
          if (ia /= i .or. ja /= j) then
             print *,i,j,ia,ja
             stop 'readhorizons: index mismatch'
          endif
          ! store only what is needed
          if (.not.crop .or. &
               & (i>=Mx1 .and. i<=Mx2 .and. j>=My1 .and. j<=My2)) then
             s(i,j,1:naz) = sread(:)
             s(i,j,naz+1) = s(i,j,1) ! wrap around
          endif
       enddo
    enddo
    close(20)
    
    s = atan(s)  ! slope -> angle
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

    if (azSun/=azSun) then ! azSun can be NaN for zenith
       getonehorizon = 0.
       return
    endif
    
    daz = 2*pi/real(naz)
    k = floor(modulo(azSun,2*pi)/daz)
    a = modulo(azSun,2*pi)/daz-k
    !if (k<0 .or. k>=naz) then 
       !print *,'azSun=',azSun,'k=',k   ! impure
       !error stop 'getonehorizon: impossible k value'
    !endif
    smax = s(i0,j0,k+1)*(1.-a) + s(i0,j0,k+2)*a
    getonehorizon = smax
  end function getonehorizon
  
  
  elemental function getoneskysize(i,j)
    ! calculates sky size (steradian) from horizons
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
    ! calculates sky size (steradian) from horizons
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
    ! calculate land view factor from horizons
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
       ! daz should be zero when line of sight coincides with direction of
       !                                                  steepest descent
       emin(k) = atan( -tan(alpha)*cos(daz(k)) )
       ! self-shadowing higher than distant horizon
       if (emin(k)>ssub(k)) ssub(k)=emin(k)
    end do
    G1 = sum(sin(ssub)**2) - sum(sin(emin)**2)
    G2 = sum(cos(daz)*( ssub+sin(ssub)*cos(ssub) - emin-sin(emin)*cos(emin) ))
    G1 = G1/naz 
    G2 = G2/naz
    
    getoneGterm = cos(alpha)*G1 + sin(alpha)*G2
    ! getoneGterm should never be negative
  end function getoneGterm

END MODULE newhorizons



pure subroutine difftopo(NSx,NSy,h,dx,dy,SlopeAngle,azFac)
  ! calculate slopes and azimuths of surface elements
  ! azFac= 0 sloped toward south
  ! azFac=+/-pi sloped toward north
  ! azFac=+pi/2 sloped toward west
  ! azFac=-pi/2 sloped toward east
  implicit none
  integer, intent(IN) :: NSx, NSy
  real(8), intent(IN) :: h(NSx,NSy), dx, dy
  real(8), intent(OUT), dimension(NSx,NSy) :: SlopeAngle, azFac
  integer i,j
  real(8) sx,sy

  do i=2,NSx-1
     do j=2,NSy-1
        sx = (h(i+1,j)-h(i-1,j))/(2.*dx)
        sy = (h(i,j+1)-h(i,j-1))/(2.*dy)
        SlopeAngle(i,j) = atan(sqrt(sx**2+sy**2))
        azFac(i,j) = atan2(sx,-sy)  ! north is up, clockwise
     enddo
  enddo
end subroutine difftopo



pure subroutine difftopo2(h,SlopeAngle,azFac,Mx1,Mx2,My1,My2)
  ! calculate slopes and azimuths of surface elements
  ! like difftopo but with different input arguments (for cropped domains)
  use filemanager, only : NSx,NSy,dx,dy
  implicit none
  integer, intent(IN) :: Mx1, Mx2, My1, My2
  real(8), intent(IN) :: h(NSx,NSy)
  real(8), intent(OUT), dimension(Mx1:Mx2,My1:My2) :: SlopeAngle, azFac
  integer i, j
  real(8) sx, sy

  do i=max(2,Mx1),min(NSx-1,Mx2)
     do j=max(2,My1),min(NSy-1,My2)
        sx = (h(i+1,j)-h(i-1,j))/(2.*dx)
        sy = (h(i,j+1)-h(i,j-1))/(2.*dy)
        SlopeAngle(i,j) = atan(sqrt(sx**2+sy**2))
        azFac(i,j) = atan2(sx,-sy)  ! north is up, clockwise
     enddo
  enddo
end subroutine difftopo2



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
  real(8), parameter :: pi=3.1415926535897932
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



subroutine getfieldofview(NSx,NSy,ffn,cc,ia,ja,dOh,landsize,CCMAX)
  ! reads subtended spherical angles from file
  implicit none
  integer, intent(IN) :: NSx, NSy
  character(len=*), intent(IN) :: ffn
  integer, intent(IN) :: CCMAX
  integer, intent(OUT) :: cc(NSx,NSy) ! number of cells in field of view
  integer(2), intent(OUT), dimension(NSx,NSy,CCMAX) :: ia, ja
  real(4), intent(OUT), dimension(NSx,NSy,CCMAX) :: dOh
  real(8), intent(OUT) :: landsize(NSx,NSy)
  integer i, j, k, i0_2, j0_2, ierr

  open(unit=20,file=ffn,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getfieldofview: input file not found'
  do i=2,NSx-1
     do j=2,NSy-1
        read(20,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') & ! format must match
             & i0_2,j0_2,cc(i,j),landsize(i,j) 
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



subroutine getviewfactors(NSx,NSy,vfn,cc,ia,ja,VF,viewsize,CCMAX)
  ! reads viewfactors from file
  implicit none
  integer, intent(IN) :: NSx, NSy
  character(len=*), intent(IN) :: vfn
  integer, intent(IN) :: CCMAX
  integer, intent(OUT) :: cc(NSx,NSy) ! number of cells in field of view
  integer(2), intent(OUT), dimension(NSx,NSy,CCMAX) :: ia, ja
  real(4), intent(OUT), dimension(NSx,NSy,CCMAX) :: VF
  real(8), intent(OUT) :: viewsize(NSx,NSy)
  integer i, j, k, i0_2, j0_2, ierr
  real(8) landsize(NSx,NSy)
  
  open(unit=21,file=vfn,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getviewfactors: input file not found'
  do i=2,NSx-1
     do j=2,NSy-1
        read(21,'(2(i5,1x),i6,1x,2(f7.5,1x))',advance='no') & !format must match
             & i0_2,j0_2,cc(i,j),landsize(i,j),viewsize(i,j)
        if (i/=i0_2 .or. j/=j0_2) stop 'getviewfactors: wrong data order'
        if (cc(i,j)>CCMAX) stop 'getviewfactors: not enough memory allocated'
        do k=1,cc(i,j)
           read(21,'(2(i5,1x),g10.4,1x)',advance='no') & ! format must match
                & ia(i,j,k),ja(i,j,k),VF(i,j,k)
        enddo
        read(21,'()')
     enddo
  enddo
  close(21)
end subroutine getviewfactors



subroutine getviewfactors_full(NSx,NSy,vfn,viewsize,VF)
  ! reads full viewfactor matrix from file
  implicit none
  integer, intent(IN) :: NSx, NSy
  character(len=*), intent(IN) :: vfn
  real(8), intent(OUT) :: viewsize(NSx,NSy)
  real(4), intent(OUT) :: VF(NSx,NSy,(NSx-2)*(NSy-2))
  integer cc, i, j, i0_2, j0_2, ierr
  real(8) landsize(NSx,NSy)
  
  open(unit=21,file=vfn,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getviewfactors_full: input file not found'
  do i=2,NSx-1
     do j=2,NSy-1
        read(21,*) i0_2,j0_2,cc,landsize(i,j),viewsize(i,j),VF(i,j,:)
        if (i/=i0_2 .or. j/=j0_2) stop 'getviewfactors_full: wrong data order'
     enddo
  enddo
  close(21)
end subroutine getviewfactors_full



integer function getmaxfieldsize(NSx,NSy,ffn)
  ! compatible with fieldofviews.dat and viewfactors.dat
  implicit none
  integer, intent(IN) :: NSx,NSy
  character(len=*), intent(IN) :: ffn
  integer maxsize
  integer cc, i, j, i0_2, j0_2, ierr

  print *,'reading ',ffn
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



subroutine readtruncSVD(U,w,V,Nflat,T,pth,ext)
  ! read truncated SVD from file
  implicit none
  integer, intent(IN) :: Nflat, T
  real(8), intent(OUT) :: U(Nflat,T), w(Nflat), V(Nflat,T)
  character(len=*), intent(IN) :: pth, ext
  integer k, ierr

  if (T<1 .or. T>Nflat) stop 'impossible value for T'

  open(20,file=pth//'svd_U'//ext,action='read',iostat=ierr)
  if (ierr>0) then
     print *, 'readtruncSVD: input file for U not found','svd_U'//ext
     stop
  endif
  do k=1,Nflat
     read(20,*) U(k,1:T) ! matrix U
  end do
  close(20)
  
  open(21,file=pth//'svd_sigma'//ext,action='read',iostat=ierr)
  if (ierr>0) stop 'readtruncSVD: input file for w not found'
  read(21,*) w(1:T)  ! diagonal of matrix W
  close(21)
  
  open(22,file=pth//'svd_V'//ext,action='read',iostat=ierr)
  if (ierr>0) stop 'readtruncSVD: input file for V not found'
  do k=1,T
     read(22,*) V(1:Nflat,k)  ! matrix V-transpose
  end do
  close(22)
end subroutine readtruncSVD



integer function countcolumns()
  ! counts the number of azimuth rays in horizons file
  ! incorporates code from Léo https://stackoverflow.com/users/10478255/leo
  use filemanager, only : sfn
  implicit none
 
  integer i, io, ierr
  integer, parameter :: MAX_NUM_OF_COLS=999
  integer, parameter :: MAX_LINE_LENGTH=10000
  character(len=MAX_LINE_LENGTH) line
  real(8), dimension(MAX_NUM_OF_COLS) :: test_array
       
  open(unit=20,file=sfn,status='old',action='read',iostat=ierr)
  if (ierr>0) then
     print *,sfn
     stop 'countcolumns: Input file not found'
  endif

  ! Get first line of file.
  DO
    READ(20,'(A)',iostat=io) line
    IF (io/=0) stop "Error reading file."
    exit 
  ENDDO

  CLOSE(20)

  do i=1,MAX_NUM_OF_COLS
    READ(line,*,iostat=io) test_array(1:i)
    if (io==-1) exit
  enddo

  !write(*,*) 'number of columns = ', (i-1) 
  countcolumns = i-1
  
end function countcolumns



subroutine update_terrain_irradiance(VF,cc,ii,jj,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
  ! multiply SW and LW irradiances with viewfactors
  use filemanager, only : NSx, NSy
  implicit none
  real(4), intent(IN) :: VF(:,:,:)
  integer, intent(IN) :: cc(NSx,NSy)
  integer(2), intent(IN) :: ii(:,:,:), jj(:,:,:)
  real(8), intent(IN), dimension(NSx,NSy) :: Qn, albedo, Tsurf
  real(8), intent(IN) :: emiss
  real(8), intent(INOUT), dimension(NSx,NSy) :: Qrefl, QIRin
  real(8), parameter :: sigSB = 5.6704e-8  
  integer i, j, iii, jjj, k
  real(8), dimension(NSx,NSy) :: Qvis, QIRout

  Qvis(:,:) = albedo * (Qn + Qrefl)
  QIRout(:,:) = emiss*sigSB*Tsurf**4 + (1-emiss)*QIRin
  
  do i=2,NSx-1
     do j=2,NSy-1
        Qrefl(i,j)=0.
        QIRin(i,j)=0.
        do k=1,cc(i,j)
           iii = ii(i,j,k); jjj = jj(i,j,k)
           Qrefl(i,j) = Qrefl(i,j) + VF(i,j,k)*Qvis(iii,jjj)
           QIRin(i,j) = QIRin(i,j) + VF(i,j,k)*QIRout(iii,jjj)
        enddo
     enddo
  enddo
end subroutine update_terrain_irradiance



subroutine update_terrain_irradiance_full(VF,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
  ! multiply SW and LW irradiances using full (Nflat x Nflat) viewfactor matrix
  use filemanager, only : NSx, NSy
  implicit none
  integer, parameter :: Nflat = (NSx-2)*(NSy-2)
  real(4), intent(IN) :: VF(NSx,NSy,Nflat)
  real(8), intent(IN), dimension(NSx,NSy) :: Qn, albedo, Tsurf
  real(8), intent(IN) :: emiss
  real(8), intent(INOUT), dimension(NSx,NSy) :: Qrefl, QIRin
  real(8), parameter :: sigSB = 5.6704e-8  
  integer i, j, ii, jj, k
  real(8), dimension(NSx,NSy) :: Qvis, QIRout

  Qvis(:,:) = albedo * (Qn + Qrefl)
  QIRout(:,:) = emiss*sigSB*Tsurf**4 + (1-emiss)*QIRin

  do i=2,NSx-1
     do j=2,NSy-1
        Qrefl(i,j)=0.
        QIRin(i,j)=0.
        do k=1,Nflat
           call k2ij(k,NSy,ii,jj)
           Qrefl(i,j) = Qrefl(i,j) + VF(i,j,k)*Qvis(ii,jj)
           QIRin(i,j) = QIRin(i,j) + VF(i,j,k)*QIRout(ii,jj)
        end do
     end do
  end do
end subroutine update_terrain_irradiance_full



subroutine update_terrain_irradiance_SVD(T,U,w,V,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
  ! multiply SW and LW irradiances using truncated SVD
  use filemanager, only : NSx, NSy
  implicit none
  integer, parameter :: Nflat = (NSx-2)*(NSy-2)
  integer, intent(IN) :: T  ! number of modes retained
  real(8), intent(IN) :: U(:,:), w(T), V(:,:)
  real(8), intent(IN), dimension(NSx,NSy) :: Qn, albedo, Tsurf
  real(8), intent(IN) :: emiss
  real(8), intent(INOUT), dimension(NSx,NSy) :: Qrefl, QIRin
  real(8), parameter :: sigSB = 5.6704e-8  
  integer iii, jjj, k, l, i, j
  real(8), dimension(NSx,NSy) :: Qvis, QIRout
  real(8) tmpv1(T), tmpv2(T)
  integer, external :: ij2k

  Qvis(:,:) = albedo * (Qn + Qrefl)
  QIRout(:,:) = emiss*sigSB*Tsurf**4 + (1-emiss)*QIRin
  
  do k=1,T ! (V^Transpose)*RHS
     tmpv1(k) = 0.
     tmpv2(k) = 0.
     do iii=2,NSx-1
        do jjj=2,NSy-1
           l = ij2k(iii,jjj,NSy)
           tmpv1(k) = tmpv1(k) + v(l,k)*Qvis(iii,jjj)
           tmpv2(k) = tmpv2(k) + v(l,k)*QIRout(iii,jjj)
        enddo
     enddo
  enddo
  do k=1,T ! w*(V^Transpose)*RHS
     tmpv1(k) = w(k)*tmpv1(k)
     tmpv2(k) = w(k)*tmpv2(k)
  enddo
  Qrefl(:,:)=0.; QIRin(:,:)=0.
  do k=1,Nflat ! u*w*(V^Transpose)*RHS
     call k2ij(k,NSy,i,j)
     Qrefl(i,j) = Qrefl(i,j) + sum( u(k,1:T) * tmpv1(1:T) )
     QIRin(i,j) = QIRin(i,j) + sum( u(k,1:T) * tmpv2(1:T) )
  end do

end subroutine update_terrain_irradiance_SVD
