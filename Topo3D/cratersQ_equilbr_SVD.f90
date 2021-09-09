program cratersQ_equilbr_SVD
!***********************************************************************
!   cratersQ_equilbr_SVD: program to calculate surface energy balance
!             for constant sun position; equilibrium solution
!             uses truncated SVD of view factor matrix  
!***********************************************************************
  use filemanager
  use allinterfaces
  use newhorizons
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8  

  integer n, i, j
  integer, parameter :: Nflat = (NSx-2)*(NSy-2)
  integer, parameter :: T=200  ! truncate after T terms
  real(8), dimension(NSx,NSy) :: h, SlopeAngle, azFac
  real(8), dimension(NSx,NSy) :: Qn, QIRin, Qrefl  ! incoming
  real(8) R, azSun, smax, emiss, betaSun, Qshadow1, Qshadow2
  real(8), dimension(NSx,NSy) :: Tsurf, Qabs, albedo
  real(8) w(Nflat)
  real(8), allocatable :: u(:,:), v(:,:)
  integer, external :: ij2k

  ! azimuth in degrees east of north, 0=north facing
  albedo(:,:) = 0.12
  emiss = 0.95

  R=1.
  azSun = 180.*d2r
  betaSun = 10.*d2r

  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,SlopeAngle,azFac)

  print *,'...reading horizons file...'
  call readhorizons

  ! read outputs of xsvdcmp
  print *,'...reading SVD files...'
  allocate(u(Nflat,T), v(Nflat,T))
  call readtruncSVD(U,w,V,Nflat,T,'./','.dat')

!!$  test: block
!!$    integer ierr, nelem
!!$    real(8) buf, VF(Nflat,Nflat)
!!$    print *,'...testing SVD decomposition...'
!!$    open(7,file='viewfactors.dat',status='old',action='read',iostat=ierr)
!!$    if (ierr>0) stop 'input file not found'
!!$    do k=1,Nflat
!!$       read(7,*) i,j,nelem,buf,buf,(VF(k,l), l=1,Nflat)
!!$       if (nelem/=Nflat) stop 'wrong input'
!!$    enddo
!!$    close(7)
!!$    ! Product U*W*(V^Transpose):
!!$    do k = 1,Nflat
!!$       do l = 1,Nflat
!!$          b = sum( u(k,:) * w(:) * v(l,:) )  ! = b(k,l)
!!$          if ( abs(a(k,l)-b)>1d-6 ) print *,'A-SVD(A)',a(k,l),b
!!$       enddo
!!$    enddo
!!$  end block test
  
  print *,'...calculating...'  
  print *,'beta=',betaSun/d2r,'azSun=',azSun/d2r

  do i=2,NSx-1
     do j=2,NSy-1
        smax = getonehorizon(i,j,azSun)
        Qn(i,j) = flux_wshad(R,sin(betaSun),azSun,SlopeAngle(i,j),azFac(i,j),smax)
     enddo
  enddo

  Tsurf(:,:) = ((1-albedo)*Qn/sigSB)**0.25
  Qrefl=0.; QIRin=0.

  do n=1,7
     !print *,'Iteration',n

     call update_terrain_irradiance_SVD(T,U,w(1:T),V,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
     !call update_terrain_irradiance_full(VF,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
       
     Qshadow1 = 0.; Qshadow2 = 0.
     do i=2,NSx-1
        Qshadow1 = Qshadow1 + sum(Qrefl(i,:),1,Qn(i,:)==0)
        Qshadow2 = Qshadow2 + sum(QIRin(i,:),1,Qn(i,:)==0)
     enddo
     
     ! Qabs is Q absorbed
     Qabs(:,:) = (1.-albedo(:,:))*(Qn+Qrefl) + emiss*QIRin
    
     Tsurf(:,:) = (Qabs/sigSB/emiss)**0.25

     print *,'Iteration',n,Qshadow1,Qshadow2
  enddo

  deallocate(u, v)

  open(unit=21,file='qinst.dat',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,4(1x,f6.1),1x,f5.1)') &
             & i,j,h(i,j),SlopeAngle(i,j), &
             & Qn(i,j),Qirin(i,j),Qrefl(i,j),Qabs(i,j),Tsurf(i,j)
     enddo
  enddo
  close(21)

end program cratersQ_equilbr_SVD

