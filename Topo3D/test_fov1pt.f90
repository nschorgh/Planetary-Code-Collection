PROGRAM test_fieldofviews1pt
!************************************************************************
! 1. calculates horizons for one point
! 2. determines visibility from one point
! 3. calculates view factors
!************************************************************************
  use allinterfaces
  use findvisibletopo, naz_s => naz, azRay_s => azRay
  implicit none
  real(8), parameter :: pi=3.1415926535897932, r2d=180./pi
  integer i0, j0, ii, jj, k
  real(8) h(NSx,NSy), smax(naz_s)
  logical visibility(NSx,NSy)
    
  call readdem(h)
  print *,'...finished reading topography...'
  print *,'# Nsx=',NSx,'Nsy=',NSy

  i0=14; j0=29
  !i0=41; j0=41
  print *,i0,j0
  print *

  print *,'# azimuths=',naz_s,'(single grid)'
  call findallhorizon_wsort_v3(h,i0,j0,smax,visibility)
  do k=1,naz_s
     write(20,'(f6.2,1x,f6.4)') azRay_s(k)*r2d,smax(k)
  enddo
  call refinevisibility_cart(i0,j0,h,visibility)
  open(unit=25,file='singlegrid.dat',action='write')
  do ii=1,NSx
     do jj=1,NSy
        if (visibility(ii,jj)) write(25,'(2(2x,i0))') ii,jj
     enddo
     !write(26,*) visibility(ii,:)
  enddo
  close(25)
  open(unit=23,file='viewfactors.dat',action='write')
  call findviewfactors(h,i0,j0,23,visibility)
  close(23)
  
END PROGRAM test_fieldofviews1pt

