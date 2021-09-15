PROGRAM fieldproperties
!***********************************************************************
! read horizons file and calculate skysizes from it
! optionally read file with facet areas
! optionally read file with viewfactors
! used for diagnostics only
!***********************************************************************
  use filemanager
  use allinterfaces
  use newhorizons

  implicit none
  integer i, j, cc(NSx,NSy)
  real(8), dimension(NSx,NSy) :: h, skysize1, skysize2, viewsize
  ! skysize + landsize = 2*pi
  ! skyviewfactor + landviewfactor = 1
  real(8), dimension(NSx,NSy) :: gterm, SlopeAngle, azFac
  real(8), parameter :: pi=3.1415926535897932

  logical, parameter :: fieldofview=.false.
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,SlopeAngle,azFac)
  
  print *,'...reading horizons file...'
  call readhorizons
  do i=2,NSx-1
     do j=2,NSy-1
        skysize1(i,j) = getoneskysize(i,j)
        skysize2(i,j) = getoneskysize_v2(i,j)
        gterm(i,j) = getoneGterm(i,j,SlopeAngle(i,j),azFac(i,j))
     enddo
  enddo
  
  if (fieldofview) then
     
     !call readviewfactors_mini(NSx,NSy,ffn,cc,landsize)
     call readviewfactors_mini(NSx,NSy,vfn,cc,viewsize)
     
     do i=2,NSx-1
        do j=2,NSy-1
           !call check_symmetry(i,j,h,VF)
        end do
     end do
  end if

  print *,'...writing output...'
  open(unit=21,file='tmp.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        !write(21,'(2(i4,1x),f9.2,3(1x,f6.3))') i,j,h(i,j), &
        !     & 2*pi-skysize1(i,j),2*pi-skysize2(i,j),landsize(i,j)
        write(21,'(2(i4,1x),f9.2,2x,3(2x,f5.3))') &
             & i,j,h(i,j),2*pi-skysize2(i,j),gterm(i,j) ! ,viewsize
     enddo
  enddo
  close(21)

END PROGRAM fieldproperties



subroutine readviewfactors_mini(NSx,NSy,vfn,cc,viewsize)
  ! reads viewfactors from file
  implicit none
  integer, intent(IN) :: NSx, NSy
  character(len=*), intent(IN) :: vfn
  integer, intent(OUT) :: cc(NSx,NSy) ! number of facets with direct line of sight
  real(8), intent(OUT) :: viewsize(NSx,NSy)
  integer i, j, i0_2, j0_2, ierr
  
  open(unit=21,file=vfn,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'getviewfactors: input file not found'
  do i=2,NSx-1
     do j=2,NSy-1
        read(21,'(2(i5,1x),i6,1x,f7.5,1x)',advance='no') & ! format must match
             & i0_2,j0_2,cc(i,j),viewsize(i,j)
        if (i/=i0_2 .or. j/=j0_2) stop 'readviewfactors_mini: wrong data order'
        read(21,'()')
     enddo
  enddo
  close(21)
end subroutine readviewfactors_mini
