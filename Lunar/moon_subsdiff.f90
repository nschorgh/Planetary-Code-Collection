program moon_subsdiff
  ! microphysical model of subsurface diffusion on airless body
  ! (improved version of what I wrote in 2006)
  implicit none
  integer :: nrmax= 1000  ! initial maximum number of molecules
  integer :: idum = -299  ! seed for random number generator
  integer :: MAXDEPTHIDX 
  integer, parameter :: nadd = 1   ! molecules added on surface
  integer, parameter :: GONE = -9999
  integer, parameter :: STEPSPERSOL = 96
  integer, parameter :: vMONO = STEPSPERSOL*4*nadd  ! =1/4 monolayer/lunation
  real(8), parameter :: secyear = 365.24*86400, lunation = 29.52*86400
  real(8) maxtime ! [s]
  real(8) dz      ! mean-free path [m]
  real(8), parameter :: dt = lunation/STEPSPERSOL ! time step [s]
  parameter (MAXDEPTHIDX = 1000, dz = 0.1e-3, maxtime = 1000.*secyear)
  !character(len=*), parameter :: fn = 'Tprofiles_finer_moonranger_100um.dat'
  logical, allocatable :: mobile(:), mobileswap(:)
  integer :: j, k, jj, nr, avcc=0, nrsource=0, avcc2=0, jold=0, zi
  integer, allocatable :: zidx(:), zidxswap(:)
  integer, dimension(0:MAXDEPTHIDX) :: binc, avbinc, avbinc2
  real(8), dimension(0:MAXDEPTHIDX) :: vvm, sigma 
  real(8) time, T(0:MAXDEPTHIDX,STEPSPERSOL)
  
  allocate(zidx(nrmax), source=GONE)
  allocate(mobile(nrmax), source=.TRUE.)

  print *,'Surface molecules added at each time step',nadd
  print *,'Monolayer',vMONO
  print *,'Time step = ',dt,'sec'
  print *,'Maximum run time = ',maxtime/secyear,'years'
  write(*,'(a,1x,f6.3)') ' Domain depth =',MAXDEPTHIDX*dz
  print *,'Mean free path and grain size = ',dz
  print *,'Number of vertical grid levels = ',MAXDEPTHIDX
  open(23,file='profiles.dat',action='write')
  open(24,file='longseries.dat',action='write')

  ! read all temperature profiles
  call temperatureprofile(MAXDEPTHIDX,STEPSPERSOL,T,dz)
  !call temperatureprofile_fromdata(MAXDEPTHIDX,STEPSPERSOL,T,fn)

  binc(:) = 0
  avbinc(:) = 0

  do j=1,nint(maxtime/dt)  ! begin time loop
     time = j*dt
     jj = mod(j-1,STEPSPERSOL)+1  ! 1...STEPSPERSOL

     ! add molecules to surface
     if (nadd>0) then
        nr = 0
        do k=1,nrmax
           if (zidx(k)<0) then ! empty slot
              zidx(k) = 0  ! added molecule to surface
              nr = nr+1
           end if
           if (nr==nadd) exit
        enddo
        if (nr<nadd) then  ! increase array size
           
           allocate(zidxswap(nrmax))
           zidxswap(:) = zidx(:)
           deallocate(zidx)
           allocate(zidx(2*nrmax))
           zidx(1:nrmax) = zidxswap(:)
           zidx(nrmax+1:2*nrmax) = GONE
           deallocate(zidxswap)
           zidx(nrmax+1:nrmax+nadd-nr) = 0 ! add remaining molecules

           allocate(mobileswap(nrmax))
           mobileswap(:) = mobile(:)
           deallocate(mobile)
           allocate(mobile(2*nrmax))
           mobile(1:nrmax) = mobileswap(:)
           mobile(nrmax+1:2*nrmax) = .TRUE.
           deallocate(mobileswap)
           
           print *,'doubled array length from ',nrmax
           nrmax = 2*nrmax
        endif
        nrsource = nrsource + nadd        
     endif

     ! update adsorbate profile
     vvm(:) = binc(:) / real(vMONO,8)  ! number of monolayers

     ! subsurface vapor diffusion
     do k = 1,nrmax
        if (mobile(k)) then
           call randomwalk1(zidx(k), dt, idum, MAXDEPTHIDX, T(:,jj), vvm(:) )
        end if
     enddo

     ! update concentration profile and identify immobilized molecules
     binc(:) = 0
     do k=1,nrmax
        zi = zidx(k)
        if (zi<0 .or. zi>MAXDEPTHIDX) cycle
        binc(zi) = binc(zi)+1
        if (binc(zi) <= vMONO) then
           mobile(k) = .TRUE.
        else
           mobile(k) = .FALSE.
        end if
     enddo
     !if (sum(binc) /= count(zidx>=0)) stop 'discrepancy'

     block  ! average concentration profile
       real(8), parameter :: outinterv = 100.*secyear
       if (mod(time,outinterv) < lunation*1.0001 .and. avcc2<STEPSPERSOL ) then
          avbinc2 = avbinc2(:) + binc(:)
          avcc2 = avcc2 + 1
          jold = j
       endif
       if (j == jold+1) then
          if (avcc2 /= STEPSPERSOL) then
             print *,'WARNING: bad average',time/secyear,avcc2,STEPSPERSOL
          endif
          sigma = avbinc2(:) / real(avcc2,8)
          write(23,('(f9.3,*(1x,f0.3))')) (time-lunation/2.)/secyear, sigma
          avbinc2(:) = 0;  avcc2 = 0
       endif
     end block
     
     !if (mod(j,12)==0) then
     if (mod(j,137)==0) then  ! reduce volume of output
        write(24,('(f11.4,2(1x,i0))')) time/secyear,count(zidx>0),nrsource
     endif
     
     if (time>maxtime) exit

     if (time >= maxtime-lunation) then
        avbinc = avbinc(:) + binc(:)
        avcc = avcc + 1
     endif
     
  enddo

  sigma = avbinc(:) / real(avcc,8)
  if (avcc /= STEPSPERSOL) print *,'Warning: averaging is off',avcc,STEPSPERSOL
  write(23,('(f9.3,*(1x,f0.3))')) (time-lunation/2.)/secyear, sigma
  
  print *,'Total time',time,time/secyear
  print *,'Source molecules',nrsource
  close(23)
  close(24)

end program moon_subsdiff



subroutine temperatureprofile(nz,STEPSPERSOL,T,dz)
  ! assign theoretical temperature profile
  implicit none
  integer, intent(IN) :: nz, STEPSPERSOL
  real(8), intent(IN) :: dz
  real(8), intent(INOUT) :: T(0:nz,STEPSPERSOL)
  real(8), parameter :: delta = 0.05  ! thermal skin depth
  !real(8), parameter :: g=0.  ! geothermal gradient [K/m]  
  real(8), parameter :: g=18.  ! geothermal gradient [K/m]  g = 0.018/0.001
  real(8), parameter :: pi = 3.1415926535
  integer i, j
  real(8) Tm, Ta, z, phi

  Tm=113.
  Ta=40.

  do j=1,STEPSPERSOL
     ! time = j*lunation/STEPSERSOL
     ! w*time = 2*pi*j/STEPSPERSOL
     phi = 2*pi*real(j,8)/STEPSPERSOL
     do i=0,nz
        z = i*dz
        T(i,j) = Tm + Ta*exp(-z/delta)*sin(z/delta-phi) + g*z
     end do
  end do
end subroutine temperatureprofile



subroutine temperatureprofile_fromdata(nz,STEPSPERSOL,T,fn)
  ! assign temperature profiles from file
  implicit none
  integer, intent(IN) :: nz, STEPSPERSOL
  real(8), intent(OUT) :: T(0:nz,STEPSPERSOL)
  character(*), intent(IN) :: fn
  integer j, ierr
  
  T(:,:) = 0.
  
  open(41,file=fn,action='read',iostat=ierr)
  if (ierr>0) stop 'file not found'

  ! it's okay if input is longer than nz
  do j=1,STEPSPERSOL
     read(41,*) T(:,j)
  enddo
  close(41)
end subroutine temperatureprofile_fromdata

