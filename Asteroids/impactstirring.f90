real(8) function pareto(idum,mean)
  ! generate Pareto distribution 
  ! f(x) = a/(x+1)^(a+1)
  ! \int f = 1, \int x f = 1/(a-1), \int x^2 f = 2/(a^2-3a+2)
  implicit none
  integer, intent(INOUT) :: idum
  real(8), intent(IN) :: mean
  real(8), parameter :: a=2.
  real(8) x
  real(8), external :: ran2

  x=ran2(idum)
  x=1./x**(1./a) - 1.
  x=x*(a-1.)*mean  ! changes 1st moment to mean
  pareto = x
end function pareto


subroutine impactstirring(nz,z,dt,sigma)
  ! statistical model of impact stirring/gardening
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), dt  ! dt in years
  real(8), intent(INOUT) :: sigma(nz)
  real(8), parameter :: Dgarden = 5e-11 ! gardening parameter (m^2/yr)
  integer j,k,NR,turns(nz),NT,i
  real(8) x,rho(nz),eavrho(nz),meanz,rhomean,zeff
  real(8), external :: pareto, colint
  integer, save :: idum = -948

  if (maxval(sigma)==0.) return  ! no ice
  if (dt<1.) return  ! almost nothing is going to happen in 1 year

  NR=1000
  eavrho = 0.
  turns = 0

  ! Set NT and meanz
  ! NT*meanz^2 should be proportional to dt
  ! this proportionality constant, Dgarden, is set by impactor flux
  ! there should be many grid points < meanz
  NT=1
  meanz=sqrt(dt*Dgarden/4.)   ! factor 4 makes 1st moment grow as sqrt(2*D*t)
  do while (NT<100 .and. meanz>z(15))  ! try to increase NT
     NT = 4*NT
     meanz = meanz/2.
  end do
  print *,'# impact stirring: mean=',meanz,'NT=',NT
  if (meanz<z(5)) print *,'Warning: Stirring depth small'

  do j=1,NR  ! realizations
     rho = sigma

     do i=1,NT
        x=pareto(idum,meanz)
        do k=1,nz
           if (z(k)>x) exit
        enddo
        k = k-1
        ! 0<k<=nz
        if (k>0) then
           zeff=colint(spread(real(1.,8),1,nz),z,nz,1,k)
           rhomean=colint(rho,z,nz,1,k)/zeff

           ! test mass conservation
           !buf=colint(spread(rhomean,1,nz),z,nz,1,k)
           !print *,colint(rho,z,nz,1,k),buf)
           rho(1:k)=rhomean
           turns(k) = turns(k)+1
        endif
     enddo
     eavrho(:) = eavrho(:) + rho(:)
  enddo

  sigma = eavrho/NR

  ! optional diagnostics
  !write(30,'(999(f8.3,1x))') sigma(:)
  !write(32,*) turns
end subroutine impactstirring
