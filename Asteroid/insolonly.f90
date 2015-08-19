subroutine insolonly(latitude,a,omega,ecc,eps,Trot,Q0mean,Qmean,Q4mean)
  ! insolation only
  ! structure parallels oneasteroid_thermal1d.f90
  implicit none
  real(8), intent(IN) :: latitude, a, omega, eps, ecc, Trot
  real(8), intent(OUT) :: Qmean, Q0mean, Q4mean
  real(8), parameter :: pi=3.1415926535897932
  integer n, nr
  real(8) edays, Rau, Ls, decl, Torb
  real(8) Q, HA, Q0
  real(8), external :: flux_noatm, a2Torb, flux2T

  Torb = a2Torb(a)
  Q0mean=0; Qmean=0.; Q4mean=0.
  nr = 100*nint(Torb/Trot)
  print *,'# steps=',nr
  do n=1,nr
     edays = n*Torb/real(nr)
     call generalorbit(edays,a,ecc,omega,eps,Ls,decl,Rau) 
     Q0 = flux_noatm(Rau,decl,latitude,0d0,0d0,0d0)
     HA = mod(edays/Trot,1.d0)*2.*pi
     Q = flux_noatm(Rau,decl,latitude,HA,0d0,0d0)
     if (Q>0.) then
        Q0mean = Q0mean + Q0
        Qmean = Qmean + Q   ! thIn = infinity
        Q4mean = Q4mean + Q**0.25  ! thIn = 0
     endif
  enddo
  Q0mean = Q0mean/nr
  Qmean = Qmean/nr
  Q4mean = (Q4mean/nr)**4

  ! scale Qmean with (1-albedo)
  ! scale Q4mean with (1-albedo)
end subroutine insolonly

