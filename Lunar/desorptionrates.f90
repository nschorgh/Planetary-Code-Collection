module desorp
  implicit none
  real(8), parameter :: kBeV=8.617333262e-5 ! [eV/K]
  !real(8), parameter :: nu=1e16
  real(8), parameter :: thetam=1e19
  real(8), parameter :: Eice=0.529, Ec=0.65, W=0.22
  !real(8), parameter :: Eice=0.65, Ec=0.65, W=0.0
end module desorp


subroutine sitedistribution(NS,maxE,nArea)
  use desorp, only: Ec, W
  implicit none
  integer, intent(IN) :: NS
  real(8), intent(IN) :: maxE
  real(8), intent(OUT) :: nArea(NS)
  integer k
  real(8) dE, E(NS)
  
  dE = maxE/real(NS)
  do concurrent (k=1:NS)
     E(k) = k*dE
     if (E(k)<Ec) then
        nArea(k) = 0.
     else
        nArea(k) = 1./W * exp(-(E(k)-Ec)/W)
        !print *,k,n(k),(ThetaS(k)<=n(k)*thetam)
     end if
  enddo
  !write(*,'(a,1x,f0.3)') 'sum(n)=',sum(n)*dE
end subroutine sitedistribution


elemental function desorptionrate(T,theta)
  ! Desorption rate of water molecules from lunar grains as a function of
  ! temperature and coverage according to Schorghofer, PSJ (2023), PSJ 4, 164.
  ! assumes efficient grain surface diffusion
  use desorp
  implicit none
  real(8), intent(IN) :: T ! [K]
  real(8), intent(IN) :: theta ! [molecules/m^2]
  real(8) desorptionrate  ! [molecules/m^2/s]
  real(8) nu, v, S, gamma, S1, S2, b
  
  v = theta/thetam
  nu = 2.2e17/sqrt(T)
  b = 1./(kBeV*T)
  if (v<=0.) then
     S = 0.
  elseif (v<=1) then  ! sub-monolayer
     S = nu * theta * exp(-Ec*b) *v**(W*b) / (1.+W*b)
  else ! multi-layer
     ! option 1 (expansion of BET formula around infinite v)
     !S = nu * thetam * ( exp(-Eice*b)*(1-1./v) + exp(-Ec*b) / (1.+W*b) /v**2 )

     ! option 2 (based on rescaled distribution of desorption energies)
     gamma = 1. / ( 1 + (v-1)*(Ec-Eice+W)*b )
     S2 = nu * thetam * exp(-Eice*b) * exp(-gamma*(Ec-Eice)*b) / (1+gamma*W*b)

     ! option 3 (composite)
     gamma = 1./v**2
     S1 = nu*thetam * exp(-Eice*b) * &
          & exp(-gamma*(Ec-Eice)*b) / (1 + gamma* W*b)
     !S = min(S1,S2)
     S = S1
  end if

  desorptionrate = S
  !if (S/=S) print *,T,theta,'S undefined'
end function desorptionrate


elemental function desorptionrate_ice(T)
  ! sublimation rate of ice
  use desorp
  implicit none
  real(8), intent(IN) :: T   ! [K]
  real(8) desorptionrate_ice ! [molecules/m^2/s]
  real(8) Sice, nu

  nu = 2.2e17/sqrt(T)
  Sice = thetam*nu*exp(-Eice/kBeV/T)
  desorptionrate_ice = Sice
end function desorptionrate_ice


elemental subroutine desorptionrateXi1(Ek,T,ThetaS,S,dSdTheta)
  ! spectral version of function desorptionrate
  ! only for sub-monolayer coverage
  use desorp, only: kBeV
  implicit none
  real(8), intent(IN) :: Ek     ! desorption energy [eV]
  real(8), intent(IN) :: T      ! [K]
  real(8), intent(IN) :: ThetaS ! [molecules/m^2/eV]
  real(8), intent(OUT) :: S, dSdTheta
  real(8) nu
  nu = 2.2e17/sqrt(T)
  dSdTheta = nu * exp(-Ek/(kBeV*T))
  S = ThetaS * dSdTheta
end subroutine desorptionrateXi1


elemental function desorptionrateXi1_multi(T,Ek,nArea,ThetaS)
  use desorp, only: thetam, kBeV, Eice
  implicit none
  real(8), intent(IN) :: T ! [K]  
  real(8), intent(IN) :: Ek  ! desorption energy [eV]
  real(8), intent(IN) :: nArea
  real(8), intent(IN) :: ThetaS ! [molecules/m^2/eV]
  real(8) Xi, desorptionrateXi1_multi  ! [molecules/m^2/s/eV]
  real(8) nu, nmono, gamma, b, v
  
  if (nArea==0.) then
     Xi=0.
  else
     nu = 2.2e17/sqrt(T)     
     nmono = nArea*thetam
     b = 1./(kBeV*T)
     if (ThetaS<=nmono) then
        Xi = nu * ThetaS * exp(-Ek*b)
     else
        v = ThetaS/nmono
        
        gamma = 1./v**2
        Xi = nu * nmono * exp(-Eice*b) * exp(-gamma*(Ek-Eice)*b)

        !gamma = 1. / ( 1 + (v-1)*(Ek-Eice+W)*b )
        !Xi2 = nu * nmono * exp(-Eice*b) * exp(-gamma*(Ek-Eice)*b)

        !Xi = min(Xi1,Xi2)
     end if
     !write(*,'(f0.3,1x,g0.3)') ThetaS/(nArea*thetam),S
  end if
  desorptionrateXi1_multi = Xi
end function desorptionrateXi1_multi


subroutine desorptionrate_inverse_Xi(T,S,NS,maxE,ThetaS)
  ! inverse of function desorptionrateXi_multi (end-member 2)
  use desorp
  implicit none
  integer, intent(IN) :: NS
  real(8), intent(IN) :: T, S, maxE
  real(8), intent(OUT) :: ThetaS(NS)
  integer k
  real(8) nu, Eprime, Ek, dE, nmono, nArea(NS)

  call sitedistribution(NS,maxE,nArea)
  dE = maxE/real(NS)
  nu = 2.2e17/sqrt(T)
  do k=1,NS
     nmono = nArea(k)*thetam
     Ek = k*dE
     if (Ek<=kBeV*T*log(nu*thetam/S)) then
        ThetaS(k) = S/nu*nArea(k)*exp(Ek/kBeV/T)
     else
        Eprime = -kBeV*T*log(S/nu/thetam)
        ThetaS(k) = nmono * sqrt( (Ek-Eice) / (Eprime-Eice) )
     endif
  end do
end subroutine desorptionrate_inverse_Xi


subroutine distribute_pseudo(NS,nArea,maxE,ThetaS)
  ! redistribute adsorbate mass into highest energies to mimic endmember with
  ! maximum surface diffusion
  use desorp, only: thetam
  implicit none
  integer, intent(IN) :: NS
  real(8), intent(IN) :: nArea(NS)
  real(8), intent(IN) :: maxE
  real(8), intent(INOUT) :: ThetaS(NS)
  integer k
  real(8) theta, dE

  dE = maxE/NS
  theta = sum(ThetaS)*dE ! total adsorbate mass
  !before = theta
  ThetaS(:) = 0.
  k = NS
  do while(theta>0. .and. k>0)
     if ( theta > nArea(k)*thetam*dE ) then
        ThetaS(k) = nArea(k)*thetam
        theta = theta - nArea(k)*thetam*dE
     else
        ThetaS(k) = theta/dE
        theta = 0.
     end if
     k=k-1
  end do
  !write(*,'(a,2(1x,e0.4))') '#',before,sum(ThetaS)*dE
  !if ((k+1)*dE<Ec) error stop ! more than monolayer
end subroutine distribute_pseudo
