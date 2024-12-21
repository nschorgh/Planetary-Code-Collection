module desorp
  implicit none
  real(8), parameter :: kBeV=8.617333262e-5 ! [eV/K]
  !real(8), parameter :: nu=1e16
  real(8), parameter :: thetam=1e19
  real(8), parameter :: Eice=0.529, Ec=0.65, W=0.22
end module desorp


elemental function desorptionrate(T,theta)
  ! Desorption rate of water molecules from lunar grains as a function of
  ! temperature and coverage according to Schorghofer, PSJ (2023), PSJ 4, 164.
  ! assumes efficient grain surface diffusion
  use desorp
  implicit none
  real(8), intent(IN) :: T ! [K]
  real(8), intent(IN) :: theta ! [molecules/m^2]
  real(8) desorptionrate  ! [molecules/m^2/s]
  real(8) nu, v, S, gamma, S1, S2,b
  
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
     gamma = 1. / ( 1 + (v-1)*(Ec-Eice+W)/(kBeV*T) )
     S2 = nu * thetam * exp(-Eice*b) * exp(-gamma*(Ec-Eice)*b) / (1+gamma*W*b)

     ! option 3 (composite)
     gamma = 1./v**2
     S1 = nu*thetam * exp(-Eice*b) * &
          & exp(-gamma*(Ec-Eice)*b) / (1 + gamma* W*b)
     S = min(S1,S2)
  end if

  desorptionrate = S
  !if (S/=S) print *,T,theta,'S undefined'
end function desorptionrate


elemental function desorptionrate_ice(T)
  ! sublimation rate of ice
  use desorp
  implicit none
  real(8), intent(IN) :: T ! [K]
  real(8) desorptionrate_ice
  real(8) Sice, nu

  nu = 2.2e17/sqrt(T)
  Sice = thetam*nu*exp(-Eice/kBeV/T)  ! [molecules/m^2/s]

  desorptionrate_ice = Sice
end function desorptionrate_ice


subroutine desorptionrateXi(NS,T,ThetaS,S)
  ! spectral version of function desorptiorate
  use desorp
  implicit none
  integer, intent(IN) :: NS ! number of energy bands
  real(8), intent(IN) :: T ! [K]
  real(8), intent(IN) :: ThetaS(NS) ! [molecules/m^2]
  real(8), intent(OUT) :: S(NS)  ! [molecules/m^2/s/eV]
  real(8), parameter :: maxE = 1.5  ! maximum desorption energy considered [eV]
  integer k
  real(8) nu, gamma, b, dE, vtotal, Eprime(NS), ThetaSprime(NS)
  real(8) n(NS) ! distribution of sites
  real(8) E(NS) ! array of energies

  nu = 2.2e17/sqrt(T)
  dE = maxE/real(NS)
  do k=1,NS
     E(k) = k*dE  ! 0...1.5 eV
  end do
  vtotal = sum(ThetaS)*dE/thetam
  !print *,'vtotal=',vtotal
  
  b = 1./(kBeV*T) ! frequently used factor

  do k=1,NS
     if (E(k)<Ec) then
        n(k) = 0.
     else
        n(k) = 1./W * exp(-(E(k)-Ec)/W)
        !v = ThetaS(k) / (n(k)*thetam)
     end if
  enddo
  
  if (vtotal<=1.001) then
     S(:) = nu * ThetaS(:) * exp(-E(:)*b)
  else ! multi-layer, untested
     ! (Eprime-Eice)=gamma*(E-Eice)
     gamma = 1./vtotal**2
     Eprime(:) = Eice + gamma*(E(:)-Eice)
     ThetaSprime(:) = mod(ThetaS(:),n*thetam)  ! ????
     S(:) = nu * ThetaSprime(:) * exp(-Eprime*b)
  end if

  !print *,'sum(Xi)=',sum(S)*dE
end subroutine desorptionrateXi


subroutine adsorptionXi(NS,E,ThetaS,nArea)
  ! end-member where adsorbed H2O undergoes no surface diffusion
  ! adsorption rate spectral density [molecules/m^2/s]
  ! adsorption uniform in area, i.e., proportional to n(E)
  implicit none
  integer, intent(IN) :: NS
  real(8), intent(IN) :: E(NS)
  real(8), intent(IN) :: ThetaS(NS)  ! [molecules/m^2/eV]
  real(8), intent(OUT) :: nArea(NS)  ! [sites/m^2/eV]
  real(8), parameter :: thetam=1e19, Eice=0.529, Ec=0.65, W=0.22
  integer k
  real(8) dE, vtotal, n(NS), gamma, Eprime(NS)
  
  dE = 1.5/NS   ! maxE/NS
  vtotal = sum(ThetaS)*dE/thetam
  
  if (vtotal<1.001) then
     do k=1,NS
        if (E(k)<Ec) then
           n(k)=0.
        else
           n(k) = 1./W * exp(-(E(k)-Ec)/W)
        end if
     end do
     !assert( ThetaS(k) <= 1.001*n(k)*thetam
  else ! gamma<1
     gamma = 1./vtotal**2
     Eprime = Eice + gamma*(E-Eice)
     do k=1,NS
        if (E(k) < Eice+gamma*(Ec-Eice)) then
           n(k) = 0.
        else
           n(k) = 1./(W*gamma) * exp( (Eice-E(k)+gamma*(Ec-Eice) ) / (W*gamma))
        end if
     end do
  end if
  
  nArea(:) = n
end subroutine adsorptionXi


subroutine distributor(NS,nmono,ThetaS)
  ! redistribute adsorbate mass when it exceeds monolayer
  ! for end-member with otherwise no surfacee diffusion
  implicit none
  integer, intent(IN) :: NS
  real(8), intent(IN) :: nmono(NS)  ! n*thetam
  real(8), intent(INOUT) :: ThetaS(NS)
  integer k
  real(8) abovem
  
  abovem=0.
  do k=1,NS
     if (ThetaS(k)>nmono(k)) then
        abovem = abovem + (Thetas(k)-nmono(k))
        ThetaS(k) = nmono(k)
     end if
  end do
  if (abovem>0.) print *,'excess redistribution'
  k = NS
  do while (abovem>0. .and. k>=1)
     if (thetaS(k)<nmono(k)) then
        if (thetaS(k)+abovem < nmono(k)) then
           thetaS(k) = thetaS(k)+abovem
           abovem = 0.
        else
           abovem = abovem - (nmono(k)-thetaS(k))
           thetaS(k) = nmono(k)
        end if
     endif
     k = k-1
  end do
  if (abovem>0.) then
     stop 'maxed everything out'
  end if
end subroutine distributor
