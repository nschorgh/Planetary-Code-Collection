!*****************************************************
! Commonly used subroutines for fast method
! written by Norbert Schorghofer 2007-2011
!*****************************************************

pure function zint(y1,y2,z1,z2)
  ! interpolate between two values, y1*y2<0
  implicit none
  real(8), intent(IN) :: y1,y2,z1,z2
  real(8) zint
  zint = (y1*z2 - y2*z1)/(y1-y2)
end function zint



pure function equildepth(nz, z, rhosatav, rhosatav0, avrho1)
!***********************************************************************
!  returns equilibrium depth for given ice content
!  this is not the true (final) equilibrium depth
!***********************************************************************
  use allinterfaces, except_this_one => equildepth
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz), rhosatav(nz)
  real(8), intent(IN) :: rhosatav0, avrho1
  integer i, typeE
  real(8) equildepth  ! =zdepthE
  !real(8), external :: zint  ! defined in allinterfaces.mod
  
  typeE = -9; equildepth = -9999.
  do i=1,nz 
     if (rhosatav(i) <= avrho1) then
        typeE=i
        exit
     endif
  enddo
  if (typeE>1) then ! interpolate
     equildepth=zint(avrho1-rhosatav(typeE-1),avrho1-rhosatav(typeE),z(typeE-1),z(typeE))
  end if
  if (typeE==1) equildepth=zint(avrho1-rhosatav0,avrho1-rhosatav(1),0.d0,z(1))
  if (equildepth>z(nz)) equildepth=-9999.   ! happens very rarely
end function equildepth



subroutine depths_avmeth(typeT, nz, z, rhosatav, rhosatav0, rlow, avrho1,  &
     & porefill, typeF, zdepthF, B, ypp, typeG, zdepthG)
!***********************************************************************
!  returns interface depth and ypp
!  also returns lower geothermally limited boundary, if applicable
!
!  this is a crucial subroutine
!
!  B = Diff*bigstep/(porosity*icedensity)  [SI units]
!***********************************************************************
  use allinterfaces, except_this_one => depths_avmeth
  implicit none
  integer, intent(IN) :: nz, typeT
  real(8), intent(IN), dimension(nz) :: z, rhosatav, porefill
  real(8), intent(IN) :: rhosatav0, rlow, avrho1
  integer, intent(OUT) :: typeF  ! index of depth below which filling occurs
  real(8), intent(INOUT) :: zdepthF
  real(8), intent(IN) :: B 
  real(8), intent(OUT) :: ypp(nz), zdepthG
  integer, intent(INOUT) :: typeG  ! positive on input when Fgeotherm>0

  integer i, typeP, nlb, newtypeG
  real(8) eta(nz), Jpump1, help(nz), yp(nz), zdepthFold, ap_one, ap(nz)
  real(8) leak, cumfill, cumfillabove

  if (typeT<0) then
     nlb = nz
     do i=1,nz
        eta(i) = constriction(porefill(i))
     enddo
  else
     !nlb = typeT-1
     nlb = typeT ! changed 2010-09-29
     do i=1,typeT-1
        eta(i) = constriction(porefill(i))
     enddo
     do i=typeT,nz
        eta(i)=0.
     enddo
  end if

!-fill depth
  zdepthFold = zdepthF
  typeF = -9;  zdepthF = -9999.
  call deriv1(z,nz,rhosatav,rhosatav0,rlow,yp)  ! yp also used below
  do i=1,nlb
     Jpump1 = (rhosatav(i)-avrho1)/z(i)  ! <0 when stable
     ! yp is always <0
     help(i) = Jpump1 - eta(i)*yp(i)
     leak = porefill(i)/B*(z(i)-zdepthFold)/(18./8314.)
     !help(i) = help(i)-leak   ! optional
     if (help(i) <= 0.) then
        typeF=i
        !print *,'#',typeF,Jpump1,eta(typeF)*yp(typeF),leak
        exit
     endif
  enddo
  if (typeF>1) zdepthF = zint(help(typeF-1),help(typeF),z(typeF-1),z(typeF))
  if (typeF==1) zdepthF=z(1)


!-depth to shallowest perennial ice
  typeP = -9 
  do i=1,nz
     if (porefill(i)>0.) then
        typeP = i  ! first point with ice
        exit
     endif
  enddo

!-calculate ypp
  !call deriv1(z,nz,rhosatav,rhosatav0,rlow,yp)
  call deriv1(z,nz,eta(:),1.d0,eta(nz-1),ap)
  if (typeP>0 .and. typeP<nz-2) then
     ap_one=deriv1_onesided(typeP,z,nz,eta(:))
     ! print *,typeP,ap(typeP),ap_one
     ap(typeP)=ap_one
  endif
  call deriv2_simple(z,nz,rhosatav(1:nz),rhosatav0,rlow,ypp(:))
  !call deriv2_full(z,nz,eta(:),rhosatav(:),1.d0,rhosatav0,rhosatav(nz-1),ypp(:))
  !write(40,*) rhosatav
  !write(41,*) yp
  !write(42,*) ypp
  ypp(:) = ap(:)*yp(1:)+eta(:)*ypp(:)
  !write(43,*) ypp
  !write(44,*) eta(1:nz)
  !write(45,*) (rhosatav(:)-avrho1)/z(:)
  ypp(:) = ypp(:)*18./8314.
  ! ypp values beyond nlb should not be used

!-geothermal stuff
  if (typeT>0) typeG=-9
  if (typeG<0) zdepthG=-9999.
  if (typeG>0 .and. typeT<0) then
     typeG=-9
     do i=2,nz
        if (yp(i)>0.) then  ! first point with reversed flux
           typeG=i
           zdepthG=zint(yp(i-1),yp(i),z(i-1),z(i))
           !zdepthG=z(typeG)
           exit
        endif
     enddo
  else
     typeG = -9
  endif
  if (typeG>0 .and. typeT<0) then
     cumfillabove = colint(porefill(:)/eta(:),z,nz,typeG-1,nz) 
     newtypeG = -9
     do i=typeG,nz
        if (minval(eta(i:nz))<=0.) cycle  ! eta=0 means completely full
        cumfill=colint(porefill(:)/eta(:),z,nz,i,nz)
        if (cumfill<yp(i)*18./8314.*B) then  ! usually executes on i=typeG
           if (i>typeG) then
              write(34,*) '# adjustment to geotherm depth by',i-typeG
              zdepthG = zint(yp(i-1)*18./8314.*B-cumfillabove, &  ! no good
                   &        yp(i)*18./8314.*B-cumfill,z(i-1),z(i))
              if (zdepthG>z(i) .or. zdepthG<z(i-1)) write(34,*) &
                   & '# WARNING: zdepthG interpolation failed',i,z(i-1),zdepthG,z(i)
              newtypeG=i
           endif
           ! otherwise leave zdepthG untouched
           exit
        endif
        cumfillabove = cumfill
     enddo
     if (newtypeG>0) typeG=newtypeG
  end if
  ! if typeG>0, then all ice at and below typeG should be erased 
end subroutine depths_avmeth



pure function constriction(porefill)
! specify constriction function here, 0<=eta<=1
  implicit none
  real(8), intent(IN) :: porefill
  real(8) eta, constriction
  if (porefill<=0.) eta = 1.
  if (porefill>0. .and. porefill<1.) then
     ! eta = 1.
     ! eta = 1-porefill
     eta = (1-porefill)**2  ! Hudson et al., JGR, 2009
  endif
  if (porefill>=1.) eta = 0.
  constriction = eta
end function constriction



pure subroutine icechanges_poreonly(nz,z,typeF,typeG,avdrhoP,ypp,B,porefill)
  use allinterfaces, except_this_one => icechanges_poreonly
  implicit none
  integer, intent(IN) :: nz, typeF, typeG
  real(8), intent(IN) :: z(nz), ypp(nz), avdrhoP, B
  real(8), intent(INOUT) :: porefill(nz)
  integer j, erase, newtypeP, ub
  real(8) integ
  
  !----retreat
  ! avdrhoP>0 is outward loss from zdepthP
  ! avdrhoP<0 means gain at zdepthP or no ice anywhere
  if (avdrhoP>0.) then
     erase=0
     do j=1,nz
        if (typeF>0 .and. j>=typeF) exit ! don't retreat beyond typeF
        integ = colint(porefill(1:nz)*z(1:nz),z(1:nz),nz,1,j)
        erase = j
        if (integ>B*avdrhoP*18./8314.) exit
     end do
     if (erase>0) porefill(1:erase)=0.
  endif

  ! new depth
  newtypeP = -9 
  do j=1,nz
     if (porefill(j)>0.) then
        newtypeP = j  ! first point with ice
        exit
     endif
  enddo

  !----diffusive filling
  ub = typeF
  if (newtypeP>0 .and. typeF>0 .and. newtypeP<ub) ub=newtypeP
  if (ub>0) then  
     do j=ub,nz
        ! B=Diff/(porosity*icedensity)*86400*365.24*bigstep
        porefill(j) = porefill(j) + B*ypp(j)
        if (porefill(j)<0.) porefill(j)=0.
        if (porefill(j)>1.) porefill(j)=1.
     enddo
  end if
  
  !----enact bottom boundary
  if (typeG>0) porefill(typeG:nz)=0.
  
end subroutine icechanges_poreonly



pure subroutine icechanges(nz,z,typeF,avdrho,avdrhoP,ypp, &
     & Diff,porosity,icefrac,bigstep,zdepthT,porefill,typeG)
!***********************************************************
! advances ice table, advances interface, grows pore ice
!
! a crucial subroutine
!***********************************************************
  use miscparameters, only : icedensity
  use allinterfaces, except_this_one => icechanges
  implicit none
  integer, intent(IN) :: nz, typeF, typeG
  real(8), intent(IN) :: z(nz), ypp(nz), avdrho, avdrhoP
  real(8), intent(IN) :: Diff, porosity, icefrac, bigstep
  real(8), intent(INOUT) :: zdepthT, porefill(nz)
  integer j, erase, newtypeP, ub, typeP, typeT
  real(8) B, beta, integ

  B = Diff*bigstep*86400.*365.24/(porosity*icedensity)

  ! advance ice table, avdrho>0 is retreat
  if (zdepthT>=0. .and. avdrho>0.) then 
     typeP=-9999; typeT=-9999
     do j=1,nz
        if (z(j)>zdepthT) then
           typeT=j
           exit
        endif
     enddo
     do j=1,nz
        if (porefill(j)>0.) then
           typeP=j
           exit
        endif
     enddo
     if (typeP==typeT) then   ! new 2011-09-01
        beta = (1-icefrac)/(1-porosity)/icefrac
        beta = Diff*bigstep*beta*86400*365.24/icedensity
        zdepthT = sqrt(2*beta*avdrho*18./8314. + zdepthT**2)
     endif
  endif
  if (zdepthT>z(nz)) zdepthT=-9999.
  
  ! advance interface, avdrhoP>0 is loss from zdepthP
  if (avdrhoP>0.) then
     erase=0
     do j=1,nz
        if (typeF>0 .and. j>=typeF) exit ! don't retreat beyond typeF
        if (zdepthT>=0. .and. z(j)>zdepthT) exit 
        integ = colint(porefill(1:nz)*z(1:nz),z(1:nz),nz,1,j)
        erase = j
        if (integ>B*avdrhoP*18./8314.) exit
     end do
     if (erase>0) porefill(1:erase)=0.
  endif

  ! new depth
  newtypeP = -9 
  do j=1,nz
     if (zdepthT>=0. .and. z(j)>zdepthT) exit
     if (porefill(j)>0.) then
        newtypeP = j  ! first point with pore ice
        exit
     endif
  enddo

  ! diffusive filling
  ub = typeF
  if (newtypeP>0 .and. typeF>0 .and. newtypeP<ub) ub=newtypeP
  if (ub>0) then  
     do j=ub,nz
        porefill(j) = porefill(j) + B*ypp(j)
        if (porefill(j)<0.) porefill(j)=0.
        if (porefill(j)>1.) porefill(j)=1.
        if (zdepthT>=0. .and. z(j)>zdepthT) exit
     enddo
  end if

  ! below icetable
  if (zdepthT>=0.) then
     do j=1,nz
        if (z(j)>zdepthT) porefill(j) = icefrac/porosity
     enddo
  else
     ! geothermal lower boundary
     if (typeG>0) porefill(typeG:nz)=0.
  end if
end subroutine icechanges


subroutine assignthermalproperties(nz,thIn,rhoc, &
     &    ti,rhocv,typeT,icefrac,porosity,porefill)
!*********************************************************
! assign thermal properties of soil
!*********************************************************
  implicit none
  integer, intent(IN) :: nz
  integer, intent(IN), optional :: typeT
  real(8), intent(IN), optional :: icefrac
  real(8), intent(IN) :: thIn, rhoc
  real(8), intent(IN), optional :: porosity, porefill(nz)
  real(8), intent(OUT) :: ti(nz), rhocv(nz)
  integer j
  real(8) newrhoc, newti, fill
  real(8), parameter :: NULL=0.

  ti(1:nz) = thIn
  rhocv(1:nz) = rhoc
  if (typeT>0) then 
     call soilthprop(porosity,NULL,rhoc,thIn,2,newrhoc,newti,icefrac)
     rhocv(typeT:nz) = newrhoc
     ti(typeT:nz) = newti
  endif
  do j=1,nz 
     fill = porefill(j)   ! off by half point
     if (fill>0. .and. (typeT<0 .or. (typeT>0 .and. j<typeT))) then
        call soilthprop(porosity,fill,rhoc,thIn,1,rhocv(j),ti(j),NULL)
     endif
  enddo
end subroutine assignthermalproperties



subroutine compactoutput(unit,porefill,nz)
  implicit none
  integer, intent(IN) :: unit,nz
  real(8), intent(IN) :: porefill(nz)
  integer j
  do j=1,nz
     if (porefill(j)==0.) then
        write(unit,'(1x,f2.0)',advance='no') porefill(j)
     else
        write(unit,'(1x,f7.5)',advance='no') porefill(j)
     endif
  enddo
  write(unit,"('')")
end subroutine compactoutput

