      subroutine setgrid(nz,z,zmax,zfac)
C     construct regularly or geometrically spaced 1D grid
C     z(n)-z(1) = 2*z(1)*(zfac**(n-1)-1)/(zfac-1)
      implicit none
      integer NMAX
      parameter (NMAX=1000)
      integer nz, i
      real*8 zfac, zmax, z(NMAX), dz
      do i=1,nz
         dz = zmax/nz
         z(i) = (i-0.5)*dz
      enddo
      if (zfac>1.) then
         dz=zmax/(3.+2.*zfac*(zfac**(nz-2)-1.)/(zfac-1.))
         z(1)=dz
         z(2)=3*z(1)
         do i=3,nz
            z(i)=(1.+zfac)*z(i-1)-zfac*z(i-2)
         enddo
      endif
      if (z(1)<1.e-5) then
      !   print *,'WARNING: first grid point is too shallow'
      endif
      end



      real*8 function smartzfac(nz_max,zmax,nz_delta,delta)
C     output can be used as input to setgrid
C     produces zfac with desired number of points within depth delta
      implicit none
      integer nz_max, nz_delta
      real*8 zmax, delta
      integer j
      real*8 f,g,gprime
      if (nz_max<nz_delta .or. nz_delta<=1) then
         stop 'inappropriate input to smartzfac'
      endif
      f=1.05
      do j=1,7  ! Newton iteration
         !print *,j,f
         g = (f-3+2*f**(nz_max-1))/(f-3+2*f**(nz_delta-1))-zmax/delta
         gprime=(1+2*(nz_max-1)*f**(nz_max-2))/(f-3+2*f**(nz_delta-1)) - 
     &        (f-3+2*f**(nz_max-1))*(1+2*(nz_delta-1)*f**(nz_delta-2))/ 
     &        (f-3+2*f**(nz_delta-1))**2
         f = f-g/gprime
      enddo
      smartzfac=f
      if (smartzfac<1. .or. smartzfac>2.) then
         print *,'zfac=',smartzfac
         stop 'unwanted result in smartzfac'
      endif
      end


      
      subroutine smartgrid(nz,z,zdepth,thIn,rhoc,porosity,ti,rhocv, 
     &     layertype,icefrac)
C***********************************************************************
C     smartgrid: returns intelligently spaced grid and appropriate 
C                values of thermal inertia ti and rho*c in icy layer
C                  
C     INPUTS: 
C             nz = number of grid points
C             z = grid spacing for dry regolith 
C                 (will be partially overwritten)
C             zdepth = depth (in meters) where ice table starts
C                      negative values indicate no ice
C             rhoc = heat capacity per volume of dry regolith [J/m^3]
C             thIn = thermal inertia of dry regolith [SI-units]
C             porosity = void space / total volume
C
C     OUTPUTS: z = grid coordinates
C              vectors ti and rhocv
C
C     Earlier versions of this subroutine used layertype=1
C***********************************************************************

      implicit none
      integer NMAX
      parameter (NMAX=1000)

      integer nz, j, b, layertype
      real*8 z(NMAX), zdepth, thIn, rhoc, porosity
      real*8 ti(NMAX), rhocv(NMAX), stretch, newrhoc, newti
      real*8 NULL, icefrac
      parameter (NULL=0.)
      
      if (zdepth>0.and.zdepth<z(nz)) then
         
         if (layertype==2) then ! mostly ice
            call soilthprop(porosity,NULL,rhoc,thIn,2,newrhoc,newti,
     &           icefrac)
         elseif (layertype==3) then ! all ice
            call soilthprop(NULL,NULL,NULL,NULL,3,newrhoc,newti,NULL)
         else ! layertype =1
            call soilthprop(porosity,1.d0,rhoc,thIn,1,newrhoc,newti,
     &           NULL)
         endif
         stretch = (newti/thIn)*(rhoc/newrhoc)
         
         b=1
         do j=1,nz
            if (z(j)<=zdepth) then 
               b=j+1
               rhocv(j)=rhoc
               ti(j)=thIn
            else
               rhocv(j)=newrhoc
               ti(j)=newti
            endif
c            print *,j,z(j),ti(j),rhocv(j)
         enddo
         do j=b+1,nz
            z(j)=z(b)+stretch*(z(j)-z(b))
         enddo

c        print *,'zmax=',z(nz),' b=',b,' stretch=',stretch
c        print *,'depth at b-1, b ',z(b-1),z(b)
c        print *,'ti1=',thIn,' ti2=',newti
c        print *,'rhoc1=',rhoc,' rhoc2=',newrhoc 
      endif
      end



      subroutine smartgrid_allice(nz,z,zdepth,thIn,rhoc,ti,rhocv)
C***********************************************************************
C     smartgrid_allice: returns intelligently spaced grid for layer
C                       of dry soil on top of layer of pure ice
C                  
C     INPUTS: 
C             nz = number of grid points
C             z = grid spacing for dry regolith 
C                 (will be partially overwritten)
C             zdepth = depth (in meters) where ice table starts
C                      negative values indicate no ice
C             rhoc = heat capacity per volume of dry regolith [J/m^3]
C             thIn = thermal inertia of dry regolith [SI-units]
C
C     OUTPUTS: z = grid coordinates
C              vectors ti and rhocv
C
C     provided for backward compatibility
C***********************************************************************

      implicit none
      integer NMAX
      parameter (NMAX=1000)

      integer nz
      real*8 z(NMAX), zdepth, thIn, rhoc, ti(NMAX), rhocv(NMAX), NULL
      parameter (NULL=0.d0)

      call smartgrid(nz,z,zdepth,thIn,rhoc,NULL,ti,rhocv,3,NULL)

      end



      subroutine soilthprop(porosity,fill,rhocobs,tiobs,layertype,
     &     newrhoc,newti,icefrac)
C***********************************************************************
c     soilthprop: assign thermal properties of icy soil or dirty ice
c
c     porositiy = void space / total volume
c     rhof = density of free ice in space not occupied by regolith [kg/m^3]
c     fill = rhof/icedensity <=1 (only relevant for type 1)
c     rhocobs = heat capacity per volume of dry regolith [J/m^3]
c     tiobs = thermal inertia of dry regolith [SI-units]
c     layertype: 1=interstitial ice, 2=pure ice or ice with dirt
c     icefrac: fraction of ice in icelayer (only relevant for type 2)
c     output are newti and newrhoc
C***********************************************************************
      implicit none

      integer layertype
      real*8 porosity, fill, rhocobs, tiobs, icefrac
      real*8 newti, newrhoc
      real*8 kobs, cobs, cice, icedensity, kice
      parameter (cobs=800.d0)   ! heat capacity
      !parameter (cice=2000.d0, icedensity=926.d0, kice=2.4d0) ! uneffected by scaling
      parameter (cice=1540.d0, icedensity=927.d0, kice=3.2d0) ! at 198 Kelvin
      real*8 fA, ki0, ki, kw, k
      parameter (kw=3.)  ! Mellon et al. 1997

      kobs = tiobs**2/rhocobs
c     k, rhoc, and ti are defined in between grid points
c     rhof and T are defined on grid points

      newrhoc=-9999.
      newti  =-9999.

      if (layertype==1) then  ! interstitial ice
         !fill(1) = rhof(1)/icedensity
         !do j=2,nz
         !   fill(j) = (rhof(j)+rhof(j-1))/2./icedensity ! check indices
         !enddo
         newrhoc = rhocobs + porosity*fill*icedensity*cice
         if (fill>0.) then
c-----------linear addition (option 1)
            k = porosity*fill*kice + kobs
c-----------Mellon et al. 1997 (option 2)
            ki0 = porosity/(1/kobs-(1-porosity)/kw)
            fA = sqrt(fill)
            ki = (1-fA)*ki0 + fA*kice
            !k = kw*ki/((1-porosity)*ki+porosity*kw)
         else
            k = kobs
         endif
         newti = sqrt(newrhoc*k)
      endif

      if (layertype==2) then ! massive ice (pure or dirty ice)
         newrhoc = rhocobs*(1.-icefrac)/(1.-porosity) + 
     &        icefrac*icedensity*cice
         k = icefrac*kice + (1.-icefrac)*kw
         newti = sqrt(newrhoc*k)
      endif

      if (layertype==3) then ! all ice, special case of layertype 2, doesn't use porosity
         newrhoc = icedensity*cice
         k = kice 
         newti = sqrt(newrhoc*k)
      endif  

      if (layertype==4) then ! all pore ice, special case of layertype 1
         newrhoc = rhocobs + porosity*icedensity*cice
         k = porosity*kice + kobs
         !k = kw*kice/((1-porosity)*kice+porosity*kw)
         newti = sqrt(newrhoc*k)
      endif  

      if (layertype==5) then ! custom values
         ! values from Mellon et al. (2004) for ice-cemented soil
         newrhoc = 2018.*1040.
         k = 2.5
         newti = sqrt(newrhoc*k)
      endif  

      end




