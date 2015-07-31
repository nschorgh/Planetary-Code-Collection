      subroutine setgrid(nz,z,zmax,zfac)
C     construct regularly or geometrically spaced 1D grid
      implicit none
      integer NMAX
      parameter (NMAX=2000)
      integer nz, i
      real*8 zfac, zmax, z(NMAX), dz
      dz = zmax/nz
      do i=1,nz
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


      
      subroutine smartgrid(nz,z,zdepth,thIn,rhoc,porosity,ti,rhocv)
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
C***********************************************************************

      implicit none
      integer NMAX
      parameter (NMAX=2000)

      integer nz, j, b
      real*8 z(NMAX), zdepth, thIn, rhoc, porosity
      real*8 ti(NMAX), rhocv(NMAX), stretch, newrhoc, newti

      if (zdepth>0.and.zdepth<z(nz)) then

         !newrhoc = rhoc + porosity*icedensity*cice
         !newti = sqrt(newrhoc*(porosity*kice+thIn**2/rhoc))
         call soilthprop(porosity,1.d0,rhoc,thIn,1,newrhoc,newti,0.d0)
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
C***********************************************************************

      implicit none
      integer NMAX
      parameter (NMAX=2000)

      integer nz, j,b
      real*8 z(NMAX), stretch, zdepth, thIn, rhoc
      real*8 ti(NMAX), rhocv(NMAX), newrhoc, newti, NULL
c      real*8 cice, icedensity, kice
c      parameter (cice=1540.d0, icedensity=927.d0, kice=3.2d0)
c      parameter (cice=1040.d0, icedensity=2018.d0, kice=2.5d0) ! Mellon's values
      parameter (NULL=0.d0)

      if (zdepth>0.and.zdepth<z(nz)) then

         !newrhoc = icedensity*cice
         !newti = sqrt(icedensity*cice*kice)
         call soilthprop(NULL,NULL,NULL,NULL,2,newrhoc,newti,1.d0)
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
         enddo
         do j=b+1,nz
            z(j)=z(b)+stretch*(z(j)-z(b))
         enddo

      endif
      end



      subroutine soilthprop(porosity,fill,rhocobs,tiobs,layertype,
     &     newrhoc,newti,icefrac)
C***********************************************************************
c     soilthprop: assign thermal properties of icy soil or dirty ice
c
c     porositiy = void space / total volume
c     rhof = density of free ice in space not occupied by regolith [kg/m^3]
c     fill = rhof/icedensity (only relevant for type 1)
c     rhocobs = heat capacity per volume of dry regolith [J/m^3]
c     tiobs = thermal inertia of dry regolith [SI-units]
c     layertype: 1=interstitial ice, 2=ice with dirt
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
      !parameter (cice=1040.d0, icedensity=2018.d0, kice=2.5d0) ! Mellon's values
      real*8 fA, ki0, ki, kw, k
      parameter (kw=3)  ! Mellon et al. 1997

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
            !k = porosity*fill*kice + kobs
c-----------Mellon et al. 1997 (option 2)
            ki0 = porosity/(1/kobs-(1-porosity)/kw)
            fA = sqrt(fill)
            ki = (1-fA)*ki0 + fA*kice
            k = kw*ki/((1-porosity)*ki+porosity*kw)
         else
            k = kobs
         endif
         newti = sqrt(newrhoc*k)
      endif

      if (layertype==2) then 
         newrhoc = rhocobs*(1-icefrac)/(1-porosity) + 
     &        icefrac*icedensity*cice
         !k = icefrac*kice + (1-icefrac)/(1-porosity)*kobs ! my old way
         k = icefrac*kice + (1-icefrac)*kw
         newti = sqrt(newrhoc*k)
      endif  

      end




