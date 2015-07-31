      real*8 function flux(R,decl,latitude,HA,albedo,
     &     fracir,fracdust,surfaceSlope,azFac)
C***********************************************************************
C   flux:  program to calculate solar insolation
C
C     R: distance from sun (AU)
C     decl: planetocentric solar declination (radians)
C     latitude: (radians)
C     HA: hour angle (radians from noon, clockwise)
C     fracir: fraction of absorption
C     fracdust: fraction of scattering
C     surfaceSlope: >0, (radians) 
C     azFac: azimuth of gradient (radians east of north)
C
C   written by Samar Khatiwala and Norbert Schorghofer
C***********************************************************************
      implicit none
      real*8 pi, So, d2r
      parameter (pi=3.1415926535897931, So=1365., d2r=pi/180.)
      real*8 R,decl,latitude,HA,albedo,fracIR,fracDust
      real*8 surfaceSlope,azFac
      real*8 c1,s1,Qo,solarAttenuation,Q
      real*8 sinbeta,sinbetaNoon,cosbeta,sintheta,azSun,buf

!      solarAttenuation = (1.-albedo)*(1.-fracDust) 
      c1=cos(latitude)*cos(decl)
      s1=sin(latitude)*sin(decl)
C     beta = 90 minus incidence angle for horizontal surface
C     beta = elevation of sun above (horizontal) horizon 
      sinbeta = c1*cos(HA) + s1
      sinbetaNoon = c1 + s1
      sinbetaNoon = max(sinbetaNoon,0.d0)  ! horizon
      
      cosbeta = sqrt(1-sinbeta**2)
C     ha -> az (option 1)
!      azSun=asin(-cos(decl)*sin(ha)/cosbeta)
C     ha -> az (option 2)
      buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta)
      if (buf>+1.) buf=1.d0   ! damn roundoff
      if (buf<-1.) buf=-1.d0  ! damn roundoff
      azSun = acos(buf)
      if (sin(HA)>=0) azSun=2*pi-azSun
C      ha -> az (option 3)  without beta
!      azSun=sin(latitude)*cos(decl)*cos(ha)-cos(latitude)*sin(decl)
!      azSun=atan(sin(ha)*cos(decl)/azSun)
!      print *,asin(sinbeta)/d2r,azSun/d2r

C     theta = 90 minus incidence angle for sloped surface
      sintheta = cos(surfaceSlope)*sinbeta +
     &     sin(surfaceSlope)*cosbeta*cos(azSun-azFac)
      sintheta = max(sintheta,0.d0)  ! horizon
      if (sinbeta<0.) sintheta=0.  ! horizontal horizon at infinity

      Qo = So/(R**2)
      solarAttenuation = (1.-albedo)*
     &     (1.- fracIR - fracDust)**(1./max(sinbeta,0.04))  ! new since May '05
C     net flux: solar insolation + IR
      Q = sintheta*solarAttenuation +   
     &     cos(surfaceSlope/2.)**2*sinbetaNoon*fracIR ! new since May '05
      if (sinbeta>0.d0) then
         Q = Q + 0.5*cos(surfaceSlope/2.)**2*fracDust ! new since May '05
      endif
      flux=Qo*Q

!      write(6,'(99(1x,f6.2))') decl/d2r,HA/d2r,flux,
!     &     asin(sintheta)/d2r,asin(sinbeta)/d2r,azSun/d2r,buf
      end

