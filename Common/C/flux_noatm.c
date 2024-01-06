#include <math.h>


double flux_noatm(double R, double decl, double latitude, double HA,
		  double surfaceSlope, double azFac) {
/**********************************************************************
   flux_noatm: calculates incoming solar flux without atmosphere
     R: distance from sun (AU)
     decl: planetocentric solar declination (radians)
     latitude: (radians)
     HA: hour angle (radians from noon, clockwise)
     surfaceSlope: >0, (radians) 
     azFac: azimuth of topographic gradient (radians east of north)
***********************************************************************/
  const double  So=1365.;  // solar constant
  const double  pi=3.1415926535897931;
  double c1,s1,sinbeta,cosbeta,sintheta,azSun,buf;
    
  c1 = cos(latitude)*cos(decl);
  s1 = sin(latitude)*sin(decl);
  // beta = 90 minus incidence angle for horizontal surface
  // beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1;
  
  cosbeta = sqrt(1-sinbeta*sinbeta);
  /* ha -> az (option 1) */
  // azSun=asin(-cos(decl)*sin(HA)/cosbeta);
  /* ha -> az (option 2) */
  buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta);
  /* buf can be NaN if cosbeta = 0 */
  if (buf>+1.) buf=+1.0;  // damn roundoff
  if (buf<-1.) buf=-1.0;  // damn roundoff
  azSun = acos(buf);
  if (sin(HA)>=0) azSun = 2*pi-azSun;
  /* ha -> az (option 3)  without beta */
  //azSun=sin(latitude)*cos(decl)*cos(HA)-cos(latitude)*sin(decl);
  //azSun=atan(sin(HA)*cos(decl)/azSun);

  /* theta = 90 minus incidence angle for sloped surface */
  sintheta = cos(surfaceSlope)*sinbeta + 
    sin(surfaceSlope)*cosbeta*cos(azSun-azFac);
  if (cosbeta==0.) sintheta = cos(surfaceSlope)*sinbeta;
  if (sintheta<0.) sintheta = 0.;  // horizon
  if (sinbeta<0.) sintheta = 0.;  // horizontal horizon at infinity
  
  return sintheta*So/(R*R);
}
