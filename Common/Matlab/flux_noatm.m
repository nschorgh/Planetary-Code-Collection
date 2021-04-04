function Q = flux_noatm(R,decl,latitude,HA,SlopeAngle,azFac)
%**********************************************************************
%   flux_noatm: calculates incoming solar flux without atmosphere
%     R: distance from sun (AU)
%     decl: planetocentric solar declination (radians)
%     latitude: (radians)
%     HA: hour angle (radians from noon, clockwise)
%     SlopeAngle: >0, (radians) 
%     azFac: azimuth of topographic gradient (radians east of north)
%            azFac=0 is south-facing  
%**********************************************************************
  So=1365.;  % solar constant [W/m^2]
  
  c1 = cos(latitude)*cos(decl);
  s1 = sin(latitude)*sin(decl);
  % beta = 90 minus incidence angle for horizontal surface
  % beta = elevation of sun above (horizontal) horizon 
  sinbeta = c1*cos(HA) + s1;
  
  cosbeta = sqrt(1-sinbeta**2);
  % ha -> az
  buf = (sin(decl)-sin(latitude)*sinbeta) / (cos(latitude)*cosbeta);
  % buf can be NaN if cosbeta = 0
  buf = min(buf,+1.0); % roundoff
  buf = max(buf,-1.0); % roundoff
  azSun = acos(buf);
  if sin(HA)>=0,
    azSun = 2*pi-azSun;
  end

  % theta = 90 minus incidence angle for sloped surface
  sintheta = cos(SlopeAngle)*sinbeta - sin(SlopeAngle)*cosbeta*cos(azSun-azFac);
  if cosbeta==0.,
    sintheta = cos(SlopeAngle)*sinbeta;
  end
  sintheta = max(sintheta,0.);  % horizon
  sinbeta = max(sinbeta,0.);  % horizontal horizon at infinity
  
  Q = sintheta*So/(R**2);
  



