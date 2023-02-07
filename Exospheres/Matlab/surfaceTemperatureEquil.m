function Tsurf = surfaceTemperatureEquil(dtsec,HAi,time)
  % equilibrium surface temperatures

  global nlon nlat
  
  % Ceres
  global semia; % = 2.7675;
  solarDay = 9.076*3600.;
  %albedo = 0.034; % Bond albedo from Li et al. (2016)
  albedo = 0.090; % V-band
  emiss = 0.95; 

  sigSB = 5.6704e-8;

  decl=0.;  % toy orbit
  sunR = semia;
  
  Qn = ones(nlon,nlat);
  
  for i=1:nlon
    lon = (i-0.5)*360/nlon;
    % longitude of morning terminator = -time/solarDay + lon + HAi ??
    HA = 2.*pi*mod( time/solarDay+(lon-HAi)/360., 1.);  % hour angle    
    for j=1:nlat
      lat = 90-(j-0.5)*180/nlat;
      Qn(i,j) = (1-albedo)*flux_noatm(sunR,decl,lat,HA,0.,0.);
    end
  end
  Tsurf(:,:) = ( Qn(:,:) / (1-albedo) / emiss / sigSB ).^0.25;

end
