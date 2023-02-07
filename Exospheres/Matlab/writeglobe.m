function writeglobe(unit,Tsurf)
  % ouptput variables on geographic grid

  global nlon nlat
  % must match temperature grid   
  dlon = 360/nlon;
  longitude = [0.5:1:nlon]*dlon;
  dlat = 180/nlat;
  latitude = 90 - [0.5:1:nlat]*dlat;
  
  for j=1:nlat
    for i=1:nlon
      fprintf(unit,"%5.1f %6.2f %7.3f\n",longitude(i),latitude(j),Tsurf(i,j));
    end
  end
end
