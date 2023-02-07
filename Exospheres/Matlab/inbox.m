function [kx,ky] = inbox(r) 
  % determine what lon/lat box coordinate r is in
  global nlon nlat
  
  dlon = 360./nlon;
  dlat = 180./nlat;
  kx = round( r(1)/dlon+0.5 );
  ky = round( (90.-r(2))/dlat+0.5 );
  if r(2)==-90.,
    ky=nlat;
  end  % roundoff issue
  %if (kx>nlon) kx=kx-nlon
  while kx>nlon,
    kx = kx-nlon;
  end
  %inbox = kx + (ky-1)*nlon;
  if (ky<1 || ky>nlat),
    error(sprintf('inbox: Index ky is out of bound %d %f\n',ky,r(2)))
  end
  if (kx<1 || kx>nlon),
    error(sprintf('inbox: Index kx is out of bound %d %f\n',kx,r(1)))
  end
end 
