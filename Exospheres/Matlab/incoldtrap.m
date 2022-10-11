function inside = incoldtrap(p_r)

  inside = false;

  % approx. relative area of spherical cap (a*pi/180)**2/2, a=sqrt(2*F)*180/pi
  % relative area of spherical cap 1-cos(a*pi/180), a=acos(1-F)*180/pi

  % CERES
  if p_r(2)> +90-2.92, % 0.13% of hemisphere
    inside = true;
  end
  if p_r(2)< -90+2.92, % 0.13% of hemisphere
    inside = true;
  end
  % dlat = rad2deg( 0.13e-2/cosd(80.) )
  %if (abs(p_r(2))>80.-dlat/2 .and. abs(p_r(2))<80.+dlat/2.) incoldtrap = true

end

