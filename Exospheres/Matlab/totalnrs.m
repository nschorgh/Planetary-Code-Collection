function cc = totalnrs(p_s)

  cc(1) = sum(p_s==0);   % on surface
  cc(2) = sum(p_s==1);   % inflight
  cc(3) = sum(p_s==-1);  % destroyed, photo
  cc(4) = sum(p_s==-2);  % destroyed, escape
  cc(5) = sum(p_s==-3);  % coldtrapped, north
  cc(6) = sum(p_s==-4);  % coldtrapped, south

