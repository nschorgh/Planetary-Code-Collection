function [p_r, p_s, p_t] = hop1(p_r, p_t, Tsurf, Q)
  % ballistic flight of one particle
  %    p_r(2)   (1)=longitude  (2)=latitude
  %    p_s      status
  %    p_t      time
  CORIOLIS = true;  % turn Coriolis effect on or off

  % photodissociation time scale at 1 AU
  %taudissoc = 20.*3600.;  % Potter & delDuca (1964)
  taudissoc = 1/12.6e-6;  % Crovisier (1989)

  mmass = 18.015;   % H2O
  global semia g Rbody vescape siderealDay
  
  % Maxwell-Boltzmann launch velocities
  sigma = sqrt(Tsurf*8314.46/mmass);  % standard deviation
  v = randn(3,1)*sigma;  % randn has unit variance
  v(3) = abs(v(3)); % vertical velocity component always upward
  
  % Armand distribution launch velocities; P(vz)=vz/sigma^2*exp(-vz^2/2*sigma^2)
  % see e.g. Devroye, p29
  v(3) = sqrt(2.)*sigma*sqrt(-log(rand));
  % use the same v(1), v(2), and sigma as for Maxwell-Boltzmann
  
  if CORIOLIS,
    v(1) = v(1) - 2*Rbody*pi/siderealDay*cosd(p_r(2));
  end
  vspeed = sqrt(sum(v(:).^2));
  if vspeed>vescape, % gravitational escape
    p_s = -2;  % lost due to escape
    p_t = 1e100;  % never use again
    return
  end

  % molecule moves on plane that goes through center of sphere/body
  % ground track is thus part of a great circle
  flighttime = 2*v(3)/g;   % time of flight for constant g
  d = 2/g*v(3)*sqrt(v(1)^2+v(2)^2);  % distance for constant g
  if vspeed>0.4*vescape,  % use non-uniform gravity formula
  %if d>0.1*Rbody then
    alpha = atan(sqrt(v(1)^2+v(2)^2)/v(3));  % angle from zenith
    [d,flighttime] = nonuniformgravity(vspeed,alpha);
  end
  %disp([flighttime,d])  % for flighttime statistics
  %az = atan2(v(2),v(1))   % for statistical tests
  cosaz = v(2) / sqrt(v(1)^2+v(2)^2);
  lat = deg2rad( p_r(2) );
  sinph2 = sin(d/Rbody)*cos(lat)*cosaz + sin(lat)*cos(d/Rbody);
  p_r(2) = asind(sinph2);
  cosph2 = sqrt(1.-sinph2^2);
  if cosph2~=0,  % not on pole
    cosdlon = (cos(d/Rbody)*cos(lat)-sin(lat)*sin(d/Rbody)*cosaz)/cosph2;
    if (cosdlon>+1.), cosdlon=+1.; end  % roundoff
    if (cosdlon<-1.), cosdlon=-1.; end  % roundoff
    dlon = acos(cosdlon);
    if (d/Rbody>pi)
      dlon = 2*pi-dlon;
    end
    if v(1)<0., dlon=-dlon; end
    p_r(1) = p_r(1) + rad2deg(dlon);
  else   % on pole
     % longitude does not matter 
    p_r(1) = 0.;  % just in case
  end
  if CORIOLIS,
    p_r(1) = p_r(1) + flighttime/siderealDay*360.;
  end
    
  if p_r(2)>90. || p_r(2)<-90,
    error(['hop1: this cannot happen',p_r(2)])
  end
  p_r(1) = mod(p_r(1),360.);   % 0 <= p_r(1) < 360.
  
  p_s = 1;
  p_t = p_t + flighttime;

  % in-flight destruction
  if Q>0 && taudissoc>0.,
    u = rand;  % between 0 and 1
    destr_rate = flighttime/(taudissoc*semia^2);  % photodissociation
    if destr_rate>0.2,
      destr_rate = 1-exp(-destr_rate);
    end
    if u < destr_rate,
      p_s = -1;  % destroyed by photodissociation
      p_t = 1e100;  % never use again
    end
  end
end  % end of function hop1



function [d,t] = nonuniformgravity(vspeed,alpha)
  % ballistic travel distance and flighttime for non-uniform gravity
  % not suitable for small velocities due to roundoff
  % vspeed ... launch speed (m/s)
  % alpha ... zenith angle of launch velocity (radian)
  % d ... flight distance (m)
  % t ... flighttime (s)
  global vescape Rbody
  
  gamma = (vspeed/vescape)^2;
  %a = Rbody/2./(1-gamma);
  ecc = sqrt(1-4*(1-gamma)*gamma*sin(alpha)^2);
  %theta = 1/ecc*(1-2*gamma*sin(alpha)^2);
  d = 2*Rbody*acos(1/ecc*(1-2*gamma*sin(alpha)^2));
  %d = 2*Rbody*acos(theta);
  %Ep = pi - 2*atan(sqrt((1-ecc)/(1+ecc))/tan(d/(4*Rbody)));
  Ep = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(d/(4*Rbody))) ;
  %Ep = 2*atan(sqrt((1+ecc)/(1-ecc)*(1-theta)/(1+theta))) ;
  if (ecc>1.-1e-5),
    d = Rbody*4*gamma*sin(alpha);
    %Ep = pi - 2*atan(sqrt((1-gamma)/gamma));
    Ep = 2*atan(sqrt(gamma/(1-gamma)));
  end
  %t = 2*sqrt(2*a^3/Rbody/vescape^2)*(Ep+ecc*sin(Ep))
  t = (Rbody/vescape)/(1-gamma)^1.5*(Ep+ecc*sin(Ep));
  if 1-2*gamma*sin(alpha)^2 > ecc, % otherwise d=NaN
    d = 0.;
    t = 0.;
  end
end
