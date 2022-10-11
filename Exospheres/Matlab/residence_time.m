function tau  = residence_time(T)
  % returns molecular residence time [s]
  
  sigma0 = 1e19;  % monolayer H2O [molecules/m^2]
  kB = 1.38065e-23;
  mu = 18.015*1.66054e-27;  % [kg]

  % vapor pressure of ice
  % eq. (2) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
  A=-6143.7; B=28.9074;
  psv = exp(A./T+B);  % Clapeyron

  % sublimation rate of H2O [#molecules/m^2/s]
  sublrate = psv ./ sqrt(2*pi*kB*T*mu);

  tau = sigma0 ./ sublrate;
  
  if T==0.,
    tau = 1e32;
  end
  
end 




