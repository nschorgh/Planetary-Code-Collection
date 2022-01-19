function [T, Fsurf] = conductionT(nz, z, dt, T, Tsurf, Tsurfp1, ti, rhoc, Fgeotherm)
%***********************************************************************
%   conductionT:  program to calculate the diffusion of temperature 
%                 into the ground with prescribed surface temperature 
%                 and variable thermal properties on irregular grid
%   Crank-Nicholson scheme, flux conservative
%
%   Eqn: rhoc*T_t = (k*T_z)_z 
%   BC (z=0): T=T(t)
%   BC (z=L): heat flux = Fgeotherm
%
%   nz = number of grid points
%   dt = time step
%   T = vertical temperature profile [K]  (in- and output)
%   Tsurf, Tsurfp1 = surface temperatures at times n and n+1  
%   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
%   rhoc = rho*c  heat capacity per volume [J m^-3 K^-1]  VECTOR
%   ti and rhoc are not allowed to vary in the layers immediately
%               adjacent to the surface or the bottom
%   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
%   Fsurf = heat flux at surface [W/m^2]  (output)
%
%   Grid: surface is at z=0
%         T(1) is at z(1); ...; T(i) is at z(i)
%         k(i) is midway between z(i-1) and z(i)
%         rhoc(i) is midway between z(i-1) and z(i)
%***********************************************************************

  % set some constants
  k = ti.^2 ./ rhoc;  % thermal conductivity
  rho2 = ( rhoc + shift(rhoc,-1) )/2;
  alpha = shift(k,-1) *dt ./ (shift(z,-1)-shift(z,+1)) ./ rho2(:) ./ (shift(z,-1)-z(:));
  gamma = k(:) *dt ./ (shift(z,-1)-shift(z,+1)) ./ rho2(:) ./ (z(:)-shift(z,+1));
  
  alpha(1) = k(2)*dt/rhoc(1)/(z(2)-z(1))/z(2);
  gamma(1) = k(1)*dt/rhoc(1)/z(1)/z(2);
  alpha(nz) = 0.;
  gamma(nz) = k(nz)*dt/(rhoc(nz)+rhoc(nz))/(z(nz)-z(nz-1))^2; % assumes rhoc(nz+1)=rhoc(nz)
  
  % elements of tridiagonal matrix
  a = -gamma(:);   %  a(1) is not used
  b = 1. + alpha(:) + gamma(:);
  c = -alpha(:);   %  c(nz) is not used

  % Set RHS
  r = gamma.*shift(T,+1) + (1-alpha-gamma).*T + alpha.*shift(T,-1);
  r(1) = alpha(1)*T(2) + (1.-alpha(1)-gamma(1))*T(1) + gamma(1)*(Tsurf+Tsurfp1);
  r(nz) = gamma(nz)*T(nz-1) + (1.-gamma(nz))*T(nz) + ...
          dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1)); % assumes rhoc(nz+1)=rhoc(nz)

  % Solve for T at n+1
  D = [ [a(2:nz);0], b, [0;c(1:nz-1)] ];
  %D = [ [-gamma(2:nz);0], 1+alpha+gamma, [0;-alpha(1:nz-1)] ];
  A = spdiags(D,-1:1,nz,nz);
  T = A\r;  % update by tridiagonal inversion
  
  Fsurf = -k(1) * (T(1)-Tsurfp1) / z(1); % heat flux into surface
  
