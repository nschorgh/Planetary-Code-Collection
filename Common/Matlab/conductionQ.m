function [T, Tsurf, Fsurf] = conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhoc,emiss,...
					 Tsurf,Fgeotherm)
%***********************************************************************
%   conductionQ:  wrapper for cranknQ, which improves stability
%   Arguments and restrictions are the same as for subroutine cranknQ below.
%   created wrapper using flux smoothing 12/2023
%***********************************************************************
 
  Ni=5;  % for flux smoothing

  Tsurfold = Tsurf;
  Told = T(:);

  [T, Tsurf, Fsurf] = cranknQ(nz,z,dt,Qn,Qnp1,T,ti,rhoc,emiss,Tsurf,Fgeotherm);
  
  % artificial flux smoothing
  if ( Tsurf>1.2*Tsurfold || Tsurf<0.8*Tsurfold )  % linearization error
     Tsurf = Tsurfold;
     T(:) = Told(:);
     avFsurf = 0.;
     for j=1:Ni
        Qartiold = ( (Ni-j+1)*Qn + (j-1)*Qnp1 ) / Ni;
        Qarti    = ( (Ni-j)*Qn + j*Qnp1 ) / Ni;
        [T,Tsurf,Fsurf] = cranknQ(nz,z,dt/Ni,Qartiold,Qarti,T,ti,rhoc,emiss,...
				  Tsurf,Fgeotherm);
        avFsurf = avFsurf + Fsurf;
     end
     Fsurf = avFsurf/Ni;

  end
end



function [T, Tsurf, Fsurf] = cranknQ(nz,z,dt,Qn,Qnp1,T,ti,rhoc,emiss,...
				     Tsurf,Fgeotherm)
%***********************************************************************
%   crankn:  program to calculate the diffusion of temperature into the
%            ground and thermal emission at the surface with variable
%            thermal properties on irregular grid
%   Crank-Nicolson scheme, flux conservative, uses Samar's radiation formula
%
%   Eqn: rhoc*T_t = (k*T_z)_z 
%   BC (z=0): Q(t) + kT_z = em*sig*T^4
%   BC (z=L): heat flux = Fgeotherm
%
%   nz = number of grid points
%   dt = time step
%   Qn,Qnp1 = net solar insolation at time steps n and n+1 [W/m^2]
%   T = vertical temperature profile [K]  (in- and output)  
%   ti = thermal inertia [J m^-2 K^-1 s^-1/2]  VECTOR
%   rhoc = rho*c  VECTOR where rho=density [kg/m^3] and 
%                              c=specific heat [J K^-1 kg^-1]
%   ti and rhoc are not allowed to vary in the layers immediately 
%               adjacent to the surface or the bottom
%   emiss = emissivity
%   Tsurf = surface temperature [K]  (in- and output)
%   Fgeotherm = geothermal heat flux at bottom boundary [W/m^2]
%   Fsurf = heat flux at surface [W/m^2]  (output)
%
%   Grid: surface is at z=0
%         z(0)=0, z(2)=3*z(1), i.e., the top layer has half the width
%         T(1) is at z(1); ...; T(i) is at z(i)
%         k(i), rhoc(i), ti(i) are midway between z(i-1) and z(i)
%     
%   originally written by Samar Khatiwala, 2001
%   extended to variable thermal properties
%         and irregular grid by Norbert Schorghofer
%   converted to Matlab/Octave 4/2021
%   added Volterra predictor 1/2024
%***********************************************************************

   % set some constants
   k = ti.^2 ./ rhoc;  % thermal conductivity

   rho2 = zeros(nz,1);
   alpha = zeros(nz,1);
   gamma = zeros(nz,1);
  
   dz = 2.*z(1);
   beta = dt/rhoc(1)/(2.*dz^2);
   alpha(1) = beta*k(2);
   gamma(1) = beta*k(1);

   rho2(2:nz-1) = ( rhoc(2:nz-1) + rhoc(3:nz) )/2;
   alpha(2:nz-1) = k(3:nz) *dt ./ (z(3:nz)-z(1:nz-2)) ./ rho2(2:nz-1) ./ (z(3:nz)-z(2:nz-1));
   gamma(2:nz-1) = k(2:nz-1) *dt ./ (z(3:nz)-z(1:nz-2)) ./ rho2(2:nz-1) ./ (z(2:nz-1)-z(1:nz-2));
  
   alpha(nz) = 0.;
   gamma(nz) = k(nz)*dt/(2*rhoc(nz))/(z(nz)-z(nz-1))^2;
  
   k1dz = k(1)/dz;
  
   % elements of tridiagonal matrix
   a = -gamma(:);   %  a(1) is not used
   b = 1. + alpha(:) + gamma(:); %  b(1) has to be reset at every timestep
   c = -alpha(:);   %  c(nz) is not used
   a(nz) = -2*gamma(nz);
   b(nz) = 1 + 2*gamma(nz);
  
   % Volterra predictor (optional)
   sigSB=5.6704e-8;
   Fsurf = - k(1) * ( T(1)-Tsurf ) / z(1);  % heat flux
   seb = -Fsurf -emiss*sigSB*Tsurf**4 + (2*Qnp1 + Qn)/3.;
   % Tpred = Tsurf + sqrt(4*dt/pi) / ti(1) * seb;  ! 1st order
   Tpred = Tsurf + seb / ( sqrt(pi/(4.*dt))*ti(1) + 8./3.*emiss*sigSB*Tsurf**3 );
   Tr = (Tsurf+Tpred)/2.;  % better reference temperature

   % Emission
   %Tr = Tsurf;            % 'reference' temperature
   arad = -3.*emiss*sigSB*Tr^4;
   brad = 2.*emiss*sigSB*Tr^3;
   ann = (Qn-arad)/(k1dz+brad);
   annp1 = (Qnp1-arad)/(k1dz+brad);
   bn = (k1dz-brad)/(k1dz+brad);
   b(1) = 1. + alpha(1) + gamma(1) - gamma(1)*bn;
  
   % Set RHS         
   r = zeros(nz,1);
   r(1) = gamma(1)*(annp1+ann) + ...
	  (1.-alpha(1)-gamma(1)+gamma(1)*bn)*T(1) + alpha(1)*T(2);
   r(2:nz-1) = gamma(2:nz-1).*T(1:nz-2) + ...
	       (1-alpha(2:nz-1)-gamma(2:nz-1)).*T(2:nz-1) + ...
	       alpha(2:nz-1).*T(3:nz);
   r(nz) = 2*gamma(nz)*T(nz-1) + (1.-2*gamma(nz))*T(nz) ...
	   + 2*dt/rhoc(nz)*Fgeotherm/(z(nz)-z(nz-1)); 
  
   % Solve for T at n+1
   D = [ [a(2:nz);0], b, [0;c(1:nz-1)] ];
   A = spdiags(D,-1:1,nz,nz);
   T = A\r;  % update by tridiagonal inversion
   
   Tsurf = 0.5*(annp1 + bn*T(1) + T(1)); % (T0+T1)/2

   Fsurf = - k(1) * (T(1)-Tsurf) / z(1); % heat flux into surface

end




