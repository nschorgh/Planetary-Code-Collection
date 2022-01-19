%***********************************************************************
% test Crank-Nicolson subroutine
%***********************************************************************
nz = 60;
Period = 88775.244*670
emiss=1.;
Fgeo = 0.2

T = zeros(nz,1);
rhocv = zeros(nz,1);
ti = zeros(nz,1);
Tmean = zeros(nz,1);

NSTEPS = 50000;
STEPSPERSOL = 120;
dt = Period / STEPSPERSOL  % time step [s]
thIn = 120.  % thermal Inertia
albedo = 0.2;
latitude = 5.  % [degree]
zmax = 2.5
zfac = 1.05;

% rhoc=thIn*sqrt(Period/pi)  % skin depth = 1
rhocv(:) = 1200.*800.;
delta = thIn/rhocv(1)*sqrt(Period/pi)
ti(:) = thIn;
  
Rau = 1.52; Decl = 0.;
  
fout1 = fopen('Tsurface','w'); % surface temperature
fout2 = fopen('Tprofile','w'); % temperature profile

% Initialize
Tsurf = 210.;
T(:) = 210.;
z = setgrid(nz,zmax,zfac);
f2 = fopen('z','w');
fprintf(f2,'%f',z(:));
fclose(f2);

latitude = latitude*pi/180.; % [deg] -> [rad]

Fmean = 0.;
  
HA = 0.;
  
% solar insolation 
Qn = (1-albedo) * flux_noatm(Rau,Decl,latitude,HA,0.,0.);

for n=0:NSTEPS
  % disp([time/Period,Qn])
  time = (n+1)*dt;   %   time at n+1; 
  HA = 2*pi*mod(time/Period,1.); %  hour angle
  Qnp1 = (1-albedo) * flux_noatm(Rau,Decl,latitude,HA,0.,0.);
  [T, Tsurf, Fsurf] = conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeo);
  Qn = Qnp1;

  if mod(n,3)==0,
    fprintf(fout1,'%12.6f %9.3f %9.3f\n',time/Period,Tsurf,T(nz) );
  end
  if n > NSTEPS-STEPSPERSOL,
    if mod(n,10)==0,
      fprintf(fout2,[repmat('%7.2f ',1,nz+1),'\n'],Tsurf,T(:) );
    end
    Fmean = Fsurf + Fmean;
    Tmean(:) = Tmean(:) + T(:);
  end
     
end

fclose(fout1);
fclose(fout2);

