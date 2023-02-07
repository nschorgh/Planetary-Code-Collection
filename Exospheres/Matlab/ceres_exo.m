%***********************************************************************
% H2O exosphere of airless body 
%***********************************************************************

global semia   % semi-major axis (AU)
global g       % gravitational acceleration on surface (m/s^2)
global Rbody   % nominal radius of Ceres (m) 
global vescape
global siderealDay  % for Coriolis effect

% define geographic grid for surface temperature
global nlon    % number of longitude bins
global nlat    % number of latitude bins

% Ceres
semia = 2.7675; 
g = 0.284;      % based on mass from Park et al. (2016)
Rbody = 470e3;  
vescape = sqrt(2*g*Rbody);
siderealDay = 655.728*3600.;

nlon=72; % 5 degree resolution
nlat=90; % 2 degree resolution

dtsec=600.; % time step for temperatures
NP=1000; % maximum number of computational particles
tmax=25;  % total run time in hours
HAi=0; % hour angle at time zero and subsolar longitude zero

p_r = nan(NP,2); % longitude(1) and latitude(2)
p_s = int8(zeros(NP,1))-9; % status 0=on surface, 1=inflight, <0= destroyed or trapped
p_t = zeros(NP,1); % time
p_n = int16(zeros(NP,1)); % # of hops (diagnostic only)
Qdissoc = 1.; % anything larger than zero turns on photo-destruction

disp('Model parameters:')
disp(sprintf('Time step = %f sec',dtsec))
disp(sprintf('Maximum time = %f hours or %f years',tmax,tmax/(24*365.242)))
disp(sprintf('Number of molecules = %d',NP))

% initialize random number generator to make results repeatable
if exist('OCTAVE_VERSION', 'builtin') == 0,
  rng(752)  % Matlab
else  % Octave
  randn('state',752)  % for randn in hop1
  rand('state',451)   % for rand in residence_time
end

% initialize counters
ccc = zeros(4,1);
cc = zeros(6,1);

% initial configuration
p_t(:)=0.;  p_s(:)=0;     % all particles on surface
p_r(:,1)=0.; p_r(:,2)=0.; % initial coordinates
%p_t(:)=1e99; p_s(:)=-9;   % no particles

% open output files
u30 = fopen('series','wt');
u50 = fopen('particles','wt');
u51 = fopen('particles_end','wt');
u61 = fopen('Tsurface_end','wt');

% loop over time steps
n = 0;
time = 0.;
while (time < tmax*3600.)
  time = (n+1)*dtsec;   % time at n+1 

  % Temperature
  %Tsurf = zeros(nlon,nlat) + 200;  % uniform temperature
  Tsurf = surfaceTemperatureEquil(dtsec,HAi,time);
  
  % some output
  cc = totalnrs(p_s);
  disp(sprintf('%f Next time step %d',time/3600.,sum(cc(1:2))))
  fprintf(u30,'%f %d %d %d %d %d %d\n',time/3600.,cc(1:2),ccc(1:4));

  % update residence times for new temperature
  for k = 1:NP
     if p_s(k)==0,  % molecule resides on surface
       [i,j] = inbox(p_r(k,:));
       p_t(k) = residence_time(Tsurf(i,j));  % relative to time
     end
  end
       
  if n==0, % write out initial distribution
    writeparticles(u50,NP,p_r,p_s,p_t,p_n)
  end

  % hopping for a time period of dtsec
  [p_r,p_s,p_t,p_n,ccc] = montecarlo(NP,p_r,p_s,p_t,p_n,Tsurf,dtsec,ccc,Qdissoc);
  
  n = n+1;
end

% write out final distribution
writeparticles(u51,NP,p_r,p_s,p_t,p_n)
writeglobe(u61,Tsurf)

cc = totalnrs(p_s);
disp(sprintf('# particles on surface %d',cc(1)))
disp(sprintf('# particles in flight %d',cc(2)))
disp(sprintf('# particles destroyed, photodissociation %d',ccc(1)))
disp(sprintf('# particles destroyed, escape %d',ccc(2)))
disp(sprintf('# particles coldtrapped %d %d %d',ccc(3),ccc(4),ccc(3)+ccc(4)))
disp(sprintf('# particles produced %d',0))
disp(sprintf('# active particles %d %d',sum(cc(1:6)),NP))

fclose(u30);
fclose(u50); fclose(u51);
fclose(u61);
