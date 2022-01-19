%***********************************************************************
% test Crank Nicolson subroutine 
%***********************************************************************

Period = 88775.244*670   % [seconds]
Fgeo = 0.; Ta=30.; Tm=190.;
nz = 70;
  
T = zeros(nz,1);
rhocv = zeros(nz,1); % volumetric heat capacity
ti = zeros(nz,1);
z = zeros(nz,1);

NSTEPS = 50000;
STEPSPERSOL = 120;
dt = Period/STEPSPERSOL; % Time step [seconds]
zmax = 2.5
zfac=1.02;
thIn = 120. % thermal inertia
  
rhocv(:) = 1200.*800.;
delta = thIn/rhocv(1)*sqrt(Period/pi)  % skin depth

fout = fopen('Tprofile','w'); % temperature profiles

T(:) = 0.;
ti(:) = thIn;
  
z = setgrid(nz,zmax,zfac);
f2=fopen('z','w');
fprintf(f2,'%f ',z(:))
fprintf(f2,'\n')
fclose(f2);
  
time = 0.;
  
Tsurf = Tm + Ta*sin(2*pi*time/Period);

for n = 0:NSTEPS
  time = (n+1)*dt;     %  time at n+1
  Tsurfp1 = Tm + Ta*sin(2*pi*time/Period);
  [T, Fsurf] = conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhocv,Fgeo);
  Tsurf = Tsurfp1;
     
  % write 12 profiles from the last sol
  if n > NSTEPS-STEPSPERSOL,
    if mod(n,10) == 9,  % synched with test_Tprofile.m
      disp([time/Period, Tsurf])
      fprintf(fout,[repmat('%7.2f ',1,nz+1),'\n'],Tsurf,T(1:nz))
    end
  end
     
end         % end the time loop

fclose(fout);

