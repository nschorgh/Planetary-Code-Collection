% read output of testcrankT.f90
wpth='./';

a=load([wpth,'Tprofile_testcrankT']);
nz=size(a,2)-1;
t=a(:,1);  % time
T=a(:,2:end);

fid=fopen([wpth,'z_testcrankT'],'r');
z=fscanf(fid,'%f');
fclose(fid);
%z=[0;z];

% plot temperature profiles
clf;
set(gcf,'defaultaxesfontsize',12)  % Octave

h1=plot(T(:,:),z,'k-');
axis ij
xlabel('Temperature (K)')
ylabel('Depth, z (m)')


% compare with analytical solution for sinusoidal surface temperature
%                                oscillation and semi-infinite domain
Ta=30.; Tm=190.; P=670.*88775.244;
rhoc = 1200*800; thIn=120;
delta=thIn/rhoc*sqrt(P/pi);   % skin depth
w=2*pi/P;
hold on
dt=P/12;
for t=0:dt:P
  T=Tm+Ta*exp(-z/delta).*sin(z/delta-w*t);
  h2=plot(T,z,'r--');
end
hold off

hl = legend([h1(1), h2],'Numerical','Analytical','location','southwest');
set(hl,'fontsize',12)

% print -dpng test_Tprofile.png
% print -depsc test_Tprofile.eps
