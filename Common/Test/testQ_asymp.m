% plot output of testcrankQ_asymp.f90 with Octave

a=load('./Tsurface');
%a=load('/home/norbert/Thermal/ThermalLibF/Common/Test/Tsurface');
nz=size(a,2)-1;
t=a(:,1);  % time
Tsurf=a(:,2);



% compare with short-term analytical solution for sudden temperature drop

thIn = 500.;
T0 = 200.;
emiss = 1.;
sigSB=5.6704e-8;

k = find(t>=0);
tanaly = [0:0.2:max(t)];
Tanaly = T0 + 2/sqrt(pi)*emiss*sigSB/thIn*(-T0**4)*sqrt(tanaly);

clf;
set(gcf,'defaultaxesfontsize',14)

hold on
plot(tanaly,Tanaly,'k-','linewidth',2)
plot(t,Tsurf,'ro')
hold off

hl = legend('Small time expansion','Numerical solution','location','southwest');
set(hl,'fontsize',14)

xlabel('Time (s)')
ylabel('Temperature (K)')
box on

%print -depsc testQ_asymp.eps
