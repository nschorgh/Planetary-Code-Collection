# Octave script

clear all
a=load('qmean.topo81');
b=load('qmean.topo81_wsurfcond');

Nx=81-2; Ny=81-2; 

h=reshape(a(:,3),Ny,Nx);
Q1=reshape(a(:,6),Ny,Nx);
T1=reshape(a(:,12),Ny,Nx);
Q2=reshape(b(:,6),Ny,Nx);
T2=reshape(b(:,12),Ny,Nx);

figure(1); clf
set(gcf,'defaultlinelinewidth',2,'defaultaxesfontsize',18,'defaulttextfontsize',18)

T=[50:5:400];
ah1=hist(T1(:),T);
ah2=hist(T2(:),T);
hold on
plot(T,ah1,'r-',T,ah2,'g--');
hold off
xlabel('Peak Temperature (K)')
ylabel('# Pixels')
box on

legend('Equilibrium I=0','with subsurface I=100','location','northwest')
title('Temperature Histograms')

print -depsc qmean_compare.eps



