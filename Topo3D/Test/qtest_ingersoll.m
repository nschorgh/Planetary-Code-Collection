
% bowl test

a=load('qinst.topo81');
Nx=81-2; Ny=81-2;

x=reshape(a(:,1),Ny,Nx);  x=x(1:Nx:end);
y=reshape(a(:,2),Ny,Nx);  y=y(1:Ny);
h=reshape(a(:,3),Ny,Nx);
Qdirect=reshape(a(:,5),Ny,Nx);
Qirre=reshape(a(:,6),Ny,Nx);
Qabs=reshape(a(:,7),Ny,Nx);
Qir=reshape(a(:,8),Ny,Nx);
Qrefl=reshape(a(:,9),Ny,Nx);
Tsurf=reshape(a(:,10),Ny,Nx);

set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultlinelinewidth',2)

y=ceil((Ny+2)/2)-1;

figure(1); clf
hold on
plot(Qdirect(:,y)/100,'k-')
plot(Qir(:,y),'ro-')
plot(Qrefl(:,y),'b-^')
hold off
xlabel('Distance (Pixel)')
ylabel('Flux (W/m^2)')
box on
%set(gca,'yscale','log')


figure(2); clf
plot(Tsurf(:,y),'go-','markersize',3)
xlabel('Distance (Pixel)')
ylabel('Temperature (K)')

S0=1365;
A=0.12;
emiss=0.95;
beta=10*pi/180;  % elevation of sun
D2d=5;
% analytical solution from Ingersoll et al. (1992), Icarus 100, 40-47.
f=1/(1+D2d^2/4);
Tingersoll = (S0*sin(beta)*f*(1-A)/(1-A*f)*(1+A*(1-f)/emiss)/5.67e-8)^0.25;
hold on
k=find(Tsurf(:,y)<200);
plot([min(k) max(k)],Tingersoll*[1 1],'r-','linewidth',2)
hold off
legend('Numerical','Analytical')

%print -depsc qtest_ingersoll.eps
