% Compare numerical solution for bowl-shaped crater with analytical solution

a=load('qinst.topo81');
Nx=81-2; Ny=81-2;

x=reshape(a(:,1),Ny,Nx);  x=x(1:Nx:end);
y=reshape(a(:,2),Ny,Nx);  y=y(1:Ny);
h=reshape(a(:,3),Ny,Nx);
Qdirect=reshape(a(:,5),Ny,Nx);
Qir=reshape(a(:,6),Ny,Nx);
Qrefl=reshape(a(:,7),Ny,Nx);
Qabs=reshape(a(:,8),Ny,Nx);
Tsurf=reshape(a(:,9),Ny,Nx);

set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12,'defaultlinelinewidth',2)

y=ceil((Ny+2)/2)-1;


S0=1365;
A=0.12;
emiss=0.95;
beta=10*pi/180;  % elevation of sun
D2d=5;


figure(1); clf
hold on
plot(Qdirect(:,y)*(1-A),'k-')
plot(Qir(:,y),'ro-')
plot(Qrefl(:,y),'b-^')
Qtotal=(1-A)*(Qdirect+Qrefl)+emiss*Qir;
plot(Qtotal(:,y),'go-','markersize',3)
%plot(sinbeta(:,y)*S0*(1-A),'--','color',[.4 .4 .4])
hold off
xlabel('Distance (Pixel)')
ylabel('Flux (W/m^2)')
box on
%set(gca,'yscale','log')
legend('direct absorbed','IR','scattered vis','total absorbed')



figure(2); clf
plot(Tsurf(:,y),'go-','markersize',3)
hold on
%plot(Tpred(:,y),'--*','color',[.4 .4 .4])
hold off
xlabel('Distance (Pixel)')
ylabel('Temperature (K)')


% analytical solution from Ingersoll et al. (1992), Icarus 100, 40-47.
f=1/(1+D2d^2/4);
Tingersoll = (S0*sin(beta)*f*(1-A)/(1-A*f)*(1+A*(1-f)/emiss)/5.67e-8)^0.25;
hold on
k=find(Tsurf(:,y)<200);
plot([min(k) max(k)],Tingersoll*[1 1],'r-','linewidth',2)
hold off
legend('Numerical','Analytical')

%print -depsc qtest_ingersoll.eps




figure(3); clf
colormap jet;
imagesc(Tsurf)
axis equal off
barh=colorbar;
set(get(barh,'ylabel'),'string','Temperature (K)')


Toutside = ( (1-A)*S0*sin(beta)/emiss/5.67e-8 )^0.25;

dx=5; dy=5;
[sy,sx] = gradient(h);
sx=sx/dx;
sy=sy/dy;
% get surface normals
snorm = sqrt(1 + sx.^2 + sy.^2);
sx = sx./snorm;
sy = sy./snorm;
sz = ones(size(sx))./snorm;
% calculate sin(local sun elevation)
sinbeta = -sx*cos(beta) + sy*0 + sz*sin(beta);
sinbeta(sinbeta<0)=0;

% construct temperatures for the entire surface
b = f/(1-A*f)*(emiss+A*(1-f));
Epred = zeros(size(Qabs));
k = find(Qdirect>0); % areas of direct sunlight
Epred(k) = sinbeta(k);
k = find(Qir>0); % inside crater (terrain irradiance)
Epred(k) = Epred(k) + b*sin(beta);
Tpred = ( (1-A)*S0*Epred/emiss/5.67e-8 ).^0.25;



figure(4); clf
colormap rainbow;
%Tsurf(Tsurf>200)=NaN;
%err = Tsurf-Tingersoll;
err = Tsurf - Tpred;
imagesc(err);
set(gca,'clim',[prctile(err(:),5) prctile(err(:),95)])
axis equal off
barh=colorbar;
set(get(barh,'ylabel'),'string','Temperature Difference (K)')


sum(abs(err(:)))/length(err(:))


