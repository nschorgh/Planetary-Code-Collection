
a=load('qmean.topo81');


Nx=81; Ny=81;

Ny=Ny-2; Nx=Nx-2;

h=reshape(a(:,3),Ny,Nx);
Qm=reshape(a(:,6),Ny,Nx);
Qx=reshape(a(:,11),Ny,Nx);
T=reshape(a(:,12),Ny,Nx);

hlevs=[-70:10:0];
hmin=min(h(:));
hmax=max(h(:));

figure(1); clf
set(0,'defaultaxesfontsize',11,'defaulttextfontsize',9)
set(0,'defaultlinelinewidth',2,'defaultaxesfontsize',9)

subplot(2,2,1)
imagesc(h)
shading flat
axis equal ij
hold on
[c,hcont]=contour(h,hlevs,'-k');
set(hcont,'linewidth',2)
hold off
barh1=colorbar;
set(get(barh1,'ylabel'),'string','Elevation (m)','fontsize',11)
axis off


subplot(2,2,2)
imagesc(Qm);
shading flat
axis equal ij
barh=colorbar;
set(get(barh,'ylabel'),'string','Mean Absorption (W/m^2)','fontsize',11)
axis off
Qmin=min(Qm(:));
Qmax=max(Qm(:));
hold on
r=(Qmax-Qmin)/(hmax-hmin);
[c,hcont]=contour(Qmin+(h-hmin)*r,Qmin+(hlevs-hmin)*r,'-k');
set(hcont,'linewidth',2)
hold off



subplot(2,2,3)
imagesc(Qx);
shading flat
axis equal ij
barh=colorbar;
set(get(barh,'ylabel'),'string','Peak Absorption (W/m^2)','fontsize',11)
axis off
Qmin=min(Qx(:));
Qmax=max(Qx(:));
hold on
r=(Qmax-Qmin)/(hmax-hmin);
[c,hcont]=contour(Qmin+(h-hmin)*r,Qmin+(hlevs-hmin)*r,'-k');
set(hcont,'linewidth',2)
hold off


subplot(2,2,4)
%T = (Qx/(0.95*5.67e-8)).^0.25;
imagesc(T);
shading flat
axis equal ij
barh=colorbar;
set(get(barh,'ylabel'),'string','Peak Temperature (K)','fontsize',11)
axis off
Qmin=min(T(:));
Qmax=max(T(:));
hold on
r=(Qmax-Qmin)/(hmax-hmin);
[c,hcont]=contour(Qmin+(h-hmin)*r,Qmin+(hlevs-hmin)*r,'-k');
set(hcont,'linewidth',2)
hold off


set(gcf,'inverthardcopy','off','color','w')
print -dpng -r720 qmean.png


