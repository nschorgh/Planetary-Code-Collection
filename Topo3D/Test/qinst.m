clear all
a=load('q.dat');
%a=load('qinst.topo81');
%a=load('qinst.dat');

Nx=81-2; Ny=81-2; 

h=reshape(a(:,3),Ny,Nx);
Q=reshape(a(:,5),Ny,Nx);
Qabs=reshape(a(:,7),Ny,Nx);
Qir=reshape(a(:,8),Ny,Nx);
Qre=reshape(a(:,9),Ny,Nx);
T=reshape(a(:,10),Ny,Nx);

hlevs=[-50:10:0];
hmin=min(h(:));
hmax=max(h(:));

figure(1); clf
set(gcf,'defaultaxesfontsize',11,'defaulttextfontsize',10)
set(gcf,'defaultlinelinewidth',2,'defaultaxesfontsize',10)
colormap('jet')


subplot(2,2,1)
imagesc(h);
shading flat
axis equal ij
barh=colorbar;
set(get(barh,'ylabel'),'string','Elevation (m)')
axis off
hold on
[c,hcont]=contour(h,hlevs,'k-');
set(hcont,'linewidth',2)
hold off


subplot(2,2,2)
imagesc(Q);
shading flat
axis equal ij
barh=colorbar;
set(get(barh,'ylabel'),'string','Direct Insolation (W/m^2)')
axis off
Qmin=min(Q(:));
Qmax=max(Q(:));
hold on
r=(Qmax-Qmin)/(hmax-hmin);
[c,hcont]=contour(Qmin+(h-hmin)*r,Qmin+(hlevs-hmin)*r,'k-');
set(hcont,'linewidth',2)
hold off


subplot(2,2,3)
imagesc(Qir);
shading flat
axis equal ij
barh=colorbar;
set(get(barh,'ylabel'),'string','Q_{IR} (W/m^2)')
axis off
Qmin=min(Qir(:));
Qmax=max(Qir(:));
hold on
r=(Qmax-Qmin)/(hmax-hmin);
[c,hcont]=contour(Qmin+(h-hmin)*r,Qmin+(hlevs-hmin)*r,'k-');
set(hcont,'linewidth',2)
hold off


subplot(2,2,4)
imagesc(T);
shading flat
axis equal ij
barh=colorbar;
set(get(barh,'ylabel'),'string','Temperature (K)')
axis off
Qmin=min(T(:));
Qmax=max(T(:));
hold on
r=(Qmax-Qmin)/(hmax-hmin);
[c,hcont]=contour(Qmin+(h-hmin)*r,Qmin+(hlevs-hmin)*r,'k-');
set(hcont,'linewidth',2)
hold off



set(gcf,'inverthardcopy','off','color','w')
print -dpng -r720 qinst.png


