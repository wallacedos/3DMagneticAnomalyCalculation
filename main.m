%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
%  Main program for the calculation of 3D total magnetic field anomaly and its derivative.
%  Author(s): Yan Yingwei
%  Copyright: 2019-2022 
%  Revision: 2.0  Date: 4/16/2019
% 
% This is the main program.
% It involves function 
%  calc3Dmaganomaly.m:    the function can read the information from main.m and stack the 
%                         total magnetic field anomaly of each block.
%  calcGmagelement.m:     the key function of calculating the total magnetic field anomaly.
%
%  calc3Dmaganomalydx.m:  the function can read the information from main.m and stack the 
%                         total magnetic field anomaly with respect to the derivative of x 
%                         of each block.
%  calcGmagdxelement.m:   the key function of calculating the total magnetic field anomaly 
%                         with respect to the derivative of x.
%
%  calc3Dmaganomalydy.m:  the function can read the information from main.m and stack the 
%                         total magnetic field anomaly with respect to the derivative of y 
%                         of each block.
%  calcGmagdyelement.m:   the key function of calculating the total magnetic field anomaly 
%                         with respect to the derivative of y.
%
%  calc3Dmaganomalydz.m:  the function can read the information from main.m and stack the 
%                         total magnetic field anomaly with respect to the derivative of z 
%                         of each block.
%  calcGmagdzelement.m:   the key function of calculating the total magnetic field anomaly 
%                         with respect to the derivative of z.
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

%% Demo1
% tic
% [X,Y,Z]=meshgrid(0:50:1000,0:50:1000,0:50:400);  % the coordinates of geological body.
% 
% Grid.X=X;
% Grid.Y=Y;
% Grid.Z=Z;
% 
% [m,n,p]=size(X);
% MagI=zeros(m-1,n-1,p-1);  % magnetization
% MagI(9,9,4)=1;
% 
% 
% Grid.MagI=MagI;
% 
% [Xpt,Ypt]=meshgrid(-1:10:1009,-1:10:1009); % the coordinates of measured points
% Mpt.X=Xpt;
% Mpt.Y=Ypt;
% 
% I=pi/4;  %  the inclination of geomagnetic field.
% A=pi/4;  %  the deflection angle of geomagnetic field.
% I0=pi/4; %  the inclination of magnetization.
% A0=pi/4; %  the deflection angle of magnetization.
% 
% % I=pi/2;  %  the inclination of geomagnetic field.
% % A=0;  %  the deflection angle of geomagnetic field.
% % I0=pi/2; %  the inclination of magnetization.
% % A0=0; %  the deflection angle of magnetization.
% 
% deltaT = calc3Dmaganomaly(Grid,Mpt,I,I0,A,A0);  % the total magnetic field anomaly.
% 
% toc
% 
% x=Xpt(1,:);
% y=Ypt(:,1);
% 
% figure(1)
% contourf(x,y,deltaT,20);
% colorbar;
% colormap(jet);
% set(gca,'YDir','reverse');
% set(gca,'XAxisLocation','top');
% xlabel('Distance (m)');
% ylabel('Distance (m)');
% xlabel(colorbar,'\DeltaT (nT)');


%% Demo2 Single cuboid
tic
[X,Y,Z]=meshgrid(6000:8000:14000,8000:4000:12000,1000:4000:5000);  % the coordinates of geological body.

Grid.X=X;
Grid.Y=Y;
Grid.Z=Z;

[m,n,p]=size(X);
MagI=zeros(m-1,n-1,p-1);  % magnetization
MagI(1,1,1)=1;


Grid.MagI=MagI;

[Xpt,Ypt]=meshgrid(0:100:20000,0:100:20000); % the coordinates of measured points
Mpt.X=Xpt;
Mpt.Y=Ypt;

I=50*pi/180;  %  the inclination of geomagnetic field.
A=30*pi/180;  %  the deflection angle of geomagnetic field.
I0=60*pi/180; %  the inclination of magnetization.
A0=10*pi/180; %  the deflection angle of magnetization.

deltaT = calc3Dmaganomaly(Grid,Mpt,I,I0,A,A0);  % the total magnetic field anomaly.

toc

x=Xpt(1,:);
y=Ypt(:,1);

figure(1)
contour(x,y,deltaT,80);
colorbar;
colormap(jet);
title('单一长方体总磁异常等值线图');
xlabel('Distance (m)');
ylabel('Distance (m)');
xlabel(colorbar,'\DeltaT (nT)');

tic
deltaTdx = calc3Dmaganomalydx(Grid,Mpt,I,I0,A,A0); 
toc

figure(2)
contour(x,y,deltaTdx,80);
colorbar;
colormap(jet);
title('单一长方体总磁异常x方向导数等值线图');
xlabel('Distance (m)');
ylabel('Distance (m)');
xlabel(colorbar,'nT/m');

tic
deltaTdy = calc3Dmaganomalydy(Grid,Mpt,I,I0,A,A0); 
toc

figure(3)
contour(x,y,deltaTdy,80);
colorbar;
colormap(jet);
title('单一长方体总磁异常y方向导数等值线图');
xlabel('Distance (m)');
ylabel('Distance (m)');
xlabel(colorbar,'nT/m');

tic
deltaTdz = calc3Dmaganomalydz(Grid,Mpt,I,I0,A,A0); 
toc

figure(4)
contour(x,y,deltaTdz,80);
colorbar;
colormap(jet);
title('单一长方体总磁异常z方向导数等值线图');
xlabel('Distance (m)');
ylabel('Distance (m)');
xlabel(colorbar,'nT/m');
