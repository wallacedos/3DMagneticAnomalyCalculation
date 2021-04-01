function deltaT = calc3Dmaganomaly(Grid,Mpt,I,I0,A,A0)
%   Summary of this function goes here.
%   deltaT = calc3Dmaganomaly(Grid,Mpt,I,I0,A,A0)
%   Detailed explanation goes here.
%   The function is for calculating the total magnetic field anomaly.
%
%   IN   Grid: the grid information of geological body (it's a struct).
%        Mpt:  the coordinates of measured points (it's a struct).
%        I:    the inclination of geomagnetic field.
%        I0:   the inclination of magnetization.
%        A:    the deflection angle of geomagnetic field.
%        A0:   the deflection angle of magnetization.
%
%  OUT  deltaT: the total  magnetic field anomaly.
%
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2019-2022 
%  Revision: 1.0  Date: 3/16/2019
%
%  Department of Geophysics, Jilin University.

XGrid=Grid.X;
YGrid=Grid.Y;
ZGrid=Grid.Z;
MagI=Grid.MagI;  % magnetization

[m,n,p]=size(MagI);


XMpt=Mpt.X;
YMpt=Mpt.Y;
XMptvec=XMpt(:);
YMptvec=YMpt(:);
[m1,n1]=size(XMpt);

M=m*n*p;
N=m1*n1;

c1=cos(I0)*sin(A0)*sin(I)+sin(I0)*cos(I)*sin(A);
c2=cos(I0)*cos(A0)*sin(I)+sin(I0)*cos(I)*cos(A);
c3=cos(I0)*cos(A0)*cos(I)*sin(A)+cos(I0)*sin(A0)*cos(I)*cos(A);
c4=cos(I0)*cos(A0)*cos(I)*cos(A);
c5=cos(I0)*sin(A0)*cos(I)*sin(A);
c6=-sin(I0)*sin(I);

coeff=[c1; c2; c3; c4; c5; c6];
mu0=4*pi*10^(-7);

deltaT=zeros(N,1);
ulbound=zeros(6,1);
h = waitbar(0,'Please wait...');
for i=1:N
    waitbar(i/N,h)
    for j=1:M
        kind=floor(j/(m*n));
        colind=floor((j-kind*m*n)/m);
        rowind=j-kind*m*n-colind*m;
        
        if((j-kind*m*n)==0)
            colind=n;
            rowind=m;
        elseif((j-kind*m*n-colind*m)==0&&(j-kind*m*n)~=0)
            rowind=m;
            kind=kind+1;
        else
            kind=kind+1;
            colind=colind+1;
        end
    
        ulbound(1)=XGrid(rowind+1,colind+1,kind+1)-XMptvec(i);
        ulbound(2)=XGrid(rowind,colind,kind)-XMptvec(i);
        ulbound(3)=YGrid(rowind+1,colind+1,kind+1)-YMptvec(i);
        ulbound(4)=YGrid(rowind,colind,kind)-YMptvec(i);
        ulbound(5)=ZGrid(rowind+1,colind+1,kind+1);
        ulbound(6)=ZGrid(rowind,colind,kind);
        
        deltaT_sub=mu0/(4*pi)*MagI(rowind,colind,kind)*calcGmagelement(coeff,ulbound);
        deltaT(i)=deltaT(i)+deltaT_sub;
    end
end
deltaT=reshape(deltaT,m1,n1)*10^(9);
end

