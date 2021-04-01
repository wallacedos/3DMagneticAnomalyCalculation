function Gelement= calcGmagelement(coeff,ulbound)
%  Summary of this function goes here.
%  it's the key function of calculating the total magnetic field anomaly.
%  More details about the formula you can find from the follwing
%  references.
%
%  References: 
%  郭志宏,管志宁,熊盛青.长方体ΔT场及其梯度场无解析奇点理论表达式
%  [J].地球物理学报,2004(06):1131-1138.
%
%  Author(s): Yan Yingwei
%  Copyright: 2019-2022 
%  Revision: 1.0  Date: 3/16/2019
%
%  Department of Geophysics, Jilin University.

k1=coeff(1);k2=coeff(2);k3=coeff(3);k4=coeff(4);k5=coeff(5);k6=coeff(6);
x1=ulbound(1);x2=ulbound(2);y1=ulbound(3);y2=ulbound(4);z1=ulbound(5);z2=ulbound(6);

% %The calculation speed could be very slow if you adopt the following method.
% X=[x1 x2];Y=[y1 y2];Z=[z1 z2];
% 
% Gelement1=0;
% Gelement2=0;
% Gelement3=0;
% Gelement4=0;
% Gelement5=0;
% Gelement6=0;
% 
% for i=1:2
%     for j=1:2
%         for k=1:2
%             n=mod(i,2)+mod(j,2)+mod(k,2);
%             n=mod(n,2)+1;
%             r=sqrt(X(i)*X(i)+Y(j)*Y(j)+Z(k)*Z(k));
%             
%             Gelement1=Gelement1+(-1)^(n)*log(X(i)+r);
%             Gelement2=Gelement2+(-1)^(n)*log(Y(j)+r);
%             Gelement3=Gelement3+(-1)^(n)*log(Z(k)+r);
%             Gelement4=Gelement4+(-1)^(n)*atan(X(i)*Y(j)/(X(i)^2+r*Z(k)+Z(k)^2));
%             Gelement5=Gelement5+(-1)^(n)*atan(X(i)*Y(j)/(Y(j)^2+r*Z(k)+Z(k)^2));
%             Gelement6=Gelement6+(-1)^(n)*atan(X(i)*Y(j)/(r*Z(k)));
%         end
%     end
% end

r1=(x1^2 + y1^2 + z1^2)^(0.5);
r2=(x2^2 + y1^2 + z1^2)^(0.5);
r3=(x1^2 + y2^2 + z1^2)^(0.5);
r4=(x2^2 + y2^2 + z1^2)^(0.5);
r5=(x1^2 + y1^2 + z2^2)^(0.5);
r6=(x2^2 + y1^2 + z2^2)^(0.5);
r7=(x1^2 + y2^2 + z2^2)^(0.5);
r8=(x2^2 + y2^2 + z2^2)^(0.5);

Gelement1=k1*(log(x1 + r1)-log(x2 + r2)... 
- log(x1 + r3) + log(x2 + r4)... 
- log(x1 + r5) + log(x2 + r6)... 
+ log(x1 + r7) - log(x2 + r8)); 

Gelement2=k2*(log(y1 + r1)-log(y1 + r5)... 
- log(y1 + r2)+log(y1 + r6)... 
- log(y2 + r3)+log(y2 + r7)... 
+ log(y2 + r4)-log(y2 + r8));

Gelement3=k3*(log(z1 + r1)-log(z1 + r3)... 
- log(z1 + r2)+log(z1 + r4)... 
- log(z2 + r5)+log(z2 + r7)... 
+ log(z2 + r6)-log(z2 + r8));

Gelement4=k4*(atan((x1*y1)/(z1*r1 + x1^2 + z1^2))-atan((x1*y2)/(z1*r3 + x1^2 + z1^2))... 
- atan((x1*y1)/(z2*r5 + x1^2 + z2^2))-atan((x2*y1)/(z1*r2 + x2^2 + z1^2))... 
+ atan((x1*y2)/(z2*r7 + x1^2 + z2^2))+atan((x2*y2)/(z1*r4 + x2^2 + z1^2))... 
+ atan((x2*y1)/(z2*r6 + x2^2 + z2^2))-atan((x2*y2)/(z2*r8 + x2^2 + z2^2)));

Gelement5=k5*(atan((x1*y1)/(z1*r1 + y1^2 + z1^2))-atan((x2*y1)/(z1*r2 + y1^2 + z1^2))... 
- atan((x1*y1)/(z2*r5 + y1^2 + z2^2))-atan((x1*y2)/(z1*r3 + y2^2 + z1^2))... 
+ atan((x2*y1)/(z2*r6 + y1^2 + z2^2))+atan((x2*y2)/(z1*r4 + y2^2 + z1^2))... 
+ atan((x1*y2)/(z2*r7 + y2^2 + z2^2))-atan((x2*y2)/(z2*r8 + y2^2 + z2^2)));

Gelement6=k6*(atan((x1*y1)/(z1*r1))-atan((x1*y1)/(z2*r5))... 
- atan((x1*y2)/(z1*r3))-atan((x2*y1)/(z1*r2))... 
+ atan((x1*y2)/(z2*r7))+atan((x2*y1)/(z2*r6))... 
+ atan((x2*y2)/(z1*r4))-atan((x2*y2)/(z2*r8)));
 
Gelement=Gelement1+Gelement2+Gelement3+Gelement4+Gelement5+Gelement6;
end

