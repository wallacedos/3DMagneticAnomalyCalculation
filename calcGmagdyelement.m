function  Gelement = calcGmagdyelement(coeff,ulbound)
%  Summary of this function goes here.
%  it's the key function of calculating the total magnetic field anomaly  
%  with respect to the derivative of y.
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

r1=(x1^2 + y1^2 + z1^2)^(1/2);
r2=(x2^2 + y1^2 + z1^2)^(1/2);
r3=(x1^2 + y2^2 + z1^2)^(1/2);
r4=(x2^2 + y2^2 + z1^2)^(1/2);
r5=(x1^2 + y1^2 + z2^2)^(1/2);
r6=(x2^2 + y1^2 + z2^2)^(1/2);
r7=(x1^2 + y2^2 + z2^2)^(1/2);
r8=(x2^2 + y2^2 + z2^2)^(1/2);

c1=x1^2 + z1^2; c2=x2^2 + z1^2;
c3=x1^2 + z2^2; c4=x2^2 + z2^2;
c5=y1^2 + z1^2; c6=y1^2 + z2^2;
c7=y2^2 + z1^2; c8=y2^2 + z2^2;

Gelement1=k1*(y1/((x1 + r1)*r1) - y1/((x1 + r5)*r5) - ...
          y1/((x2 + r2)*r2) - y2/((x1 + r3)*r3) + ...
          y1/((x2 + r6)*r6) + y2/((x1 + r7)*r7) + ...
          y2/((x2 + r4)*r4) - y2/((x2 + r8)*r8)) ;

Gelement2=k2*((y1/r1 + 1)/(y1 + r1) - (y1/r5 + 1)/(y1 + r5) - ...
         (y1/r2 + 1)/(y1 + r2) + (y1/r6 + 1)/(y1 + r6) - ...
         (y2/r3 + 1)/(y2 + r3) + (y2/r7 + 1)/(y2 + r7) + ...
         (y2/r4 + 1)/(y2 + r4) - (y2/r8 + 1)/(y2 + r8));

Gelement3=k3*(y1/((z1 + r1)*r1) - y1/((z1 + r2)*r2) - ...
          y1/((z2 + r5)*r5) - y2/((z1 + r3)*r3) + ...
          y1/((z2 + r6)*r6) + y2/((z1 + r4)*r4) + ...
          y2/((z2 + r7)*r7) - y2/((z2 + r8)*r8));

Gelement4=k4*((x1/(z1*r1 + c1) - (x1*y1^2*z1)/((z1*r1 + c1)^2*r1))/((x1^2*y1^2)/(z1*r1 + c1)^2 + 1) - ...
         (x1/(z1*r3 + c1) - (x1*y2^2*z1)/((z1*r3 + c1)^2*r3))/((x1^2*y2^2)/(z1*r3 + c1)^2 + 1) - ...
         (x2/(z1*r2 + c2) - (x2*y1^2*z1)/((z1*r2 + c2)^2*r2))/((x2^2*y1^2)/(z1*r2 + c2)^2 + 1) - ...
         (x1/(z2*r5 + c3) - (x1*y1^2*z2)/((z2*r5 + c3)^2*r5))/((x1^2*y1^2)/(z2*r5 + c3)^2 + 1) + ...
         (x2/(z1*r4 + c2) - (x2*y2^2*z1)/((z1*r4 + c2)^2*r4))/((x2^2*y2^2)/(z1*r4 + c2)^2 + 1) + ...
         (x1/(z2*r7 + c3) - (x1*y2^2*z2)/((z2*r7 + c3)^2*r7))/((x1^2*y2^2)/(z2*r7 + c3)^2 + 1) + ...
         (x2/(z2*r6 + c4) - (x2*y1^2*z2)/((z2*r6 + c4)^2*r6))/((x2^2*y1^2)/(z2*r6 + c4)^2 + 1) - ...
         (x2/(z2*r8 + c4) - (x2*y2^2*z2)/((z2*r8 + c4)^2*r8))/((x2^2*y2^2)/(z2*r8 + c4)^2 + 1));

Gelement5=k5*((x1/(z1*r1 +c5) - (x1*y1*(2*y1 + (y1*z1)/r1))/(z1*r1 +c5)^2)/((x1^2*y1^2)/(z1*r1 +c5)^2 + 1) - ...
         (x2/(z1*r2 +c5) - (x2*y1*(2*y1 + (y1*z1)/r2))/(z1*r2 +c5)^2)/((x2^2*y1^2)/(z1*r2 +c5)^2 + 1) - ...
         (x1/(z2*r5 +c6) - (x1*y1*(2*y1 + (y1*z2)/r5))/(z2*r5 +c6)^2)/((x1^2*y1^2)/(z2*r5 +c6)^2 + 1) - ...
         (x1/(z1*r3 + c7) - (x1*y2*(2*y2 + (y2*z1)/r3))/(z1*r3 + c7)^2)/((x1^2*y2^2)/(z1*r3 + c7)^2 + 1) + ...
         (x2/(z2*r6 +c6) - (x2*y1*(2*y1 + (y1*z2)/r6))/(z2*r6 +c6)^2)/((x2^2*y1^2)/(z2*r6 +c6)^2 + 1) + ...
         (x2/(z1*r4 + c7) - (x2*y2*(2*y2 + (y2*z1)/r4))/(z1*r4 + c7)^2)/((x2^2*y2^2)/(z1*r4 + c7)^2 + 1) + ...
         (x1/(z2*r7 + c8) - (x1*y2*(2*y2 + (y2*z2)/r7))/(z2*r7 + c8)^2)/((x1^2*y2^2)/(z2*r7 + c8)^2 + 1) - ...
         (x2/(z2*r8 + c8) - (x2*y2*(2*y2 + (y2*z2)/r8))/(z2*r8 + c8)^2)/((x2^2*y2^2)/(z2*r8 + c8)^2 + 1));

Gelement6=k6*((x1/(z1*r1) - (x1*y1^2)/(z1*(x1^2 +c5)^(3/2)))/((x1^2*y1^2)/(z1^2*(x1^2 +c5)) + 1) - ...
          (x1/(z1*r3) - (x1*y2^2)/(z1*(x1^2 + c7)^(3/2)))/((x1^2*y2^2)/(z1^2*(x1^2 + c7)) + 1) - ...
          (x1/(z2*r5) - (x1*y1^2)/(z2*(x1^2 +c6)^(3/2)))/((x1^2*y1^2)/(z2^2*(x1^2 +c6)) + 1) - ...
          (x2/(z1*r2) - (x2*y1^2)/(z1*(x2^2 +c5)^(3/2)))/((x2^2*y1^2)/(z1^2*(x2^2 +c5)) + 1) + ...
          (x1/(z2*r7) - (x1*y2^2)/(z2*(x1^2 + c8)^(3/2)))/((x1^2*y2^2)/(z2^2*(x1^2 + c8)) + 1) + ...
          (x2/(z1*r4) - (x2*y2^2)/(z1*(x2^2 + c7)^(3/2)))/((x2^2*y2^2)/(z1^2*(x2^2 + c7)) + 1) + ...
          (x2/(z2*r6) - (x2*y1^2)/(z2*(x2^2 +c6)^(3/2)))/((x2^2*y1^2)/(z2^2*(x2^2 +c6)) + 1) - ...
          (x2/(z2*r8) - (x2*y2^2)/(z2*(x2^2 + c8)^(3/2)))/((x2^2*y2^2)/(z2^2*(x2^2 + c8)) + 1));

Gelement=Gelement1+Gelement2+Gelement3+Gelement4+Gelement5+Gelement6;
end

