# 3DMagneticAnomalyCalculation
 Matlab program for the calculation of 3D total magnetic field anomaly and its derivative.<br>
# 3DMagneticAnomalyCalculation程序包使用说明<br>
## 简介<br>
　　基于matlab的程序包3DMagneticAnomalyCalculation可实现对三维地质体在地表产生的总磁异常及其*x*、*y*、*z*方向的导数的计算，三维地质体的剖分单元为长方体。程序包路径下的main.m给出了单一长方体磁异常计算的例子，以验证程序的正确性。程序包下载完成后，将其放在matlab的当前路径下，运行main.m即可获得计算结果。请注意，main.m当中的注释包含了对各函数的功能说明，一定要仔细阅读。使用者可以通过更改main.m中的地质模型参数，轻松完成对其他三维地质模型的总磁异常和导数的计算。<br>
## 算例：单一长方体磁异常数值模拟<br>
　　长发体参数：<br>
　　X方向：6000-14000 m；Y方向：8000-12000 m；Z方向：1000-5000 m；地磁场倾角50度，偏角30度；磁化方向倾角60度，偏角10度。<br>
　　观测系统参数：<br>
　　观测面z=0 m；X方向：0-20000 m;方向：0-20000 m;点距和线距均为100 m。以下4幅图分别为总磁异常及其*x*、*y*、*z*方向导数的等值线图。<br>
  
  ![图 1 单一长方体总磁异常等值线图](https://github.com/wallacedos/3DMagneticAnomalyCalculation/blob/main/Figures/Fig1.jpg) <br>
  
  ![图 2 单一长方体总磁异常x方向导数等值线图](https://github.com/wallacedos/3DMagneticAnomalyCalculation/blob/main/Figures/Fig2.jpg) <br>
  
  ![图 3 单一长方体总磁异常y方向导数等值线图](https://github.com/wallacedos/3DMagneticAnomalyCalculation/blob/main/Figures/Fig3.jpg) <br>
  
  ![图 4 单一长方体总磁异常z方向导数等值线图](https://github.com/wallacedos/3DMagneticAnomalyCalculation/blob/main/Figures/Fig4.jpg) <br>
 
 ## 参考文献 <br>
郭志宏,管志宁,熊盛青.长方体ΔT场及其梯度场无解析奇点理论表达式[J].地球物理学报,2004(06):1131-1138. <br>
## 更多精彩内容
![](https://github.com/wallacedos/3DMagneticAnomalyCalculation/blob/main/Figures/mmexport1617252557610.jpg)
