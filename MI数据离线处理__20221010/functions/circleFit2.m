function [center,radius] = circleFit2(dataIn)
%CIRCLEFIT2 此处显示有关此函数的摘要
%   此处显示详细说明
xlist = real(dataIn);
ylist = imag(dataIn);
n = length(dataIn);

Mx = sum(xlist);
My = sum(ylist);
Mxx = sum(xlist.^2);
Myy = sum(ylist.^2);
Mxy = sum(xlist.*ylist);
Mxxx = sum(xlist.*xlist.*xlist);
Myyy = sum(ylist.*ylist.*ylist);
Mxxy = sum(xlist.*xlist.*ylist);
Mxyy = sum(xlist.*ylist.*ylist);

M1 = n*Mxy - Mx*My;
M2 = n*Mxx - Mx*Mx;
M3 = n*Mxxx + n*Mxyy - Mx*(Mxx+Myy);
M4 = n*Myy - My*My;
M5 = n*Myyy + n*Mxxy - My*(Mxx+Myy);

a = (M1*M5 - M3*M4)/(M2*M4 - M1*M1);
b = (M1*M3 - M2*M5)/(M2*M4 - M1*M1);
c = -(a*Mx + b*My + Mxx + Myy)/n;

center = -a/2 + 1i*(-b/2);
radius = ((a*a + b*b - 4*c)^0.5)/2;
end

