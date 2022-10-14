function [ score ] = breathingDataScore_slope( dataIn )
envelope_points = 300;
[up, down] = envelope(dataIn,envelope_points,'peak');  % 求包络线的方式1
n = 0:1:size(dataIn,2)-1;
a1 = polyfit(n,down,1); % 1阶多项式最小二乘拟合
y1 = polyval(a1,n);
score = abs( 1 / ( (n(1,2) - n(1,1)) / (y1(1,2) - y1(1,1)) ) );
end

