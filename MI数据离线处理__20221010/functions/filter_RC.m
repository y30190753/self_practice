function dataOut = filter_RC(dataIn)
% FILTER_RC RCµÍÍ¨
% RCµÍÍ¨
fs = 50;
fc = 0.1;
wc = 2*pi*fc;
Ts = 1/fs;
b = [1/wc -1/wc];
a = [1/wc+Ts -1/wc];

dataNum = length(dataIn);
xPre = 0;
yPre = 0;
dataOut = zeros(1,dataNum);

for idx = 1:dataNum
    in = dataIn(idx);
    out = (single(b(1))*single(in) + single(b(2))*single(xPre) - single(a(2))*single(yPre))/single(a(1));
    xPre = in;
    yPre = out;
    
    dataOut(idx) = yPre;
end
end

