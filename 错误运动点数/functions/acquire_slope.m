function [slope] = acquire_slope(dataIn, num)
%% 求斜率
% 对数据进行归一化
% data_dealed = mapminmax(dataIn,0,size(dataIn,2));
% n_samples = size(data_dealed,2);
% t = 0:1:n_samples-1;
% %
% coef_fit = polyfit(t(1,:),data_dealed(1,:),1);
% slope = abs(coef_fit(1,1));

% 对数据进行归一化
data_dealed = mapminmax(dataIn,0,1);
num_epoch = num;
len_epoch = floor(size(data_dealed,2) / num_epoch); % 分成四份

part1 = max(data_dealed(1,1:len_epoch));
part2 = max(data_dealed(1,end-len_epoch+1:end));

fenzi = abs(part1 - part2);
fenmu = min(part1, part2);
slope = fenzi / fenmu;

end

