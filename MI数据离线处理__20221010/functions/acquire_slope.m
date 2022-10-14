function [slope] = acquire_slope(dataIn)
%% 求斜率

% 对数据进行归一化
data_dealed = mapminmax(dataIn,0,1);
len_epoch = floor(size(data_dealed,2) / 4); %
part1 = max(data_dealed(1,1:len_epoch));
part2 = max(data_dealed(1,end-len_epoch+1:end));
slope_part1 = abs(part1 - part2) / min(part1, part2);
%
data_dealed = mapminmax(-1 * dataIn,0,1);
len_epoch = floor(size(data_dealed,2) / 4); %
part1 = max(data_dealed(1,1:len_epoch));
part2 = max(data_dealed(1,end-len_epoch+1:end));
slope_part2 = abs(part1 - part2) / min(part1, part2);
%
slope = min(slope_part1,slope_part2);
end

