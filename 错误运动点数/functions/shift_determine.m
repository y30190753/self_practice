function [shift_flag, slope_value] = shift_determine(dataIn, num_epoch, threshold_max, threshold_min)
%%  判断一个数据块是否出现基线漂移
%      输入 dataIn, 1 * n_points
%      输出 flag == 1 判断有基线漂移; flag == 0 判断没有基线漂移
%%
% 对数据进行归一化
data_dealed_max = mapminmax(dataIn, 0, 1); % 
data_dealed_min = mapminmax(-1*dataIn, 0, 1); %
num_division = num_epoch; % 分割的份数
% 关于 max 的部分
slope_max = acquire_slope(data_dealed_max, num_division);
% 关于 min 的部分
slope_min = acquire_slope(data_dealed_min, num_division);
slope_value = min(slope_max, slope_min);
%
part1 = slope_max > threshold_max;
part2 = slope_min > threshold_min;
%
shift_flag = part1 && part2;
end

