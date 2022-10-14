function [shift_flag, slope_value] = shift_determine(dataIn, num_epoch, threshold_max, threshold_min)
%%  �ж�һ�����ݿ��Ƿ���ֻ���Ư��
%      ���� dataIn, 1 * n_points
%      ��� flag == 1 �ж��л���Ư��; flag == 0 �ж�û�л���Ư��
%%
% �����ݽ��й�һ��
data_dealed_max = mapminmax(dataIn, 0, 1); % 
data_dealed_min = mapminmax(-1*dataIn, 0, 1); %
num_division = num_epoch; % �ָ�ķ���
% ���� max �Ĳ���
slope_max = acquire_slope(data_dealed_max, num_division);
% ���� min �Ĳ���
slope_min = acquire_slope(data_dealed_min, num_division);
slope_value = min(slope_max, slope_min);
%
part1 = slope_max > threshold_max;
part2 = slope_min > threshold_min;
%
shift_flag = part1 && part2;
end

