function [score] = rollingCorr(dataIn)
%%
%           输入： dataIn：1*512
%
%%
% pi/2 划分成10份 pi/20
angle_set = pi/20:pi/20:pi/2;
corr_set = zeros(1, size(angle_set,2)); % 10个相关系数

for i=1:1:length(angle_set)
    %
    phase_change = exp(1j * angle_set(1,i));
    temp_data = phase_change * dataIn; %
    % 模板
    %     temp_data_template = exp(1j * angle(temp_data)); %
    %     corr_part1 = corr(transpose(real(temp_data)), transpose(real(temp_data_template)));
    %     corr_part2 = corr(transpose(imag(temp_data)), transpose(imag(temp_data_template)));
    %     corr_set(1, i) = 0.5 * corr_part1 + 0.5 * corr_part2;
    %     corr_set(1, i) = corr(transpose(compared), transpose(ones(1,size(temp_data,2))));
    corr_set(1, i) = temp_score;
end
%
score = sum(corr_set) / size(corr_set, 2);
end

