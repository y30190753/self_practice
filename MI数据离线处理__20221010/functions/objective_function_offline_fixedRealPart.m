function [ f ] = objective_function_offline_fixedRealPart( x )
%%
global param_objectiveFun; % Nc * Np
num_channel=size(param_objectiveFun,1);
%%
% param_vector = x(1:num_channel, 1) .* exp(1j .* x((num_channel+1):(num_channel*2), 1));
% 固定实数部分的权重
param_vector = ones(num_channel, 1) .* exp(1j .* x);
data_combined = param_vector.' * param_objectiveFun;

data_combined_phaseCurve = -(angle(data_combined));
data_combined_ideal = exp(-1j * data_combined_phaseCurve);
coef = corr(data_combined.', data_combined_ideal.');
%
% part1 = corr(transpose(real(data_combined)), transpose(real(data_combined_ideal))); % 实数部分相关系数
% part2 = corr(transpose(imag(data_combined)), transpose(imag(data_combined_ideal))); % 虚数部分相关系数
part3 = abs(coef); % 模长和实部都有效果
%%
% f= -(part1 + part2);
f = -part3;
end

