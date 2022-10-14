function [score] = breathing_score_CCA(dataIn)
%% 利用 CCA 相似性评价波形的呼吸指数
[n_channels, n_samples] = size(dataIn);
%% 产生 0.1Hz - 0.6Hz 的理想模板
t = 0:1:n_samples-1;
Fs = 50; % 采样频率
t = t / Fs;
ideal_template = [];
fre_scope = 0.1:0.1:0.6; % 呼吸频率范围
for i=1:1:length(fre_scope)
    part1 = sin(2*pi*fre_scope(1,i)*t);
    part2 = cos(2*pi*fre_scope(1,i)*t);
    ideal_template = [ideal_template;part1;part2];
end
%
% weights_xx是模板的权重
[weights,weights_xx,coef,~,~,~] = canoncorr(transpose(dataIn), transpose(ideal_template));
score = coef;
end

