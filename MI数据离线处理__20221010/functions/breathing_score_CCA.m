function [score] = breathing_score_CCA(dataIn)
%% ���� CCA ���������۲��εĺ���ָ��
[n_channels, n_samples] = size(dataIn);
%% ���� 0.1Hz - 0.6Hz ������ģ��
t = 0:1:n_samples-1;
Fs = 50; % ����Ƶ��
t = t / Fs;
ideal_template = [];
fre_scope = 0.1:0.1:0.6; % ����Ƶ�ʷ�Χ
for i=1:1:length(fre_scope)
    part1 = sin(2*pi*fre_scope(1,i)*t);
    part2 = cos(2*pi*fre_scope(1,i)*t);
    ideal_template = [ideal_template;part1;part2];
end
%
% weights_xx��ģ���Ȩ��
[weights,weights_xx,coef,~,~,~] = canoncorr(transpose(dataIn), transpose(ideal_template));
score = coef;
end

