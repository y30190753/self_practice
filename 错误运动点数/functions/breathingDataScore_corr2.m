function [ score ] = breathingDataScore_corr2( dataIn )
%%                    根据不同方式计算出来的呼吸曲线相关性计算评价指标
%     dataIn： 1 * 512
%%
breathing_curve1 = -unwrap(angle(dataIn));
breathing_curve2 = acquire_breathCurve_by_another_way(dataIn);
score = corr(transpose(breathing_curve1), transpose(breathing_curve2));
% % % 画图
% figure;
% subplot(2,1,1);
% plot(breathing_curve1);
% subplot(2,1,2);
% plot(breathing_curve2);
% title(num2str(score));
end

