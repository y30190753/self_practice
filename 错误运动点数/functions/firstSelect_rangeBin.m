function [ rangeBin_indices, weights_res, scores_5 ] = firstSelect_rangeBin( data_epoch )
%%
%    输入：data_epoch   20 * len

%% 选择同一个rangebin的天线数据
rawData_reshaped = zeros(5,4,size(data_epoch,2));
for i=1:1:5
    for j=1:1:4
        rawData_reshaped(i,j,:) = squeeze(data_epoch((i-1)*4+j,:));
    end
end
breathing_scores = zeros(1,5);
weights_scores = zeros(8,5); % 记录选择rangeBin时的权重
weights_initial = zeros(8,1);
weights_initial(1,1) = 1;
for i=1:1:5
    [weights_temp, ~] = opt_algorithm(squeeze(rawData_reshaped(i,:,:)),weights_initial,1,10);
    dealed_data = [weights_temp(1,1)*exp(1j*weights_temp(5,1)),weights_temp(2,1)*exp(1j*weights_temp(6,1)),...
        weights_temp(3,1)*exp(1j*weights_temp(7,1)),weights_temp(4,1)*exp(1j*weights_temp(8,1))] * squeeze(rawData_reshaped(i,:,:));
    
    breathing_scores(1,i) = computeConfidenceMetric2(diff(-unwrap(angle(dealed_data))),50,size(dealed_data,2));
    
    weights_scores(:,i) = weights_temp(:,1);
end

[~, rangeBin_indices] = max(breathing_scores);
scores_5 = breathing_scores;
% rangeBin_indices = 5; % debug
weights_res = weights_scores(:,:);
%
%% 跨 rangeBin 选择通道
% breathing_scores = zeros(1,20);
% for i=1:1:length(breathing_scores)
%     breathing_scores(1,i) = computeConfidenceMetric2(diff(-unwrap(angle(data_epoch(i,:)))),50,size(data_epoch,2));
% end
% num_channels = 4;
% [~,indices] = sort(-1 * breathing_scores(1,:));
% rangeBin_indices = indices(1,1:num_channels);


end

