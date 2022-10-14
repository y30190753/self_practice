function algorithmForWaveopt_part(GUIdata,selected_index,subName,stage_number)
%%
%      1、利用差分信号选择rangeBin，自定义交替优化
%      2、使用的是复数权重
%%
global count_sub;
%%
GUIdata.rangeProf = GUIdata.rangeProf(:,selected_index(1):selected_index(end));
GUIdata.numFrame = 1:1:size(selected_index,2);
GUIdata.outHeart = GUIdata.outHeart(1,selected_index(1):selected_index(end));
GUIdata.outIndex = GUIdata.outIndex(1,selected_index(1):selected_index(end));
%% 旧算法1\2
old_optAlgorithm_flag = 0;
if(old_optAlgorithm_flag == 1)
    [outHeart, outPhase] = old_optAlgorithm(GUIdata, subName, selected_index);
elseif(old_optAlgorithm_flag == 0)
    outHeart=rand(1, length(GUIdata.numFrame));
    outPhase=rand(1, length(GUIdata.numFrame));
end
disp('老算法运行完毕。');
%%
rangeProf = reshape(GUIdata.rangeProf,22,4,size(GUIdata.numFrame,2));
%% 预先选择bin
% 利用差分信号的呼吸频谱指标
% rawData = rangeProf(rangeBin_scope(1,:),:,:); % 5*4*8000
rawData = rangeProf(10:1:14, :, :); % 5*4*selected
rawData_adjusted = zeros(size(rawData,1) * 4,size(rawData,3)); % 20 * 8000
num_Rx=4;
count_temp = 0;
for i=1:1:size(rawData, 1) % 5
    for ii=1:1:num_Rx
        count_temp = count_temp + 1;
        rawData_adjusted(count_temp,:) = squeeze(rawData(i,ii,:));
    end
end
scores_temp = zeros(1, size(rawData_adjusted, 1)); % 1,20，confidenceMetric2
for i=1:1:length(scores_temp)
    fft_size=1024;
    temp_data = rawData_adjusted(i,1:fft_size);
    % scores_temp(1, i) = computeConfidenceMetric2(temp_data, 50, size(temp_data, 2)); % 呼吸频谱评价指标
    % scores_temp(1, i) = breathingdata_score_slope(-unwrap(angle(temp_data)));
    temp_curve=-unwrap(angle(temp_data));
    temp_curve_diff=diff(temp_curve);
    Fs=50;
    scores_temp(1, i)=computeConfidenceMetric2(temp_curve_diff,Fs,fft_size);
end
% 通道算分
[~, scores_sorted_index] = sort(-1 * scores_temp); %
sorted_index_eight=scores_sorted_index(1,1:8); % 排名前8的通道
scores_matrix_integer=zeros(5,4); % 积分牌
scores_matrix_BS=zeros(5,4); % 积分牌
channel_selected_accordingScores=[];
count_temp=0;
for itemp=1:1:size(scores_matrix_integer,1)
    for iitemp=1:1:size(scores_matrix_integer,2)
        count_temp=count_temp+1;
        if(ismember(count_temp,sorted_index_eight)) % 0/1积分
            scores_matrix_integer(itemp,iitemp)=1;
            channel_selected_accordingScores=[channel_selected_accordingScores,count_temp];
        end
        % 能量占比原始分
        scores_matrix_BS(itemp,iitemp)=scores_temp(1,count_temp);
    end
end
% 计算5个rangeBin的scores
scores_5_rangeBins_integer=sum(scores_matrix_integer,2);
scores_5_rangeBins_BS=sum(scores_matrix_BS,2);
% [~,rangeBin_selected_index]=sort(-1*scores_5_rangeBins_BS);
[~,rangeBin_selected_index]=sort(-1*scores_5_rangeBins_integer); % 根据0/1积分排序
num_rangeBin_selected = 5; % rangeBin的选择数量
rangeBin_selected_index=rangeBin_selected_index(1:num_rangeBin_selected);
rangeBin_selected_index=sort(rangeBin_selected_index);
rangeBin_scope = 10:1:14; % rangeBin的选择范围
rangeBin_scope_selected = rangeBin_scope(rangeBin_selected_index); % 12, 13
channel_selected_index=[];
for itemp=1:1:size(rangeBin_scope_selected,2)
    temp_index=((rangeBin_scope_selected(1,itemp)-9-1)*num_Rx+1):1:((rangeBin_scope_selected(1,itemp)-9-1)*num_Rx+4);
    channel_selected_index=[channel_selected_index,temp_index];
end
% channel_selected_index=channel_selected_accordingScores; % 根据scores选择的通道
channel_selected_index=channel_selected_index; % 根据0/1积分选择的通道
% 通道选择 - 图形效果
fig_temp=figure;
count_temp=0;
for i=1:1:length(scores_temp)
    count_temp=count_temp+1;
    subplot(5,4,count_temp);
    % if(ismember(i, channel_selected_accordingScores)) % scores排序通道
    if(ismember(i, channel_selected_index)) % 0/1积分排序通道
        plot(-unwrap(angle(rawData_adjusted(i,:))),'r');
    else
        plot(-unwrap(angle(rawData_adjusted(i,:))),'b');
    end
    title(num2str(scores_temp(1, i)));
end
figure_fullscreen(fig_temp);
%% 主体变量
numFrames = size(rawData,3);
dealed_data_val_record = zeros(1,numFrames);
model_update_record = zeros(1,numFrames);
rangeBin_switch_record = zeros(1,numFrames);
opt_dectection_record = zeros(1,numFrames);
delta_dataLinear_record = zeros(1,numFrames);
len_rawData_buffer = 512; % origin len
len_opt = 512; % 算法执行的基础数据量
rawData_buffer = zeros(size(rawData_adjusted,1),len_rawData_buffer);
len_dealed_data_buffer = len_rawData_buffer;
dealed_data_buffer = zeros(1,len_dealed_data_buffer);
len_dealedData_phase_buffer = len_rawData_buffer;
dealed_data_phase_buffer = zeros(1,len_dealedData_phase_buffer); 
slope = 0; % 判断大漂移的斜率
pre_phase = 0;
cur_phase = 0;
delta_phase = cur_phase - pre_phase;
dealed_data_phase_record = zeros(1,numFrames);
dealed_data_val_record_compared = zeros(1, numFrames);
pre_phase_compare = 0;
cur_phase_compare = 0;
delta_phase_compare = cur_phase_compare - pre_phase_compare;
phase_record_compare = zeros(1,numFrames);
num_channels = size(channel_selected_index,2); % 用来结合的通道数量
cur_weights = zeros(num_channels, 1);
first_data_flag = 0;
threshold_slope = 0;
count = 0; % 又是计数器
modelUpdate_flag = 0;
num_opt = 10; % 自定义优化的随机点数
for frame_index=1:1:numFrames
    count = count + 1;
    rawData_current  = rawData_adjusted(:,frame_index);
    rawData_buffer = circshift(rawData_buffer,[0,-1]);
    rawData_buffer(:,end) = rawData_current;
    if(count < len_opt || count < len_rawData_buffer) % 未满第一个数据块
        continue;
    elseif(count == len_rawData_buffer) % 第一个数据块
        % 选择 rangeBin
        rawData_source = rawData_buffer(:,end+1-len_opt:end); % 监控的数据长度和原始数据缓存器之间存在差
        rawData_source_selected=rawData_buffer(channel_selected_index,end+1-len_opt:end);
        [ cur_weights, ~ ] = opt_algorithm( rawData_source_selected,cur_weights,1,num_opt ); % 权重的求取
        % 改：自定义交替优化
        % 先 theta 为对齐相位，后 实数部分 为扩大好通道的影响
        % cur_weights = weights_opt_alternative(rawData_source_selected,cur_weights,1,num_opt); % 自定义交替优化
        %
        % reshape权重
        cur_weights = cur_weights(1:num_channels, 1) .* exp(1j .* cur_weights((num_channels+1):(num_channels*2), 1));
        %
        dealed_data_temp = cur_weights.' * rawData_source_selected;
        slope = acquire_slope(-unwrap(angle(dealed_data_temp(1,:))));
        dealed_data_val_record(1,1:frame_index) = dealed_data_temp(1,:);
        dealed_data_val_record_compared(1,1:frame_index) = dealed_data_temp(1,:);
        dealed_data_phase_record(1,1:frame_index) = -unwrap(angle(dealed_data_temp));
        phase_record_compare(1,1:frame_index) = -unwrap(angle(dealed_data_temp));
        % 存入dealed_data_phase_buffer
        dealed_data_phase_buffer(1,:) = dealed_data_phase_record(1,end+1-len_dealedData_phase_buffer:end);
        dealed_data_buffer(1,:) = dealed_data_val_record(1,max(frame_index-len_dealed_data_buffer+1,1):frame_index);
        cur_phase = dealed_data_phase_record(1,frame_index);
        cur_phase_compare = dealed_data_phase_record(1,frame_index);
        pre_weights = cur_weights;
        opt_dectection_record(1, count) = 1; %
        first_data_flag = 1;
        disp('第一个数据块处理完毕！');
        %
    else % 第2块之后的数据
        rawData_current_selected = rawData_current(channel_selected_index, 1); % Nc*1
        dealed_data_current_selected = cur_weights.' * rawData_current_selected;
        cur_phase = angle_unwrap(-angle(dealed_data_current_selected),pre_phase); % 解缠绕
        % 存入相位缓存
        dealed_data_phase_buffer = circshift(dealed_data_phase_buffer,[0,-1]);
        dealed_data_phase_buffer(1,end) = -angle(dealed_data_current_selected);
        dealed_data_phase_buffer(1,:) = unwrap(dealed_data_phase_buffer(1,:));
        cur_phase_compare = dealed_data_phase_buffer(1,end);
        if(modelUpdate_flag == 1)
            delta_phase = pre_phase - cur_phase;
            delta_phase_compare = pre_phase_compare - cur_phase_compare;
            modelUpdate_flag = 0;
        end
        %
        if(first_data_flag == 1)
            delta_phase_compare = pre_phase_compare - cur_phase_compare;
            first_data_flag = 0;
        end
        %
        % buffer 肯定是在较前面的
        dealed_data_buffer = circshift(dealed_data_buffer,[0,-1]);
        dealed_data_buffer(1,end) = dealed_data_current_selected;
        % 是否是在这里保存
        dealed_data_val_record(1,frame_index) = dealed_data_current_selected;
        dealed_data_phase_record(1,frame_index) = cur_phase + delta_phase;
        phase_record_compare(1, frame_index) = cur_phase_compare + delta_phase_compare;
        dealed_data_val_record_compared(1, frame_index) = dealed_data_current_selected;
        % 指标更新周期，计算的材料是rawData_buffer
        if(mod(count,len_opt) == 0)
            opt_dectection_record(1,frame_index) = 1;
            % 判断是否出现了大漂移需强制更新
            data_temp = -unwrap(angle(dealed_data_buffer(1,:)));
            slope = acquire_slope(data_temp); % 求得当前数据块的线性拟合斜率
            threshold_slope = 0.6;
            if(slope > threshold_slope)
                modelUpdate_flag = 1;
            end
        end
        %
        if(modelUpdate_flag == 1) % 更新模型参数
            [cur_weights_new, ~] = opt_algorithm(rawData_buffer(channel_selected_index,:),cur_weights,1,num_opt); % 权重的求取
            % cur_weights_new = weights_opt_alternative(rawData_buffer(channel_selected_index,:),cur_weights,1,num_opt); % 自定义交替优化
            pre_weights = cur_weights;
            cur_weights = cur_weights_new(1:num_channels, 1) .* exp(1j .* cur_weights_new((num_channels+1):(num_channels*2), 1));
            model_update_record(1,min(frame_index+20,numFrames)) = 1;
            % 把之前的数据也用档次权重重新优化一下，这里要考虑通道选择的问题
            dealed_data_val_record_compared(1, 1: frame_index) = cur_weights.' * rawData_adjusted(channel_selected_index,1:frame_index);
        end
    end
    pre_phase = dealed_data_phase_record(1,frame_index);
    pre_phase_compare = phase_record_compare(1,frame_index);
    delta_dataLinear_record(1,frame_index) = slope;
    %
    disp([subName,',多通道算法处理帧数: ',num2str(frame_index)]);
end
%% 画图
fig = figure;
%[ha, pos] = tight_subplot(行, 列, [上下间距,左右间距],[下边距,上边距 ], [左边距,右边距 ])
[ha, ~] = tight_subplot(6, 1, [.05 .05], [.03 .03], [.19 .19]);
% subplot(611);
axes(ha(1));
plot(outHeart/8000);
title('原算法提取呼吸信号');axis([0,numFrames,min(outHeart/8000),max(outHeart/8000)]);
% subplot(612);
axes(ha(2));
plot(filter_RC(outPhase));
title('新算法提取呼吸信号');axis([0,numFrames,min(filter_RC(outPhase)),max(filter_RC(outPhase))]);
%
% subplot(613);
axes(ha(3));
plot(delta_dataLinear_record);title('数据线性的变化');
hold on; t1 = (zeros(1,numFrames) + 1) * threshold_slope; plot(t1,'r');
axis([0,numFrames,min(delta_dataLinear_record),max(delta_dataLinear_record)]);
%
% subplot(713);
% plot(dealed_data_phase_record);
% hold on;
% for tempi=1:1:(numFrames)
%     if(model_update_record(1,tempi) == 1)
%         plot([tempi,tempi],[-10000,10000],'y'); % 黄色是更新模型参数
%         hold on;
%     end
%     if(opt_dectection_record(1,tempi) == 1)
%         plot([tempi,tempi],[-10000,10000],'b'); % 蓝色是数据块的长度
%         hold on;
%     end
%     if(rangeBin_switch_record(1,tempi) == 1)
%         plot([tempi,tempi],[-10000,10000],'r'); % 红色是切换rangeBin
%         hold on;
%     end
% end
% axis([0,numFrames,min(dealed_data_phase_record),max(dealed_data_phase_record)]);
% title(['多通道结合算法']);
%% 对比的解缠方式
% subplot(614);
axes(ha(4));
plot(phase_record_compare);title(['unwrap和buffer的解缠方式, ', subName, ', ',num2str(stage_number)]);
hold on;
for tempi=1:1:(numFrames)
    if(model_update_record(1,tempi) == 1)
        plot([tempi,tempi],[-10000,10000],'y', 'LineWidth', 0.01); % 黄色是更新模型参数
        hold on;
    end
    if(opt_dectection_record(1,tempi) == 1)
        plot([tempi,tempi],[-10000,10000],'b', 'LineWidth', 0.01); % 蓝色是数据块的长度
        hold on;
    end
    if(rangeBin_switch_record(1,tempi) == 1)
        plot([tempi,tempi],[-10000,10000],'r', 'LineWidth', 0.01); % 红色是切换rangeBin
        hold on;
    end
end
axis([0,numFrames,min(phase_record_compare),max(phase_record_compare)]);
%
% subplot(616);
axes(ha(6));
% 对信号进行偶延拓
continuation_points = 1000;
temp = -unwrap(angle(dealed_data_val_record_compared));
data_temp = temp(1,1:continuation_points); data_evenContinuation_front = zeros(1,continuation_points);
for i = 1:1:continuation_points
    data_evenContinuation_front(1, i) = data_temp(1, continuation_points-(i-1));
end
data_temp = temp(1, end-continuation_points+1:end); data_evenContinuation_end = zeros(1,continuation_points);
for i = 1:1:continuation_points
    data_evenContinuation_end(1, i) = data_temp(1, continuation_points-(i-1));
end
temp_evenContinuation = [data_evenContinuation_front, temp, data_evenContinuation_end];
% 中值滤波
medfilt_para = 55; % debug = 40
temp_medfilted = medfilt1(temp_evenContinuation, medfilt_para);
% 下包络线
envelope_points = 512; % debug = 300
[~, down] = envelope(temp_medfilted, envelope_points,'peak');
down = down(1, continuation_points+1:continuation_points+size(temp, 2));
plot(temp - down);
axis([0,numFrames,min(temp - down),max(temp - down)]);
title('仅使用最后一组权重 - 去除下包络线');
%
% subplot(615);
axes(ha(5));
temp = -unwrap(angle(dealed_data_val_record_compared));
plot(temp);
hold on;
plot(down);
axis([0,numFrames,min(temp),max(temp)]);
title('仅使用最后一组权重');
figure_fullscreen(fig);
%
saveas(fig_temp,[subName,num2str(stage_number),'_通道选择']);
saveas(fig,[subName,num2str(stage_number)]);
disp([subName, ' is done! 序号: ', num2str(count_sub),', 阶段号: ', num2str(stage_number)]);
end

