function algorithmForWaveopt_part(GUIdata,selected_index,subName,stage_number)
%%
%      1�����ò���ź�ѡ��rangeBin���Զ��彻���Ż�
%      2��ʹ�õ��Ǹ���Ȩ��
%%
global count_sub;
%%
GUIdata.rangeProf = GUIdata.rangeProf(:,selected_index(1):selected_index(end));
GUIdata.numFrame = 1:1:size(selected_index,2);
GUIdata.outHeart = GUIdata.outHeart(1,selected_index(1):selected_index(end));
GUIdata.outIndex = GUIdata.outIndex(1,selected_index(1):selected_index(end));
%% ���㷨1\2
old_optAlgorithm_flag = 0;
if(old_optAlgorithm_flag == 1)
    [outHeart, outPhase] = old_optAlgorithm(GUIdata, subName, selected_index);
elseif(old_optAlgorithm_flag == 0)
    outHeart=rand(1, length(GUIdata.numFrame));
    outPhase=rand(1, length(GUIdata.numFrame));
end
disp('���㷨������ϡ�');
%%
rangeProf = reshape(GUIdata.rangeProf,22,4,size(GUIdata.numFrame,2));
%% Ԥ��ѡ��bin
% ���ò���źŵĺ���Ƶ��ָ��
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
scores_temp = zeros(1, size(rawData_adjusted, 1)); % 1,20��confidenceMetric2
for i=1:1:length(scores_temp)
    fft_size=1024;
    temp_data = rawData_adjusted(i,1:fft_size);
    % scores_temp(1, i) = computeConfidenceMetric2(temp_data, 50, size(temp_data, 2)); % ����Ƶ������ָ��
    % scores_temp(1, i) = breathingdata_score_slope(-unwrap(angle(temp_data)));
    temp_curve=-unwrap(angle(temp_data));
    temp_curve_diff=diff(temp_curve);
    Fs=50;
    scores_temp(1, i)=computeConfidenceMetric2(temp_curve_diff,Fs,fft_size);
end
% ͨ�����
[~, scores_sorted_index] = sort(-1 * scores_temp); %
sorted_index_eight=scores_sorted_index(1,1:8); % ����ǰ8��ͨ��
scores_matrix_integer=zeros(5,4); % ������
scores_matrix_BS=zeros(5,4); % ������
channel_selected_accordingScores=[];
count_temp=0;
for itemp=1:1:size(scores_matrix_integer,1)
    for iitemp=1:1:size(scores_matrix_integer,2)
        count_temp=count_temp+1;
        if(ismember(count_temp,sorted_index_eight)) % 0/1����
            scores_matrix_integer(itemp,iitemp)=1;
            channel_selected_accordingScores=[channel_selected_accordingScores,count_temp];
        end
        % ����ռ��ԭʼ��
        scores_matrix_BS(itemp,iitemp)=scores_temp(1,count_temp);
    end
end
% ����5��rangeBin��scores
scores_5_rangeBins_integer=sum(scores_matrix_integer,2);
scores_5_rangeBins_BS=sum(scores_matrix_BS,2);
% [~,rangeBin_selected_index]=sort(-1*scores_5_rangeBins_BS);
[~,rangeBin_selected_index]=sort(-1*scores_5_rangeBins_integer); % ����0/1��������
num_rangeBin_selected = 5; % rangeBin��ѡ������
rangeBin_selected_index=rangeBin_selected_index(1:num_rangeBin_selected);
rangeBin_selected_index=sort(rangeBin_selected_index);
rangeBin_scope = 10:1:14; % rangeBin��ѡ��Χ
rangeBin_scope_selected = rangeBin_scope(rangeBin_selected_index); % 12, 13
channel_selected_index=[];
for itemp=1:1:size(rangeBin_scope_selected,2)
    temp_index=((rangeBin_scope_selected(1,itemp)-9-1)*num_Rx+1):1:((rangeBin_scope_selected(1,itemp)-9-1)*num_Rx+4);
    channel_selected_index=[channel_selected_index,temp_index];
end
% channel_selected_index=channel_selected_accordingScores; % ����scoresѡ���ͨ��
channel_selected_index=channel_selected_index; % ����0/1����ѡ���ͨ��
% ͨ��ѡ�� - ͼ��Ч��
fig_temp=figure;
count_temp=0;
for i=1:1:length(scores_temp)
    count_temp=count_temp+1;
    subplot(5,4,count_temp);
    % if(ismember(i, channel_selected_accordingScores)) % scores����ͨ��
    if(ismember(i, channel_selected_index)) % 0/1��������ͨ��
        plot(-unwrap(angle(rawData_adjusted(i,:))),'r');
    else
        plot(-unwrap(angle(rawData_adjusted(i,:))),'b');
    end
    title(num2str(scores_temp(1, i)));
end
figure_fullscreen(fig_temp);
%% �������
numFrames = size(rawData,3);
dealed_data_val_record = zeros(1,numFrames);
model_update_record = zeros(1,numFrames);
rangeBin_switch_record = zeros(1,numFrames);
opt_dectection_record = zeros(1,numFrames);
delta_dataLinear_record = zeros(1,numFrames);
len_rawData_buffer = 512; % origin len
len_opt = 512; % �㷨ִ�еĻ���������
rawData_buffer = zeros(size(rawData_adjusted,1),len_rawData_buffer);
len_dealed_data_buffer = len_rawData_buffer;
dealed_data_buffer = zeros(1,len_dealed_data_buffer);
len_dealedData_phase_buffer = len_rawData_buffer;
dealed_data_phase_buffer = zeros(1,len_dealedData_phase_buffer); 
slope = 0; % �жϴ�Ư�Ƶ�б��
pre_phase = 0;
cur_phase = 0;
delta_phase = cur_phase - pre_phase;
dealed_data_phase_record = zeros(1,numFrames);
dealed_data_val_record_compared = zeros(1, numFrames);
pre_phase_compare = 0;
cur_phase_compare = 0;
delta_phase_compare = cur_phase_compare - pre_phase_compare;
phase_record_compare = zeros(1,numFrames);
num_channels = size(channel_selected_index,2); % ������ϵ�ͨ������
cur_weights = zeros(num_channels, 1);
first_data_flag = 0;
threshold_slope = 0;
count = 0; % ���Ǽ�����
modelUpdate_flag = 0;
num_opt = 10; % �Զ����Ż����������
for frame_index=1:1:numFrames
    count = count + 1;
    rawData_current  = rawData_adjusted(:,frame_index);
    rawData_buffer = circshift(rawData_buffer,[0,-1]);
    rawData_buffer(:,end) = rawData_current;
    if(count < len_opt || count < len_rawData_buffer) % δ����һ�����ݿ�
        continue;
    elseif(count == len_rawData_buffer) % ��һ�����ݿ�
        % ѡ�� rangeBin
        rawData_source = rawData_buffer(:,end+1-len_opt:end); % ��ص����ݳ��Ⱥ�ԭʼ���ݻ�����֮����ڲ�
        rawData_source_selected=rawData_buffer(channel_selected_index,end+1-len_opt:end);
        [ cur_weights, ~ ] = opt_algorithm( rawData_source_selected,cur_weights,1,num_opt ); % Ȩ�ص���ȡ
        % �ģ��Զ��彻���Ż�
        % �� theta Ϊ������λ���� ʵ������ Ϊ�����ͨ����Ӱ��
        % cur_weights = weights_opt_alternative(rawData_source_selected,cur_weights,1,num_opt); % �Զ��彻���Ż�
        %
        % reshapeȨ��
        cur_weights = cur_weights(1:num_channels, 1) .* exp(1j .* cur_weights((num_channels+1):(num_channels*2), 1));
        %
        dealed_data_temp = cur_weights.' * rawData_source_selected;
        slope = acquire_slope(-unwrap(angle(dealed_data_temp(1,:))));
        dealed_data_val_record(1,1:frame_index) = dealed_data_temp(1,:);
        dealed_data_val_record_compared(1,1:frame_index) = dealed_data_temp(1,:);
        dealed_data_phase_record(1,1:frame_index) = -unwrap(angle(dealed_data_temp));
        phase_record_compare(1,1:frame_index) = -unwrap(angle(dealed_data_temp));
        % ����dealed_data_phase_buffer
        dealed_data_phase_buffer(1,:) = dealed_data_phase_record(1,end+1-len_dealedData_phase_buffer:end);
        dealed_data_buffer(1,:) = dealed_data_val_record(1,max(frame_index-len_dealed_data_buffer+1,1):frame_index);
        cur_phase = dealed_data_phase_record(1,frame_index);
        cur_phase_compare = dealed_data_phase_record(1,frame_index);
        pre_weights = cur_weights;
        opt_dectection_record(1, count) = 1; %
        first_data_flag = 1;
        disp('��һ�����ݿ鴦����ϣ�');
        %
    else % ��2��֮�������
        rawData_current_selected = rawData_current(channel_selected_index, 1); % Nc*1
        dealed_data_current_selected = cur_weights.' * rawData_current_selected;
        cur_phase = angle_unwrap(-angle(dealed_data_current_selected),pre_phase); % �����
        % ������λ����
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
        % buffer �϶����ڽ�ǰ���
        dealed_data_buffer = circshift(dealed_data_buffer,[0,-1]);
        dealed_data_buffer(1,end) = dealed_data_current_selected;
        % �Ƿ��������ﱣ��
        dealed_data_val_record(1,frame_index) = dealed_data_current_selected;
        dealed_data_phase_record(1,frame_index) = cur_phase + delta_phase;
        phase_record_compare(1, frame_index) = cur_phase_compare + delta_phase_compare;
        dealed_data_val_record_compared(1, frame_index) = dealed_data_current_selected;
        % ָ��������ڣ�����Ĳ�����rawData_buffer
        if(mod(count,len_opt) == 0)
            opt_dectection_record(1,frame_index) = 1;
            % �ж��Ƿ�����˴�Ư����ǿ�Ƹ���
            data_temp = -unwrap(angle(dealed_data_buffer(1,:)));
            slope = acquire_slope(data_temp); % ��õ�ǰ���ݿ���������б��
            threshold_slope = 0.6;
            if(slope > threshold_slope)
                modelUpdate_flag = 1;
            end
        end
        %
        if(modelUpdate_flag == 1) % ����ģ�Ͳ���
            [cur_weights_new, ~] = opt_algorithm(rawData_buffer(channel_selected_index,:),cur_weights,1,num_opt); % Ȩ�ص���ȡ
            % cur_weights_new = weights_opt_alternative(rawData_buffer(channel_selected_index,:),cur_weights,1,num_opt); % �Զ��彻���Ż�
            pre_weights = cur_weights;
            cur_weights = cur_weights_new(1:num_channels, 1) .* exp(1j .* cur_weights_new((num_channels+1):(num_channels*2), 1));
            model_update_record(1,min(frame_index+20,numFrames)) = 1;
            % ��֮ǰ������Ҳ�õ���Ȩ�������Ż�һ�£�����Ҫ����ͨ��ѡ�������
            dealed_data_val_record_compared(1, 1: frame_index) = cur_weights.' * rawData_adjusted(channel_selected_index,1:frame_index);
        end
    end
    pre_phase = dealed_data_phase_record(1,frame_index);
    pre_phase_compare = phase_record_compare(1,frame_index);
    delta_dataLinear_record(1,frame_index) = slope;
    %
    disp([subName,',��ͨ���㷨����֡��: ',num2str(frame_index)]);
end
%% ��ͼ
fig = figure;
%[ha, pos] = tight_subplot(��, ��, [���¼��,���Ҽ��],[�±߾�,�ϱ߾� ], [��߾�,�ұ߾� ])
[ha, ~] = tight_subplot(6, 1, [.05 .05], [.03 .03], [.19 .19]);
% subplot(611);
axes(ha(1));
plot(outHeart/8000);
title('ԭ�㷨��ȡ�����ź�');axis([0,numFrames,min(outHeart/8000),max(outHeart/8000)]);
% subplot(612);
axes(ha(2));
plot(filter_RC(outPhase));
title('���㷨��ȡ�����ź�');axis([0,numFrames,min(filter_RC(outPhase)),max(filter_RC(outPhase))]);
%
% subplot(613);
axes(ha(3));
plot(delta_dataLinear_record);title('�������Եı仯');
hold on; t1 = (zeros(1,numFrames) + 1) * threshold_slope; plot(t1,'r');
axis([0,numFrames,min(delta_dataLinear_record),max(delta_dataLinear_record)]);
%
% subplot(713);
% plot(dealed_data_phase_record);
% hold on;
% for tempi=1:1:(numFrames)
%     if(model_update_record(1,tempi) == 1)
%         plot([tempi,tempi],[-10000,10000],'y'); % ��ɫ�Ǹ���ģ�Ͳ���
%         hold on;
%     end
%     if(opt_dectection_record(1,tempi) == 1)
%         plot([tempi,tempi],[-10000,10000],'b'); % ��ɫ�����ݿ�ĳ���
%         hold on;
%     end
%     if(rangeBin_switch_record(1,tempi) == 1)
%         plot([tempi,tempi],[-10000,10000],'r'); % ��ɫ���л�rangeBin
%         hold on;
%     end
% end
% axis([0,numFrames,min(dealed_data_phase_record),max(dealed_data_phase_record)]);
% title(['��ͨ������㷨']);
%% �ԱȵĽ����ʽ
% subplot(614);
axes(ha(4));
plot(phase_record_compare);title(['unwrap��buffer�Ľ����ʽ, ', subName, ', ',num2str(stage_number)]);
hold on;
for tempi=1:1:(numFrames)
    if(model_update_record(1,tempi) == 1)
        plot([tempi,tempi],[-10000,10000],'y', 'LineWidth', 0.01); % ��ɫ�Ǹ���ģ�Ͳ���
        hold on;
    end
    if(opt_dectection_record(1,tempi) == 1)
        plot([tempi,tempi],[-10000,10000],'b', 'LineWidth', 0.01); % ��ɫ�����ݿ�ĳ���
        hold on;
    end
    if(rangeBin_switch_record(1,tempi) == 1)
        plot([tempi,tempi],[-10000,10000],'r', 'LineWidth', 0.01); % ��ɫ���л�rangeBin
        hold on;
    end
end
axis([0,numFrames,min(phase_record_compare),max(phase_record_compare)]);
%
% subplot(616);
axes(ha(6));
% ���źŽ���ż����
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
% ��ֵ�˲�
medfilt_para = 55; % debug = 40
temp_medfilted = medfilt1(temp_evenContinuation, medfilt_para);
% �°�����
envelope_points = 512; % debug = 300
[~, down] = envelope(temp_medfilted, envelope_points,'peak');
down = down(1, continuation_points+1:continuation_points+size(temp, 2));
plot(temp - down);
axis([0,numFrames,min(temp - down),max(temp - down)]);
title('��ʹ�����һ��Ȩ�� - ȥ���°�����');
%
% subplot(615);
axes(ha(5));
temp = -unwrap(angle(dealed_data_val_record_compared));
plot(temp);
hold on;
plot(down);
axis([0,numFrames,min(temp),max(temp)]);
title('��ʹ�����һ��Ȩ��');
figure_fullscreen(fig);
%
saveas(fig_temp,[subName,num2str(stage_number),'_ͨ��ѡ��']);
saveas(fig,[subName,num2str(stage_number)]);
disp([subName, ' is done! ���: ', num2str(count_sub),', �׶κ�: ', num2str(stage_number)]);
end

