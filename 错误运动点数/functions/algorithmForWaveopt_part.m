function algorithmForWaveopt_part(GUIdata,selected_index,subName,stage_number)
%%
%
%
%
%%
global count_sub;
%%
GUIdata.rangeProf = GUIdata.rangeProf(:,selected_index(1):selected_index(end));
GUIdata.numFrame = 1:1:size(selected_index,2);
GUIdata.outHeart = GUIdata.outHeart(1,selected_index(1):selected_index(end));
GUIdata.outIndex = GUIdata.outIndex(1,selected_index(1):selected_index(end));
%
%% 旧算法1\2
% [outHeart, outPhase] = old_optAlgorithm(GUIdata, subName, selected_index);
disp('老算法运行完毕。');
%%
rangeProf = reshape(GUIdata.rangeProf,22,4,size(GUIdata.numFrame,2));
rangeBin_scope = 10:1:14; % rangeBin的选择范围
num_rangeBins = size(rangeBin_scope, 2);
% rawData = rangeProf(rangeBin_scope(1,:),:,:); % 5*4*8000
rawData = rangeProf(rangeBin_scope(1,:), :, :); % 5*4*selected
rawData_adjusted = zeros(num_rangeBins * 4, size(rawData,3)); % 20 * 8000
count_temp = 0;
for i=1:1:num_rangeBins
    for ii=1:1:4
        count_temp = count_temp + 1;
        rawData_adjusted(count_temp,:) = squeeze(rawData(i,ii,:));
    end
end
numFrames = size(rawData,3);
% 缓存器
len_rawData_buffer = 512; % 数据片段长度
rawData_buffer = zeros(num_rangeBins * 4, len_rawData_buffer);
count_numPictures = 0;
count_epoch=0;
%
count = 0;
for frame_index = 1:1:numFrames
    count = count + 1;
    rawData_buffer = circshift(rawData_buffer, [0,-1]);
    rawData_buffer(:,end) = rawData_adjusted(:, frame_index);
    %
    if (mod(count, len_rawData_buffer) == 0)
        % 算呼吸频谱
        breathing_spectrum_scores = zeros(1,size(rawData_buffer,1));
        breathing_spectrum_scores_normalized = zeros(1,size(rawData_buffer,1)); % 复数域数据进行归一化的
        absSpectrum_record=zeros(size(rawData_buffer,1),len_rawData_buffer);
        sum_3_part_record = zeros(size(rawData_buffer,1),3);
        for itemp=1:1:size(rawData_buffer,1)
            temp_data = -unwrap(angle(rawData_buffer(itemp,:)));
            temp_data2 = rawData_buffer(itemp,:); % 原始复数域数据
            temp_data = diff(temp_data); % 差分
            %             breathing_spectrum_scores(1,itemp)=computeConfidenceMetric2(temp_data,50,len_rawData_buffer);
            % 复数域数据的错误圆心点数
            breathing_spectrum_scores(1,itemp)=breathingScore_circCenterPosition(temp_data2);
            % 归一化
            temp_data2_nor = temp_data2 ./ abs(temp_data2); % way1
            %             temp_data2_nor = exp(-1j * -unwrap(angle(rawData_buffer(itemp,:)))); % way2
            breathing_spectrum_scores_normalized(1, itemp) = breathingScore_circCenterPosition(temp_data2_nor);
        end
        max_breathing_scores = max(breathing_spectrum_scores);
        max_breathing_scores_nor = max(breathing_spectrum_scores_normalized);
        %
        %% 原始波形，未差分
        fig1 = figure;count_temp = 0;
        % 画布划分
        %[ha, pos] = tight_subplot(行, 列, [上下间距,左右间距],[下边距,上边距 ], [左边距,右边距 ])
        [ha, pos] = tight_subplot(5, 4*2, [.05 .02], [.06 .06], [.06 .06]);
        for itemp = 1:1:20
            count_temp = count_temp + 1;
            count_app=count_temp + 4*(floor((count_temp-1)/4));
            temp_data1 = rawData_buffer(itemp,:);
            temp_data2 = -unwrap(angle(temp_data1));
            %             temp_data2 = diff(temp_data2); % 差分
            %             subplot(4,5,count_temp);
            axes(ha(count_app));
            if(breathing_spectrum_scores(1,itemp)==max_breathing_scores)
                plot(temp_data2,'r');
            else
                plot(temp_data2,'b');
            end
            title(['sub(',num2str(count_sub),')','原波未差分, ',num2str(breathing_spectrum_scores(1,itemp))]);
            set(gca,'FontSize',8) % 设置坐标轴刻度字体名称，大小
        end
        %
        %% 复数域波形
        %         fig3=figure;
        % 画布划分
        %[ha, pos] = tight_subplot(行, 列, [上下间距,左右间距],[下边距,上边距 ], [左边距,右边距 ])
        %         [ha, pos] = tight_subplot(4, 5, [.05 .05], [.06 .06], [.06 .06]);
        %         count_temp = 0;
        for itemp = 1:1:20
            count_temp = count_temp + 1;
            count_app=ceil((count_temp-20)/4)*4+count_temp-20;
            temp_data1 = rawData_buffer(itemp,:);
            axes(ha(count_app));
            if(breathing_spectrum_scores(1,itemp)==max_breathing_scores)
                plot(temp_data1,'r');
            else
                plot(temp_data1,'b');
            end
            title(['sub(',num2str(count_sub),')','复数域, ',num2str(breathing_spectrum_scores(1,itemp))]);
            set(gca,'FontSize',8) % 设置坐标轴刻度字体名称，大小
        end
        figure_fullscreen(fig1);
        %% 复数域数据进行归一化
        fig2 = figure;
        count_temp = 0;
        % 画布划分
        %[ha, pos] = tight_subplot(行, 列, [上下间距,左右间距],[下边距,上边距 ], [左边距,右边距 ])
        [ha, pos] = tight_subplot(5, 4*2, [.05 .02], [.06 .06], [.06 .06]);
        for itemp = 1:1:20
            count_temp = count_temp + 1;
            count_app=count_temp + 4*(floor((count_temp-1)/4));
            temp_data1 = rawData_buffer(itemp,:);
            temp_data2 = -unwrap(angle(temp_data1));
            %             temp_data2 = diff(temp_data2); % 差分
            %             subplot(4,5,count_temp);
            axes(ha(count_app));
            if(breathing_spectrum_scores_normalized(1,itemp)==max_breathing_scores_nor)
                plot(temp_data2,'r');
            else
                plot(temp_data2,'b');
            end
            title(['sub(',num2str(count_sub),')','原波未差分, ',num2str(breathing_spectrum_scores_normalized(1,itemp))]);
            set(gca,'FontSize',8) % 设置坐标轴刻度字体名称，大小
        end
        %
        %% 复数域波形
        %         fig3=figure;
        % 画布划分
        %[ha, pos] = tight_subplot(行, 列, [上下间距,左右间距],[下边距,上边距 ], [左边距,右边距 ])
        %         [ha, pos] = tight_subplot(4, 5, [.05 .05], [.06 .06], [.06 .06]);
        %         count_temp = 0;
        for itemp = 1:1:20
            count_temp = count_temp + 1;
            count_app=ceil((count_temp-20)/4)*4+count_temp-20;
            temp_data1 = rawData_buffer(itemp,:);
            phase_curve = -unwrap(angle(temp_data1));
            % 归一化
            temp_data1 = temp_data1 ./ abs(temp_data1); % way1
            %             temp_data1 = exp(-1j * phase_curve); % way2
            axes(ha(count_app));
            if(breathing_spectrum_scores_normalized(1,itemp)==max_breathing_scores_nor)
                plot(temp_data1,'r');
            else
                plot(temp_data1,'b');
            end
            title(['sub(',num2str(count_sub),')','复数域, ',num2str(breathing_spectrum_scores_normalized(1,itemp))]);
            set(gca,'FontSize',8) % 设置坐标轴刻度字体名称，大小
        end
        figure_fullscreen(fig2);
        %
        %
        count_epoch=count_epoch+1;
        pause(0.0005);
        %%
        %         saveas(fig1,[num2str(count_sub),'_phase_',num2str(count_epoch)]);
        saveas(fig1,[num2str(count_sub),'_phase_val_',num2str(count_epoch)]); % 2合1
        %         saveas(fig1,[num2str(count_sub),'_diff_val_',num2str(count_epoch)]); % 2合1
        
        %
        %         saveas(fig2,[num2str(count_sub),'_diff_',num2str(count_epoch)]);
        %         saveas(fig3,[num2str(count_sub),'_val_',num2str(count_epoch)]);
        %         close all;
        disp('debug');
    end
end
disp([subName, ' is done! 序号: ', num2str(count_sub),', 阶段号: ', num2str(stage_number)]);
end