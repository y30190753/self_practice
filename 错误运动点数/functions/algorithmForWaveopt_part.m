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
%% ���㷨1\2
% [outHeart, outPhase] = old_optAlgorithm(GUIdata, subName, selected_index);
disp('���㷨������ϡ�');
%%
rangeProf = reshape(GUIdata.rangeProf,22,4,size(GUIdata.numFrame,2));
rangeBin_scope = 10:1:14; % rangeBin��ѡ��Χ
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
% ������
len_rawData_buffer = 512; % ����Ƭ�γ���
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
        % �����Ƶ��
        breathing_spectrum_scores = zeros(1,size(rawData_buffer,1));
        breathing_spectrum_scores_normalized = zeros(1,size(rawData_buffer,1)); % ���������ݽ��й�һ����
        absSpectrum_record=zeros(size(rawData_buffer,1),len_rawData_buffer);
        sum_3_part_record = zeros(size(rawData_buffer,1),3);
        for itemp=1:1:size(rawData_buffer,1)
            temp_data = -unwrap(angle(rawData_buffer(itemp,:)));
            temp_data2 = rawData_buffer(itemp,:); % ԭʼ����������
            temp_data = diff(temp_data); % ���
            %             breathing_spectrum_scores(1,itemp)=computeConfidenceMetric2(temp_data,50,len_rawData_buffer);
            % ���������ݵĴ���Բ�ĵ���
            breathing_spectrum_scores(1,itemp)=breathingScore_circCenterPosition(temp_data2);
            % ��һ��
            temp_data2_nor = temp_data2 ./ abs(temp_data2); % way1
            %             temp_data2_nor = exp(-1j * -unwrap(angle(rawData_buffer(itemp,:)))); % way2
            breathing_spectrum_scores_normalized(1, itemp) = breathingScore_circCenterPosition(temp_data2_nor);
        end
        max_breathing_scores = max(breathing_spectrum_scores);
        max_breathing_scores_nor = max(breathing_spectrum_scores_normalized);
        %
        %% ԭʼ���Σ�δ���
        fig1 = figure;count_temp = 0;
        % ��������
        %[ha, pos] = tight_subplot(��, ��, [���¼��,���Ҽ��],[�±߾�,�ϱ߾� ], [��߾�,�ұ߾� ])
        [ha, pos] = tight_subplot(5, 4*2, [.05 .02], [.06 .06], [.06 .06]);
        for itemp = 1:1:20
            count_temp = count_temp + 1;
            count_app=count_temp + 4*(floor((count_temp-1)/4));
            temp_data1 = rawData_buffer(itemp,:);
            temp_data2 = -unwrap(angle(temp_data1));
            %             temp_data2 = diff(temp_data2); % ���
            %             subplot(4,5,count_temp);
            axes(ha(count_app));
            if(breathing_spectrum_scores(1,itemp)==max_breathing_scores)
                plot(temp_data2,'r');
            else
                plot(temp_data2,'b');
            end
            title(['sub(',num2str(count_sub),')','ԭ��δ���, ',num2str(breathing_spectrum_scores(1,itemp))]);
            set(gca,'FontSize',8) % ����������̶��������ƣ���С
        end
        %
        %% ��������
        %         fig3=figure;
        % ��������
        %[ha, pos] = tight_subplot(��, ��, [���¼��,���Ҽ��],[�±߾�,�ϱ߾� ], [��߾�,�ұ߾� ])
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
            title(['sub(',num2str(count_sub),')','������, ',num2str(breathing_spectrum_scores(1,itemp))]);
            set(gca,'FontSize',8) % ����������̶��������ƣ���С
        end
        figure_fullscreen(fig1);
        %% ���������ݽ��й�һ��
        fig2 = figure;
        count_temp = 0;
        % ��������
        %[ha, pos] = tight_subplot(��, ��, [���¼��,���Ҽ��],[�±߾�,�ϱ߾� ], [��߾�,�ұ߾� ])
        [ha, pos] = tight_subplot(5, 4*2, [.05 .02], [.06 .06], [.06 .06]);
        for itemp = 1:1:20
            count_temp = count_temp + 1;
            count_app=count_temp + 4*(floor((count_temp-1)/4));
            temp_data1 = rawData_buffer(itemp,:);
            temp_data2 = -unwrap(angle(temp_data1));
            %             temp_data2 = diff(temp_data2); % ���
            %             subplot(4,5,count_temp);
            axes(ha(count_app));
            if(breathing_spectrum_scores_normalized(1,itemp)==max_breathing_scores_nor)
                plot(temp_data2,'r');
            else
                plot(temp_data2,'b');
            end
            title(['sub(',num2str(count_sub),')','ԭ��δ���, ',num2str(breathing_spectrum_scores_normalized(1,itemp))]);
            set(gca,'FontSize',8) % ����������̶��������ƣ���С
        end
        %
        %% ��������
        %         fig3=figure;
        % ��������
        %[ha, pos] = tight_subplot(��, ��, [���¼��,���Ҽ��],[�±߾�,�ϱ߾� ], [��߾�,�ұ߾� ])
        %         [ha, pos] = tight_subplot(4, 5, [.05 .05], [.06 .06], [.06 .06]);
        %         count_temp = 0;
        for itemp = 1:1:20
            count_temp = count_temp + 1;
            count_app=ceil((count_temp-20)/4)*4+count_temp-20;
            temp_data1 = rawData_buffer(itemp,:);
            phase_curve = -unwrap(angle(temp_data1));
            % ��һ��
            temp_data1 = temp_data1 ./ abs(temp_data1); % way1
            %             temp_data1 = exp(-1j * phase_curve); % way2
            axes(ha(count_app));
            if(breathing_spectrum_scores_normalized(1,itemp)==max_breathing_scores_nor)
                plot(temp_data1,'r');
            else
                plot(temp_data1,'b');
            end
            title(['sub(',num2str(count_sub),')','������, ',num2str(breathing_spectrum_scores_normalized(1,itemp))]);
            set(gca,'FontSize',8) % ����������̶��������ƣ���С
        end
        figure_fullscreen(fig2);
        %
        %
        count_epoch=count_epoch+1;
        pause(0.0005);
        %%
        %         saveas(fig1,[num2str(count_sub),'_phase_',num2str(count_epoch)]);
        saveas(fig1,[num2str(count_sub),'_phase_val_',num2str(count_epoch)]); % 2��1
        %         saveas(fig1,[num2str(count_sub),'_diff_val_',num2str(count_epoch)]); % 2��1
        
        %
        %         saveas(fig2,[num2str(count_sub),'_diff_',num2str(count_epoch)]);
        %         saveas(fig3,[num2str(count_sub),'_val_',num2str(count_epoch)]);
        %         close all;
        disp('debug');
    end
end
disp([subName, ' is done! ���: ', num2str(count_sub),', �׶κ�: ', num2str(stage_number)]);
end