function [ breathCurve ] = acquire_breathCurve_by_another_way( dataIn )
%                                              ����ʸ����λ�� ��ȡ ��������
% dataIn: 1 * 512
%%
%% ��λ����һ����ʽ
pre_dealed_data = dataIn(1,2);
cur_dealed_data = pre_dealed_data;
new_coordinate_record = zeros(1, size(dataIn,2)); % ��������
%
len_data_buffer = 50; data_buffer = zeros(1,len_data_buffer); % ���ݻ�����
data_buffer(1,end-1) = dataIn(1,1);
data_buffer(1,end) = dataIn(1,2);
polarity_change_flag = 0;
polarity_index = 1; % ��ʼ״̬��ȷ��
polarity_record = zeros(1, size(dataIn,2));
vector_pre = 0;
vector_cur = 0;
%
for cur_data_index = 3:1:size(dataIn,2)
    cur_dealed_data = dataIn(1,cur_data_index); % ��ǰ����
    % �滺����
    data_buffer = circshift(data_buffer,[0,-1]);
    data_buffer(1,end) = cur_dealed_data;
    % ���㿪ʼ�׶ε�����
    if(cur_data_index == 3)
        % _temp ��ʾԭʼ���ݼ������������
        vector_cur_temp = data_buffer(1,end) - data_buffer(1,end-1);    % ��ǰ�������
        vector_pre_temp = data_buffer(1,end-1) - data_buffer(1,end-2);  % ��һ�������
        vector_pre = vector_pre_temp * polarity_index;
        vector_cur = vector_cur_temp * polarity_index;
        polarity_change_flag = convexity_concavity_identification(vector_cur, vector_pre); % 1 �Ǹı���
        if(polarity_change_flag == 1) % ���Է����˸ı�
            polarity_index = -1 * polarity_index; % ����ָ��
        end
        vector_cur = polarity_index * vector_cur_temp;
        new_coordinate_record(1,cur_data_index-1) = vector_cur;
        vector_pre = vector_cur;
    else % ֮�������
        vector_cur_temp = data_buffer(1,end) - data_buffer(1,end-1);
        vector_cur = vector_cur_temp * polarity_index; % ���м���
        polarity_change_flag = convexity_concavity_identification(vector_cur, vector_pre); % 1 �Ǹı���
        if(polarity_change_flag == 1) % ���Է����˸ı�
            polarity_index = -1 * polarity_index; % ����ָ��
        end
        vector_cur = vector_cur_temp * polarity_index; % ���µ�ǰ����
        new_coordinate_record(1,cur_data_index-1) = vector_cur;
        vector_pre = vector_cur;
    end
    %
    %     new_coordinate_record(1,cur_data_index-1) = vector_cur;
    polarity_record(1,cur_data_index-2) = polarity_index;
end
breathCurve = -unwrap(angle(new_coordinate_record));

end

