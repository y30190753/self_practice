function [ breathCurve ] = acquire_breathCurve_by_another_way( dataIn )
%                                              利用矢量相位差 求取 呼吸波形
% dataIn: 1 * 512
%%
%% 相位的另一种形式
pre_dealed_data = dataIn(1,2);
cur_dealed_data = pre_dealed_data;
new_coordinate_record = zeros(1, size(dataIn,2)); % 最终数据
%
len_data_buffer = 50; data_buffer = zeros(1,len_data_buffer); % 数据缓存器
data_buffer(1,end-1) = dataIn(1,1);
data_buffer(1,end) = dataIn(1,2);
polarity_change_flag = 0;
polarity_index = 1; % 初始状态不确定
polarity_record = zeros(1, size(dataIn,2));
vector_pre = 0;
vector_cur = 0;
%
for cur_data_index = 3:1:size(dataIn,2)
    cur_dealed_data = dataIn(1,cur_data_index); % 当前数据
    % 存缓存器
    data_buffer = circshift(data_buffer,[0,-1]);
    data_buffer(1,end) = cur_dealed_data;
    % 计算开始阶段的向量
    if(cur_data_index == 3)
        % _temp 表示原始数据计算出来的向量
        vector_cur_temp = data_buffer(1,end) - data_buffer(1,end-1);    % 当前点的向量
        vector_pre_temp = data_buffer(1,end-1) - data_buffer(1,end-2);  % 上一点的向量
        vector_pre = vector_pre_temp * polarity_index;
        vector_cur = vector_cur_temp * polarity_index;
        polarity_change_flag = convexity_concavity_identification(vector_cur, vector_pre); % 1 是改变了
        if(polarity_change_flag == 1) % 极性发生了改变
            polarity_index = -1 * polarity_index; % 极性指数
        end
        vector_cur = polarity_index * vector_cur_temp;
        new_coordinate_record(1,cur_data_index-1) = vector_cur;
        vector_pre = vector_cur;
    else % 之后的数据
        vector_cur_temp = data_buffer(1,end) - data_buffer(1,end-1);
        vector_cur = vector_cur_temp * polarity_index; % 含有极性
        polarity_change_flag = convexity_concavity_identification(vector_cur, vector_pre); % 1 是改变了
        if(polarity_change_flag == 1) % 极性发生了改变
            polarity_index = -1 * polarity_index; % 极性指数
        end
        vector_cur = vector_cur_temp * polarity_index; % 更新当前向量
        new_coordinate_record(1,cur_data_index-1) = vector_cur;
        vector_pre = vector_cur;
    end
    %
    %     new_coordinate_record(1,cur_data_index-1) = vector_cur;
    polarity_record(1,cur_data_index-2) = polarity_index;
end
breathCurve = -unwrap(angle(new_coordinate_record));

end

