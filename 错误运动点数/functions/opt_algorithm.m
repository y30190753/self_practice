function [ weights, objective_value ] = opt_algorithm( dataIn1,dataIn2,flag,num_comopute )
%%
%   输入：data_buffer的部分数据, 5 * Np;  flag==1是初始优化；flag==2是继承参数
%   输出：模型参数，weights
%
%%
global param_objectiveFun;
global param_model;
%%
[n_channel, n_points] = size(dataIn1);
param_objectiveFun = dataIn1;
param_model = dataIn2;

lower_bound = zeros(20, 1);
upper_bound = ones(20, 1);

x_max_record = [];
y_max_record = +2.1;

if(flag == 1)
    N_computeTimes = num_comopute; % 随机点数
    count_temp = 0;
    for temp=1:1:N_computeTimes
        count_temp = count_temp + 1;
        [x,y] = fmincon('objective_function_offline',rand(20,1),[],[],[],[],lower_bound,upper_bound,'constraint_condition_offline');
        %
        if(y_max_record>y)
            y_max_record=y;
            x_max_record=x;
        end
    end
    x=x_max_record;
    y=y_max_record;
    %
    weights = x;
    objective_value = -1 * y;
elseif(flag == 2)
    [x,y] = fmincon('objective_function_offline',param_model,[],[],[],[],lower_bound,upper_bound,'constraint_condition_offline');
    weights = x;
    objective_value = -1 * y;
end

end

