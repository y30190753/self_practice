function [ weights_final ] = weights_opt_alternative( dataIn1,dataIn2,flag,num_comopute )
%%
%    dataIn1：原始数据， dataIn2：上次迭代的权重
%    flag: 1--随机初始值， num_comopute：优化次数
%
%%
global param_objectiveFun;
global param_model;
global param_fixed_phase;
%%
[n_channel, n_points] = size(dataIn1);
param_objectiveFun = dataIn1;
param_model = dataIn2;
%
weights_maxScore_record = [];
loss_max_record = +2.1;
N_computeTimes = num_comopute; % 随机点数
pre_param_part1 = [];
pre_param_part2 = [];
%
if(flag == 1) % 随机初始权重
    count_temp=0;
    for itemp=1:1:num_comopute
        count_temp=count_temp+1;
        if(count_temp > 1) % 继承权重
            initial_params_part1 = pre_param_part1;
            initial_params_part2 = pre_param_part2;
        else
            initial_params_part1 = rand(1,1);
            initial_params_part2 = rand(size(dataIn1,1),1);
        end
        
        % 交替优化_part1 -- 实数部分为1，优化相位
        lower_bound = [-pi];
        upper_bound = [ pi];
        [phase_opted, loss_1] = fmincon('objective_function_offline_fixedRealPart',...
            initial_params_part1,...
            [],[],[],[],...
            lower_bound,...
            upper_bound,...
            'constraint_conditoin_offline_fixedRealPart');
        param_fixed_phase = phase_opted;
        pre_param_part1 = param_fixed_phase;
        % 得到 优化的相位：param_fixed_phase
        %
        % 交替优化_part2 -- 固定相位，优化实数部分
        % 是否需要对实数部分做大小约束？ 需要，否则会出现负值。
        %         lower_bound = [zeros(size(dataIn1,1),1); zeros(size(dataIn1,1),1) - pi]; % [0;0;...;-pi;-pi];
        %         upper_bound = [ones(size(dataIn1,1),1);  zeros(size(dataIn1,1),1) + pi]; % [1;1;... ;pi; pi];
        lower_bound = zeros(size(dataIn1,1),1) + 0;
        upper_bound = zeros(size(dataIn1,1),1) + 1;
        [weights_real_part_opted, loss2] = fmincon('objective_function_offline_fixedPhase',...
            initial_params_part2,...
            [],[],[],[],...
            lower_bound,...
            upper_bound,...
            'constraint_condition_offline_fixedPhase');
        pre_param_part2 = weights_real_part_opted;
        % 权重整理，写成 2Nc * 1 形式
        weights_reshape = [weights_real_part_opted; ones(n_channel, 1) * param_fixed_phase];
        % 记录最优的
        if(loss_max_record>loss2) % loss2为最终损失
            loss_max_record=loss2;
            weights_maxScore_record=weights_reshape;
        end
    end
elseif(flag == 2) % 继承初始权重
    % ......
end
weights_final = weights_maxScore_record;
end

