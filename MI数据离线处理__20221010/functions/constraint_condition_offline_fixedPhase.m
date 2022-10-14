function [ h,g ] = constraint_condition_offline_fixedPhase( x )
%%
global param_objectiveFun;
%%
% 非线性不等式约束
h=[];
% 非线性等式约束
g=[sum(x(1:size(param_objectiveFun,1), 1).^2) - 1];
% g=[];
end