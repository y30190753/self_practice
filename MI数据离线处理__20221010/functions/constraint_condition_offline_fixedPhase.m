function [ h,g ] = constraint_condition_offline_fixedPhase( x )
%%
global param_objectiveFun;
%%
% �����Բ���ʽԼ��
h=[];
% �����Ե�ʽԼ��
g=[sum(x(1:size(param_objectiveFun,1), 1).^2) - 1];
% g=[];
end