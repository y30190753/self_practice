function [ f ] = objective_function_offline( x )
%%
global param_objectiveFun;
global lambda;
%%
data_combined = transpose(x) * param_objectiveFun;

% data_combined = [x(1)*exp(1j*x(5)),x(2)*exp(1j*x(6)),x(3)*exp(1j*x(7)),x(4)*exp(1j*x(8))] * param_objectiveFun;

% data_combined_phaseCurve = -unwrap(angle(data_combined));
data_combined_phaseCurve = -(angle(data_combined));

data_combined_ideal = exp(-1j * data_combined_phaseCurve);

% part1 = corr(transpose(real(data_combined)), transpose(real(data_combined_ideal)));
% part2 = corr(transpose(imag(data_combined)), transpose(imag(data_combined_ideal)));
part3 = corr(data_combined.', data_combined_ideal.');

%%
% f= -(part1+part2); % 求极小值
f = -abs(part3);
end