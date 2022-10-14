function [score] = breathingScore_circCenterPosition(dataIn)
%%
%     score 错误的点数
%     dataIn: 1 * n_points
%%
num_counter = 0;
n_samples = size(dataIn, 2);
for i = 3:1:n_samples
    data_point1 = dataIn(1, i-2);
    data_point2 = dataIn(1, i-1);
    data_point3 = dataIn(1, i-0);
    flag = circle_center_wrong([data_point1,data_point2,data_point3]);
    if (flag == -1)
        num_counter = num_counter + 1;
    end
end
score = -1 * num_counter;
end

function [flag] = circle_center_wrong(point_list)
%% 判断当前点的运动时，圆心是否处于错误的一侧
%  输入: point_list： 1 * 3 复数list
%  输出：flag == 1 没有出错，flag == -1 出错

flag_clockwise = is_cloclwise(point_list); % 1表示顺时针运动
circle_center_position = is_cloclwise([point_list(1,end-2), point_list(1,end-1), 0]); % 1 表示在右侧

% 判断两个相邻的向量之间的夹角是否为钝角
vector1_temp = point_list(1,2) - point_list(1,1);
vector2_temp = point_list(1,3) - point_list(1,2);
vector1 = [real(vector1_temp), imag(vector1_temp)];
vector2 = [real(vector2_temp), imag(vector2_temp)];
%
angle_90 = sum(vector1 .* vector2);
if(angle_90 > 0) % 锐角
    flag = flag_clockwise * circle_center_position;
else % 钝角的判断
    flag = 1;
end
end

function [ flag ] = is_cloclwise(point_list)
%%
% 输入: point_list： 1 * 3 复数list
% 输出： flag == 1 表示是顺时针运动，flag == -1 表示是逆时针运动，flag == 0 表示三点共线
%%
x1 = real(point_list(1,1)); y1 = imag(point_list(1,1));
x2 = real(point_list(1,2)); y2 = imag(point_list(1,2));
x3 = real(point_list(1,3)); y3 = imag(point_list(1,3));
%
reference = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
if(reference > 0)
    flag = -1;
elseif(reference < 0)
    flag = 1;
else
    flag = 0;
end
end

