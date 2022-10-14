function [ angle_unwrapped ] = angle_unwrap( cur_angle, pre_angle )
%% 用于相位解缠绕
if(abs(cur_angle - pre_angle) > pi)
    if(cur_angle - pre_angle > pi)
        while(cur_angle - pre_angle > pi)
            cur_angle = cur_angle - 2*pi;
        end
    else
        while(cur_angle - pre_angle < -1*pi)
            cur_angle = cur_angle + 2*pi;
        end
    end
    angle_unwrapped = cur_angle;
else
    angle_unwrapped = cur_angle;
end
end

