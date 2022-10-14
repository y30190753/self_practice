function [ flag ] = convexity_concavity_identification( rawData1, rawData2 )
%%
%    输入，当前
%    输出：夹角大于90度，输出极性变化
%%
vector1 = [real(rawData1);imag(rawData1)];
vector2 = [real(rawData2);imag(rawData2)];
%
if( sign(dot(vector1, vector2)) == -1)
    flag = 1;
else
    flag = 0;
end
end

