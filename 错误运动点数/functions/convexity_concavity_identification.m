function [ flag ] = convexity_concavity_identification( rawData1, rawData2 )
%%
%    ���룬��ǰ
%    ������нǴ���90�ȣ�������Ա仯
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

