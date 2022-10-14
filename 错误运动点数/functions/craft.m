%
clear all;
close all;
clc;

%
data1 = 10+40j;
data2 = 20+20j;
data3 = 15+35j;


% data1 = data1 / abs(data1);
% data2 = data2 / abs(data2);
% data3 = data3 / abs(data3);



data_list = [data1,data2,data3];

figure; plot(data_list);
hold on; plot(real(data_list), imag(data_list), 'ro');
hold on; plot(0,0,'o');
breathingScore_circCenterPosition(data_list);





