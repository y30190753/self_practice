function algorithmForWaveOpt( GUIdata, fileLoadPath, subName )
%%
%         改成离线数据处理
%
%%
global count_sub;
% count_sub = count_sub + 1;
%% 15 * 4
timeIndex_34_of15Sub = [6190,15190,16000,25200;
    8770,14700,14850,21770;
    4001,13000,13191,22200;
    6251,15250,15500,25800;
    6301,15300,15560,26170;
    6218,15295,15435,24708;
    966,10029,10097,19519;
    7600,17000,17400,17800;
    8100,27100,27300,36800;
    5700,14700,14950,30600;
    5050,14050,14360,27000;
    18300,27400,27530,36820;
    11710,20740,20890,30320;
    18000,27000,27260,36600;
    1100,19000,19000,19000;];
%
rawData_stage3 = timeIndex_34_of15Sub(count_sub,1):timeIndex_34_of15Sub(count_sub,2); % 第三阶段的数据
rawData_stage4 = timeIndex_34_of15Sub(count_sub,3):timeIndex_34_of15Sub(count_sub,4); % 第四阶段的数据
%
algorithmForWaveopt_part(GUIdata, rawData_stage3, subName, 3);
algorithmForWaveopt_part(GUIdata, rawData_stage4, subName, 4);
%
disp([subName, ' is done! 序号: ', num2str(count_sub)]);
end

