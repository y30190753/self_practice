function figure_fullscreen( h )
%
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame'); % 关闭相关的警告提示（因为调用了非公开接口）
jFrame = get(h,'JavaFrame'); % 获取底层 Java 结构相关句柄吧
pause(0.1); % 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
set(jFrame,'Maximized',1); %设置其最大化为真（0 为假）
pause(0.1); % 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame'); % 打开相关警告设置
end

