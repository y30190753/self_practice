function figure_fullscreen( h )
%
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame'); % �ر���صľ�����ʾ����Ϊ�����˷ǹ����ӿڣ�
jFrame = get(h,'JavaFrame'); % ��ȡ�ײ� Java �ṹ��ؾ����
pause(0.1); % �� Win 10��Matlab 2017b �����²���ͣ�ٻᱨ Java �ײ���󡣸��˸�����Ҫ���Խ���ʵ����֤
set(jFrame,'Maximized',1); %���������Ϊ�棨0 Ϊ�٣�
pause(0.1); % ����ʵ���з��������ͣ�٣����ڿ����������仯������ȡ�Ĵ��ڴ�С����ԭ���ĳߴ硣���˸�����Ҫ���Խ���ʵ����֤
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame'); % ����ؾ�������
end

