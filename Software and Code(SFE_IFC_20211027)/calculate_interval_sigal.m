function [interval] = calculate_interval_sigal(data)
%calculate_interval �����������С���ʱ��
%   �˴���ʾ��ϸ˵��
% �������㻬��ƽ��
data1=data(:,4);
for mm=1:3
    data1=smoothdata(data1,'movmean',3);
end
flow_s_sort=sort(data1);
f_threshold=flow_s_sort(floor(length(data1)*80/100));%��75%��λѡ����ֵ%80% 20200409��
[~,loc,w,~]=findpeaks(data1,'MinPeakHeight',f_threshold,'Annotate','extents');
interval(1)=ceil(1.3*median(w));
if interval(1)<5
    interval(1)=5;
end
temp=diff(loc);
% temp(temp>quantile(temp,1-0.5*size(data,1)/(365*length(loc))))=0;
% interval(1)=ceil(quantile(temp,0.50));
interval(2)=ceil(quantile(w,0.85));
end

