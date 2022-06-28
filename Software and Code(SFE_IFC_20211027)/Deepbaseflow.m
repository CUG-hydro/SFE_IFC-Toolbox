function [baseflow] = Deepbaseflow(flowdata)
%Deepbaseflow 计算深层基流
%   此处显示详细说明

data=flowdata(:,4);
data1=reshape(data(1:floor(length(data)/365)*365),365,[]);
data2=movmean(data1,91);
data3=min(data2(46:end-45,:));
baseflow=mean(data3(:));
end

