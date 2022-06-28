function [baseflow] = Deepbaseflow(flowdata)
    data = flowdata(:, 4);
    data1 = reshape(data(1:floor(length(data) / 365) * 365), 365, []);
    data2 = movmean(data1, 91);
    data3 = min(data2(46:end - 45, :));
    baseflow = mean(data3(:));
end
