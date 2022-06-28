function [peaks_datenum] = AMS_sample(flowdata)

    data = flowdata(:, 4);
    datenum_1 = datenum(flowdata(:, 1:3));
    data1 = reshape(data(1:floor(length(data) / 365) * 365), 365, []);
    datenum_2 = reshape(datenum_1(1:floor(length(data) / 365) * 365), 365, []);
    [M, I] = max(data1, [], 1, 'linear');
    out_datenum = (datenum_2(I))';
    peaks_datenum = [M', M', M', out_datenum, out_datenum, out_datenum];
end
