function [output_para] = validatation_poisson_ks(date_events, flowdata)
    %validatation_poisson_ks 使用ks检验验证序列的泊松假设
    %   date_events 事件的洪峰序列及对应时间，flowdata 为流量序列
    years_count = tabulate(flowdata(:, 1));
    years_num = years_count(years_count(:, 2) >= 365, 1);
    peaks_year = date_events(ismember(date_events(:, 1), years_num), 1);
    peaksnum_count = tabulate(peaks_year);
    peaksnum_count(peaksnum_count(:, 2) == 0 & peaksnum_count(:, 1) < years_num(1), :) = [];

    %      poidist=fitdist(peaksnum_count(:,2),'Poisson') ;
    %      [h,p,ks,~]=kstest(peaksnum_count(:,2),'CDF',poidist);
    c_var = std(peaksnum_count(:, 2))^2;
    c_mean = mean(peaksnum_count(:, 2));
    I_d = sqrt(c_var / c_mean);
    confidence_interval = [chi2inv(0.025, length(peaksnum_count) - 1), chi2inv(0.975, length(peaksnum_count) - 1)];
    Id5_95 = confidence_interval / (length(peaksnum_count) - 1);
    %% 计算Kendall tau 相关
    tau = corr([date_events(1:end - 1, 4), date_events(2:end, 4)], 'Type', 'Kendall');
    n = length(date_events);
    mean_s = -2 / (3 * (n - 1));
    var_s = (20 * n^3 - 74 * n^2 + 54 * n + 148) / (45 * (n - 1)^2 * (n - 2)^2);
    tau95up = mean_s + 1.65 * sqrt(var_s);
    tau95down = mean_s - 1.65 * sqrt(var_s);

    output_para = [I_d, Id5_95, tau(1, 2), tau95up, tau95down];
end
