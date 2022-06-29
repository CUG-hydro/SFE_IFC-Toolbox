function test_synthesis20200715()

clc
clear
load mopex423data.mat

old_corr = zeros(423, 3);
new_corr = zeros(423, 3);
old_corr_log = zeros(423, 3);
new_corr_log = zeros(423, 3);

% \theta denotes the interval time between two consecutive peaks (days), Eq. 1

for sta = 1:423
    disp(sta)
    %% Independence criteria and maximum days of start and end flood
    area = alldata_selected423{sta, 4}; % km^2
    interval = floor(5 + log(area / 1.609^2)); % min interval
    F_B_time = floor(1.5 * (5 + log(area / 1.609^2))); % max interval
    
    %% Select peaks serise
    % R2Q, mm/d to m/s 
    % 4 variables: year, month, day, Q
    flow_discharge = [alldata_selected423{sta, 6}(:, [1 2 3]), alldata_selected423{sta, 6}(:, 6) * area * 10^3/24/3600];
    flow_s_sort = sort(flow_discharge(:, 4));
    f_threshold = flow_s_sort(floor(length(flow_discharge(:, 4)) * 80/100)); % 80% quantile as the initial threshold
    peaks_datenum = selectpeaks(flow_discharge, f_threshold, interval);
    peaks_serise = peaks_datenum(:, 2);

    
    thre_can = sort(unique(peaks_serise)); % threshold candidate
    years = size(flow_discharge, 1) / 365;
    temp = length(thre_can);

    if temp <= 25
        temp = 26;
    end

    for num = 1:temp - 20
        peaks_serise1 = peaks_serise(peaks_serise > thre_can(num));
        x = (peaks_serise1 - thre_can(num))
        gpdist = fitdist(x, 'gp');
        [h, p, ad_sta, ~] = adtest(x, 'Distribution', gpdist);
        p_value(num, 1) = p;
        ad_value(num, 1) = ad_sta;
        k_shape(num, 1) = gpdist.k;
        ratio(num, 1) = length(peaks_serise1) / years; % 每年峰值的个数
    end

    index1_5 = find(ratio > 1.2 & ratio < 5); %determine PPY range
    thr_ad_ppy = [thre_can(index1_5), ad_value(index1_5), ratio(index1_5), p_value(index1_5)];
    %     [M,I]=min(abs(thr_ad_ppy(:,3)-3));
    [M, I] = min(thr_ad_ppy(:, 2));

    if isempty(I); I = 1; end
    selected_thr(sta, 1:4) = [thre_can(1), ad_value(1), ratio(1), p_value(1)];

    try
        peaks_datenum = selectpeaks(flow_discharge, selected_thr(sta, 1), interval);
    catch
        p_value = [];
        ad_value = [];
        ratio = [];
        continue
    end

    p_value = [];
    ad_value = [];
    ratio = [];
    [s_e_date_q, dura, ~] = starenddate(flow_discharge, peaks_datenum, F_B_time);
    count_jj = [];

    for jj = 1:size(peaks_datenum, 1)
        if 3/4 * peaks_datenum(jj, 2) <= max(s_e_date_q(jj, [4 8]), [], 2)
            count_jj(jj, 1) = jj;
        end
    end

    count_jj(count_jj == 0) = [];
    s_e_date_q(count_jj, :) = [];
    dura(count_jj, :) = [];
    [oldflood_p, oldflood_f, T1] = oldfloodevents(flow_discharge, s_e_date_q(:, [1 2 3 5 6 7]), dura(:, 1));

    deepbaseflow = Deepbaseflow(flow_discharge);
    deepbaseflow(deepbaseflow <= 1) = 1; %20200409revised

    for ii = 0:30
        [rece_series, ~] = extract_recession(flow_discharge, 5 + ii, 1);
        if length(rece_series) < 100
            %                     ii
            break
        end

    end

    Min_durationEditField = 5 + ii;
    [rece_series, dura] = extract_recession(flow_discharge, Min_durationEditField, 1);
    [posi, sort_num, MRCdata, COR_MRCdata] = MRC_auto(rece_series, dura);
    [fit1, fit2, fit3, gof1, gof2, gof3] = CreateFit(COR_MRCdata);

    para(1) = fit1.a;

    x = (COR_MRCdata(1, 1):0.01:COR_MRCdata(end, 1))';

    RE_COR_MRCdata = [x, fit1(x)];

    [newflood_p, newflood_f, T2] = FloodCharacteristics(flow_discharge, s_e_date_q, RE_COR_MRCdata);

    %% 计算洪水特征间的相关关系
    oldtemp = corr(oldflood_f);
    newtemp = corr(newflood_f);
    old_corr(sta, :) = [oldtemp(1, 2), oldtemp(1, 3), oldtemp(2, 3)];
    new_corr(sta, :) = [newtemp(1, 2), newtemp(1, 3), newtemp(2, 3)];

    oldtemplog = corr([log10(oldflood_f(:, 1:2)), oldflood_f(:, 3)]);
    newtemplog = corr([log10(newflood_f(:, 1:2)), newflood_f(:, 3)]);
    old_corr_log(sta, :) = [oldtemplog(1, 2), oldtemplog(1, 3), oldtemplog(2, 3)];
    new_corr_log(sta, :) = [newtemplog(1, 2), newtemplog(1, 3), newtemplog(2, 3)];

end

%% figure preview
%     yyaxis left
%     plot(thre_can(index1_5),ad_value(index1_5),'LineWidth',2); hold on
%     plot(selected_thr(1,1),selected_thr(1,2),'rp','MarkerSize',10,'MarkerFaceColor',[0.7 0 0.7])
%     % plot(thre_can(index1_5),1-p_value(index1_5))
%     ylabel 'AD statistic'
%     yyaxis right
%     plot(thre_can(index1_5),ratio(index1_5),'LineWidth',2)
%      ylabel 'Peaks per year';
%
%        xlabel 'Threshold(m^3/s)'
%     grid on
%         box on
%          set(gca,'FontSize',18,'LineWidth',2);
%     print(gcf,'-r600','-dpng',[num2str(sta),'.png']);
%     close

% out_data=[cell2mat(alldata_selected423(:,2:4)),mean_pre,selected_thr];% longitude ,latitude ,area ( cubic Km),mean annual precipitation ,Optimal threshold,AD statistics ,PPY

end
