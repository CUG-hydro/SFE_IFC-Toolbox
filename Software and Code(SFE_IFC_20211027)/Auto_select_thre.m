function [output_para, figuredata] = Auto_select_thre(flowdata, interval)
    %Auto_select_thre ����AD��������Ŷ��Զ�ѡȡ���POT��ֵ
    %   input flowdata  Ϊ�������ݣ���������������intervalΪ����������С�б�ʱ�䣬
    %   output_paraΪ���������ѡ�����ֵ��AD����ֵ���Ӧ��pֵ��ppy����figuredata����ֵ����ֵ�ı仯
    flow_s_sort = sort(flowdata(:, 4));
    f_threshold = flow_s_sort(floor(length(flowdata(:, 4)) * 95/100)); %��95 %��λѡ���ʼ��ֵ
    [peaks_datenum] = selectpeaks(flowdata, f_threshold, interval);
    peaks_serise = peaks_datenum(:, 2);
    thre_can = sort(unique(peaks_serise));
    years = size(flowdata, 1) / 365;
    temp = length(thre_can);
    temp(temp <= 25) = 26;

    for num = 1:temp - 25
        peaks_serise1 = peaks_serise(peaks_serise > thre_can(num));
        gpdist = fitdist((peaks_serise1 - thre_can(num)), 'gp');
        [~, p, ad_sta, ~] = adtest((peaks_serise1 - thre_can(num)), 'Distribution', gpdist);
        p_value(num, 1) = p;
        ad_value(num, 1) = ad_sta;
        k_shape(num, 1) = gpdist.k;
        ratio(num, 1) = length(peaks_serise1) / years;
    end

    index1_5 = find(ratio > 1.2 & ratio < 5); %determine PPY range
    thr_ad_ppy_p = [thre_can(index1_5), ad_value(index1_5), ratio(index1_5), p_value(index1_5)];
    [~, I] = min(thr_ad_ppy_p(:, 2));

    if isempty(I) == 1
        output_para = [thre_can(1), ad_value(1), ratio(1), p_value(1)];
    else
        output_para = thr_ad_ppy_p(I, :);
    end

    figuredata = thr_ad_ppy_p;
end
