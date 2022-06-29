function [peaks_datenum] = selectpeaks(Q, t, threshold, interval)
    %selectpeaks 基于POT方法选择洪峰并进行独立性判断
    %   主要是对给定的百分位阈值与独立判断条件的间隔时间进行洪峰筛选与函数POTpeak基本相同
    %% 洪峰的初步筛选
    % Q = data(:, 4);
    % 两次三点滑动平均
    % for mm=1:1
    % Q=smoothdata(Q,'movmean',3);
    % end
    % datevec = data(:, 1:3);
    % t = t(datevec);
    flow_s_sort = sort(Q);
    %20200616gai将百分比输入换为阈值直接输入
    f_threshold = threshold;
    %f_threshold=flow_s_sort(floor(length(Q)*Quantile/100));%以50%分位选择阈值
    % index_exce=find(Q>f_threshold);
    %挑选洪峰应用findpeaks
    [peaks, t_peaks] = findpeaks(Q, t, 'MinPeakHeight', f_threshold);

    %% 检验洪峰独立性，20191211改
    t_diff = diff(t_peaks);
    
    % 挑出独立的洪峰间隔, Eq. 3
    event_count = zeros(length(peaks) - 1, 1);
    for ii = 1:(length(peaks) - 1)
        index_1 = find(t == t_peaks(ii));
        index_2 = find(t == t_peaks(ii + 1));
        span = index_1:index_2;
        midqmin = min(Q(span));
        Qmin66 = 3/4 * min(peaks(ii), peaks(ii + 1));

        if t_diff(ii) >= interval && midqmin < Qmin66
            event_count(ii, 1) = ii; 
        end
    end

    %% 针对超过双峰的洪水作出判断
    for kk = 1:10
        c_num = 0;
        ind1 = find(event_count == 0);  % 挨的过近的洪峰
        ind2 = find(diff(ind1) > 1);    % 划分为不同的场次
        ind2 = [0; ind2; length(ind1)]; % 

        for nn = 1:size(ind2, 1) - 1

            i_begin = ind1(ind2(nn) + 1);
            i_end = ind1(ind2(nn + 1));

            datespan = (t_peaks(i_begin):t_peaks(i_end + 1))';
            eventflow = Q(t >= datespan(1) & t <= datespan(end));
            minvalue = min(eventflow);
            minvalue_date = datespan(eventflow == minvalue);
            % 如果是紧挨着
            % 合并紧挨着的low peaks
            if i_end - i_begin == 1
                % 这一步是在做什么？
                if length(datespan) >= interval && minvalue < 3/4 * min(peaks(i_begin), peaks(i_end + 1))
                    c_num = c_num + 1;
                    if minvalue_date < t_peaks(i_begin + 1)
                        event_count(i_begin) = i_begin;
                    else
                        event_count(i_begin + 1) = i_begin + 1;
                    end
                end
            elseif i_end - i_begin >= 2
                index_min = find(t_peaks(i_begin:i_end + 1) < minvalue_date(1));
                temp = peaks(i_begin:i_end + 1);
                
                if minvalue < 3/4 * max(temp)
                    c_num = c_num + 1;
                    event_count(i_begin + index_min(end) - 1) = i_begin + index_min(end) - 1;
                end
            end
        end

        if c_num == 0
            break
        end
    end

    event_count(event_count == 0) = [];
    event_p_info = zeros(length(event_count) + 1, 3);
    event_p_datenuminfo = zeros(length(event_count) + 1, 3);

    for jj = 0:length(event_count)

        if jj == 0
            [event_peak, loc] = max(peaks(1:event_count(jj + 1)));
            e_p_star = peaks(1);
            e_p_end = peaks(event_count(jj + 1));
            e_p_star_dnum = t_peaks(1);
            e_p_end_dnum = t_peaks(event_count(jj + 1));
            e_p_datenum = t_peaks(loc);
        elseif jj == length(event_count)
            [event_peak, loc] = max(peaks(event_count(jj) + 1:end));
            e_p_star = peaks(event_count(jj) + 1);
            e_p_end = peaks(end);
            e_p_star_dnum = t_peaks(event_count(jj) + 1);
            e_p_end_dnum = t_peaks((end));
            e_p_datenum = t_peaks(event_count(jj) + loc);
        else
            [event_peak, loc] = max(peaks(event_count(jj) + 1:event_count(jj + 1)));
            e_p_star = peaks(event_count(jj) + 1);
            e_p_end = peaks(event_count(jj + 1));
            e_p_star_dnum = t_peaks(event_count(jj) + 1);
            e_p_end_dnum = t_peaks(event_count(jj + 1));
            e_p_datenum = t_peaks(event_count(jj) + loc);
        end

        event_p_info(jj + 1, :) = [e_p_star, event_peak, e_p_end];
        event_p_datenuminfo(jj + 1, :) = [e_p_star_dnum, e_p_datenum, e_p_end_dnum];
    end

    peaks_datenum = [event_p_info, event_p_datenuminfo];
    % % 检验独立性，去掉不独立的洪峰，给定一个最小间隔时间
    % for jj=1:10
    %     peaks_1=peaks;
    %     peaksdatenum_1=t_peaks;
    %     for ii=1:length(peaks)-1
    %         index_1=find(t==t_peaks(ii));
    %         index_2=find(t==t_peaks(ii+1));
    %         span=index_1:index_2;
    %         midq=Q(span);
    %         Qmin50=3/4*min(peaks(ii),peaks(ii+1));
    %         if t_peaks(ii+1)-peaksdatenum(ii)<=interval||isempty(find(midq<=Qmin50, 1))==1
    %             if peaks(ii+1)>peaks(ii)
    %                 peaks_1(ii)=0;peaksdatenum_1(ii)=0;
    %             else
    %                 peaks_1(ii+1)=0;peaksdatenum_1(ii+1)=0;%将小于的峰值赋值为0
    %             end
    %         end
    %     end
    %     peaks=peaks_1(peaks_1>0);
    %     t_peaks=peaksdatenum_1(peaksdatenum_1>0);
    %     if length(peaks)==length(peaks_1)
    %         break
    %     end
    % end
    % peaksdatevec=datevec(t_peaks);
    % date_peaks=[peaksdatevec(:,1:3),peaks];
end

% mm/d, Q/m3
