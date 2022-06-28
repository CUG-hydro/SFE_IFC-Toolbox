function [rece_series, dura] = extract_recession(date_flow, min_dura, mov_mean)
    %extract recession
    raw_flow_s = date_flow(:, 4);
    raw_flow_s(raw_flow_s < 0.1) = 0.1;
    mov_flow_s = raw_flow_s;

    if mov_mean == 0

        for n = 1:mov_mean
            mov_flow_s = smoothdata(mov_flow_s, 'movmean', 3);
        end

    end

    %只提取十五年20200617
    % temp=min(5476,length(mov_flow_s));
    % mov_flow_s=mov_flow_s(1:temp,1);
    div_mov_flow = mov_flow_s(2:end) ./ mov_flow_s(1:end - 1);
    index1 = find(div_mov_flow > 0.98);
    index2 = find(diff(index1) >= min_dura);

    for ii = 1:size(index2, 1)
        recession = mov_flow_s(index1(index2(ii)) + 1:index1(index2(ii) + 1));
        raw_rece = raw_flow_s(index1(index2(ii)) + 1:index1(index2(ii) + 1));
        r_index = (index1(index2(ii)) + 1:index1(index2(ii) + 1))';
        rece_series{ii, 1} = [r_index, recession, raw_rece];
        dura(ii, 1) = length(r_index);
    end

end
