function [newflood_p, newflood_f, baseflow, K] = newfloodevents(data, date_peaks, s_e_date_q, dura)
    %floodfeature 根据退水曲线与基流修正场次洪水特征
    %   data为原始流量过程，date_peaks为洪水事件，newflood为修正后的洪水过程，eventsfeature为修正后的洪水特征值，包含洪峰洪量与洪水持续时间
    %% 计算退水曲线的参数
    % 算法：将所有流量序列的后一时间段的流量除以前一天的流量，取对应值为0.5—0.96间的为退水过程的流量消退系数，求均值得到K
    flow_s = data(:, 4);
    flow_datenum = datenum(data(:, 1:3));
    flow_s_sort = sort(flow_s);
    %baseflow=mean(flow_s_sort(1:floor(length(flow_s)*0.05)));%定义基流
    baseflow = quantile(s_e_date_q(:, 4), 0.1); %定义基流
    baseflow(baseflow <= 1) = 1;
    %% remove old method
    % cg=flow_s(2:end)./flow_s(1:end-1);
    % cg(cg<0.5)=0;cg(cg>0.95)=0;
    % cg=cg(cg>0);
    % cg(isempty(cg))=0.94;
    % K=-1/log(mean(cg));
    %% 应用退水曲线切割径流
    %未去掉基流

    stardatenum = datenum(s_e_date_q(:, 1:3));
    enddatenum = datenum(s_e_date_q(:, 5:7));
    newflood_p = cell(size(date_peaks, 1), 1); %洪水过程
    newflood_f = zeros(size(date_peaks, 1), 3); %洪水特征，三列分别为洪峰，洪量，持续时间

    for ii = 1:size(date_peaks, 1)
        e_process = flow_s(find(flow_datenum == stardatenum(ii)):find(flow_datenum == enddatenum(ii))); %场次洪水径流过程
        step = dura(ii) / 2;
        step(step < 10) = 10; step(step > 30) = 30;
        T = (0:ceil(step))';

        if e_process(1) <= baseflow && e_process(end) <= baseflow
            qguocheng_3 = e_process;
        elseif e_process(1) <= baseflow && e_process(end) > baseflow
            receflow = e_process(end) * exp(-T / K); %本场洪水退水曲线
            qguocheng_3 = [e_process; receflow(2:end)]; %原始洪水过程与退水曲线相接
        elseif e_process(1) > baseflow && e_process(end) <= baseflow
            lastreceflow = e_process(1) * exp(-T / K); %上一场退水曲线
            lastreceflow_1 = lastreceflow - baseflow; %退水曲线与基流的差值
            lastreceflow_1 = lastreceflow_1(lastreceflow_1 > 0); %提取大于0的部分

            if length(lastreceflow_1) > length(e_process)
                lastreceflow_1 = lastreceflow_1(1:length(e_process));
            end

            qguocheng_2 = e_process(1:length(lastreceflow_1)) - lastreceflow_1;
            qguocheng_3 = [qguocheng_2; e_process(length(lastreceflow_1) + 1:end)]; %修正后的洪水过程

        else
            %     e_process(e_process<=baseflow)=baseflow;%去掉
            lastreceflow = e_process(1) * exp(-T / K); %上一场退水曲线
            lastreceflow_1 = lastreceflow - baseflow; %退水曲线与基流的差值
            lastreceflow_1 = lastreceflow_1(lastreceflow_1 > 0); %提取大于0的部分
            receflow = e_process(end) * exp(-T / K); %本场洪水退水曲线
            qguocheng_1 = [e_process; receflow(2:end)]; %原始洪水过程与退水曲线相接
            qguocheng_1 = qguocheng_1(qguocheng_1 >= baseflow); %剔除尾部小于基流的部分

            if length(lastreceflow_1) > length(qguocheng_1)
                lastreceflow_1 = lastreceflow_1(1:length(qguocheng_1));
            end

            qguocheng_2 = qguocheng_1(1:length(lastreceflow_1)) - lastreceflow_1;
            qguocheng_3 = [qguocheng_2; qguocheng_1(length(lastreceflow_1) + 1:end)]; %修正后的洪水过程
        end

        realstartime = stardatenum(ii, 1);
        realendtime = stardatenum(ii, 1) + length(qguocheng_3) - 1;
        span = (realstartime:realendtime)';
        spandate = datevec(span);
        floodprocess = [spandate(:, 1:3), qguocheng_3]; %洪水过程未减基流
        newflood_f(ii, 1) = max(qguocheng_3 - baseflow);
        newflood_f(ii, 2) = sum(qguocheng_3 - baseflow) * 24 * 3600/10^4;
        newflood_p{ii, 1} = floodprocess;
        newflood_f(ii, 3) = realendtime - realstartime + 1; %洪水特征值均减了基流

    end
