function [oldflood_p, oldflood_f, T] = oldfloodevents(data, starendtime, dura)
    %flood_events_plot 绘制场次洪水过程线
    %   data前三列为径流序列的年月日，第四列为流量，starendtime前三列为开始时间，后三列为结束时间
    oldflood_p = cell(size(starendtime, 1), 1); %洪水过程
    oldflood_f = zeros(size(starendtime, 1), 3); %洪水特征，三列分别为洪峰，洪量，持续时间

    for ii = 1:size(starendtime, 1)
        %     temp=ceil(dura(ii)/3);
        flow_datenum = datenum(data(:, 1:3));
        star_datenum = datenum(starendtime(ii, 1:3));
        end_datenum = datenum(starendtime(ii, 4:6));
        star_index = find(flow_datenum == star_datenum);
        end_index = find(flow_datenum == end_datenum);
        %     plot(flow_datenum(star_index-temp:end_index+temp),data(star_index-temp:end_index+temp,4),'LineWidth',2)
        flooddura_t = datevec(flow_datenum(star_index:end_index));
        oldflood_p{ii, 1} = [flooddura_t(:, 1:3), data(star_index:end_index, 4)];
        oldflood_f(ii, 1) = max(data(star_index:end_index, 4));
        oldflood_f(ii, 2) = sum(data(star_index:end_index, 4)) * 24 * 3600/10^4;
        oldflood_f(ii, 3) = dura(ii, 1);

        %% 设置坐标轴与参考线
        %     ylim=get(gca,'Ylim'); % 获取当前图形的纵轴的范围
        %     hold on
        %     plot([flow_datenum(star_index),flow_datenum(star_index)],ylim,'m--','LineWidth',1.5); hold on
        %     plot([flow_datenum(end_index),flow_datenum(end_index)],ylim,'m--','LineWidth',1.5)
        %     dateaxis('x',6)
        %     title(datestr(star_datenum,26))
        %     print(gcf,'-r600','-dpng',[num2str([ii,starendtime(ii,1:3)]),'.png']);
        %     close
        % end
    end

    Number = (1:size(oldflood_f))';
    Stardate = datetime(starendtime(:, 1:3));
    Enddate = datetime(starendtime(:, 4:6));
    Peakvalue = oldflood_f(:, 1);
    Volume = oldflood_f(:, 2);
    Duration = oldflood_f(:, 3);
    T = table(Number, Stardate, Enddate, Peakvalue, Volume, Duration);
