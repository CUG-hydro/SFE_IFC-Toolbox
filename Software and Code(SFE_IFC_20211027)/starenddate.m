function [ s_e_date_q,dura,f_low] = starenddate( data,events,interval)
%starenddate 初步确定洪水事件的起始终止时间
%   data为有效径流，前三列为时间，后一列为流量，events为事件，同四列,interval为事件独立性判断的间隔时间
%   算法为，寻找洪峰前后n天低于洪峰一半的最小值点，n按照EPA的洪水独立性来判断
% 20191211修改

data1=data(:,4);
data2=data(:,4);
% for mm=1:1
% data1=smoothdata(data1,'movmean',3);
% end

qdata=data1;
qdate=data(:,1:3);
qdatenum=datenum(qdate);
% peak=events(:,4);
%  bpeakdatenum=datenum(events(:,1:3));
flow_s_sort=sort(qdata);
f_low=flow_s_sort(floor(length(qdata)*0.1));
if f_low==0
    f_low=0.1;
end
qdata(qdata<0.1)=0.1;
%% 求事件的起始终止时间
for ii=1:size(events,1)
    interval1=interval;
    interval2=interval;
    if ii>1&&events(ii,6)-events(ii-1,4)<interval1
        interval1=events(ii,6)-events(ii-1,4);
    end 
    sindex=find(qdatenum==events(ii,4));%找的时间在径流序列的位置
    ss=sindex-interval1-1;
    for day=1:interval1 %从洪峰往前遍历
        
        if ss<1
            ss=1;
            [guaidianq(ii,1),lca]=min(qdata(ss:sindex));
            guaidiand(ii,1)=qdatenum(lca);
            break
            % 倘若前15天不存在假设的点，则取前15天流量的最小值
        end
        %             处理事件前的缺测情况
        if qdatenum(sindex-day)-qdatenum(sindex-day-1)>1
            guaidianq(ii,1)=qdata(sindex-day);
            guaidiand(ii,1)=qdatenum(sindex-day);
            break
        end
        
        
        if (qdata(sindex-day)-qdata(sindex-day-1)<=0&&qdata(sindex-day)-qdata(sindex-day+1)<0)&&(qdata(sindex-day)<0.5*events(ii,2)||qdata(sindex-day)<=f_low)
%             
            guaidianq(ii,1)=qdata(sindex-day);
            guaidiand(ii,1)=qdatenum(sindex-day);
            break   %找洪峰前15天低于峰值50%的波谷，跳出循环
        else
            if day==interval1
            guaidianq(ii,1)=min(qdata(sindex-interval1:sindex));
            guaidiand(ii,1)=qdatenum(sindex+max(find(qdata(sindex-interval1:sindex)==min(qdata(sindex-interval1:sindex))))-interval1-1);    % 倘若前15天不存在假设的点，则取前15天流量的最小值
            end
        end
        
    end
    
    eindex=find(qdatenum==events(ii,6));
    ee=eindex+interval2;
    for day=1:interval2
        %处理最后一个事件
        if ee>=length(qdata)
            ee=length(qdata);
            [guaidianq(ii,2),lca]=min(qdata(eindex:ee));
            guaidiand(ii,2)=qdatenum(lca+eindex-1);
            break
        end
        
        %处理存在缺测的情况
        if qdatenum(eindex+day)-qdatenum(eindex+day+1)<-1
            guaidianq(ii,2)=qdata(eindex+day);
            guaidiand(ii,2)=qdatenum(eindex+day);
            break
        end
        
        if (qdata(eindex+day)-qdata(eindex+day-1)<0&&qdata(eindex+day)-qdata(eindex+day+1)<=0)&&(qdata(eindex+day)<0.5*events(ii,2)||qdata(eindex+day)<=f_low)
            
            guaidianq(ii,2)=qdata(eindex+day);
            guaidiand(ii,2)=qdatenum(eindex+day);
            break
        else
            if day==interval2
            guaidianq(ii,2)=min(qdata(eindex:eindex+interval2));
            guaidiand(ii,2)=qdatenum(eindex+min(find(qdata(eindex:eindex+interval2)==min(qdata(eindex:eindex+interval2))))-1);
            end
        end
        
    end
    
end
%% 三点滑动平均后微调
% for jj=1:size(events,1)
%     sindex_c=find(qdatenum==guaidiand(jj,1));
%     if sindex_c>1
%     [value,loc_c]=min(data2(sindex_c-1:sindex_c+1));
%     guaidiand(jj,1)=guaidiand(jj,1)-2+loc_c;
%     guaidianq(jj,1)=value;
%     end
%     
%     eindex_c=find(qdatenum==guaidiand(jj,2));
%     if eindex_c<length(qdatenum)
%     [value,loc_c]=min(data2(eindex_c-1:eindex_c+1));
%     guaidiand(jj,2)=guaidiand(jj,2)-2+loc_c;
%     guaidianq(jj,2)=value;
%     end
%     
% end
%% 将事件开始结束的时间输出
duration=guaidiand(:,2)-guaidiand(:,1)+1;
star_peak=events(:,4)-guaidiand(:,1)+1;
faileventsimdex=find(star_peak==interval+1);
dura=[duration,star_peak];
stardate=datevec(guaidiand(:,1));
enddate=datevec(guaidiand(:,2));
s_e_date_q=[stardate(:,1:3),guaidianq(:,1),enddate(:,1:3),guaidianq(:,2)];

end

