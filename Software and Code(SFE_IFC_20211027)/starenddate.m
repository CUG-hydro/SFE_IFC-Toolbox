function [ s_e_date_q,dura,f_low] = starenddate( data,events,interval)
%starenddate ����ȷ����ˮ�¼�����ʼ��ֹʱ��
%   dataΪ��Ч������ǰ����Ϊʱ�䣬��һ��Ϊ������eventsΪ�¼���ͬ����,intervalΪ�¼��������жϵļ��ʱ��
%   �㷨Ϊ��Ѱ�Һ��ǰ��n����ں��һ�����Сֵ�㣬n����EPA�ĺ�ˮ���������ж�
% 20191211�޸�

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
%% ���¼�����ʼ��ֹʱ��
for ii=1:size(events,1)
    interval1=interval;
    interval2=interval;
    if ii>1&&events(ii,6)-events(ii-1,4)<interval1
        interval1=events(ii,6)-events(ii-1,4);
    end 
    sindex=find(qdatenum==events(ii,4));%�ҵ�ʱ���ھ������е�λ��
    ss=sindex-interval1-1;
    for day=1:interval1 %�Ӻ����ǰ����
        
        if ss<1
            ss=1;
            [guaidianq(ii,1),lca]=min(qdata(ss:sindex));
            guaidiand(ii,1)=qdatenum(lca);
            break
            % ����ǰ15�첻���ڼ���ĵ㣬��ȡǰ15����������Сֵ
        end
        %             �����¼�ǰ��ȱ�����
        if qdatenum(sindex-day)-qdatenum(sindex-day-1)>1
            guaidianq(ii,1)=qdata(sindex-day);
            guaidiand(ii,1)=qdatenum(sindex-day);
            break
        end
        
        
        if (qdata(sindex-day)-qdata(sindex-day-1)<=0&&qdata(sindex-day)-qdata(sindex-day+1)<0)&&(qdata(sindex-day)<0.5*events(ii,2)||qdata(sindex-day)<=f_low)
%             
            guaidianq(ii,1)=qdata(sindex-day);
            guaidiand(ii,1)=qdatenum(sindex-day);
            break   %�Һ��ǰ15����ڷ�ֵ50%�Ĳ��ȣ�����ѭ��
        else
            if day==interval1
            guaidianq(ii,1)=min(qdata(sindex-interval1:sindex));
            guaidiand(ii,1)=qdatenum(sindex+max(find(qdata(sindex-interval1:sindex)==min(qdata(sindex-interval1:sindex))))-interval1-1);    % ����ǰ15�첻���ڼ���ĵ㣬��ȡǰ15����������Сֵ
            end
        end
        
    end
    
    eindex=find(qdatenum==events(ii,6));
    ee=eindex+interval2;
    for day=1:interval2
        %�������һ���¼�
        if ee>=length(qdata)
            ee=length(qdata);
            [guaidianq(ii,2),lca]=min(qdata(eindex:ee));
            guaidiand(ii,2)=qdatenum(lca+eindex-1);
            break
        end
        
        %�������ȱ������
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
%% ���㻬��ƽ����΢��
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
%% ���¼���ʼ������ʱ�����
duration=guaidiand(:,2)-guaidiand(:,1)+1;
star_peak=events(:,4)-guaidiand(:,1)+1;
faileventsimdex=find(star_peak==interval+1);
dura=[duration,star_peak];
stardate=datevec(guaidiand(:,1));
enddate=datevec(guaidiand(:,2));
s_e_date_q=[stardate(:,1:3),guaidianq(:,1),enddate(:,1:3),guaidianq(:,2)];

end

