function [ peaks_datenum ] = selectpeaks( data ,threshold,interval )
%selectpeaks ����POT����ѡ���岢���ж������ж�
%   ��Ҫ�ǶԸ����İٷ�λ��ֵ������ж������ļ��ʱ����к��ɸѡ�뺯��POTpeak������ͬ
%% ���ĳ���ɸѡ
data1=data(:,4);
%�������㻬��ƽ��
% for mm=1:1
% data1=smoothdata(data1,'movmean',3);
% end
flow_s=data1;
flow_datevec=data(:,1:3);
flow_datenum=datenum(flow_datevec);
flow_s_sort=sort(flow_s);
%20200616gai���ٷֱ����뻻Ϊ��ֱֵ������
f_threshold=threshold;
%f_threshold=flow_s_sort(floor(length(flow_s)*Quantile/100));%��50%��λѡ����ֵ
% index_exce=find(flow_s>f_threshold);
%��ѡ���Ӧ��findpeaks
[peaks,peaksdatenum]=findpeaks(flow_s,flow_datenum,'MinPeakHeight',f_threshold);

% f_exce=flow_s(index_exce);
% % ��ѡ���
% num_1=0;
% for num=1:length(f_exce)
%     if index_exce(num)==1||index_exce(num)==length(flow_s)
%         continue
%     end
%     e_1=flow_s(index_exce(num)-1);
%     e_2=flow_s(index_exce(num));
%     e_3=flow_s(index_exce(num)+1);
%     if e_1<e_2&&e_3<=e_2
%         if e_2==e_3&&e_2<=min(flow_s(index_exce(num):index_exce(num)+5))
%             continue
%         end
%         num_1=num_1+1;
%         peaks(num_1,1)=e_2;
%         peaksdate(num_1,1:3)=flow_datevec(index_exce(num),1:3);
%         peaksdatenum(num_1,1:3)=flow_datenum(index_exce(num));
%     end
% end
% peaksdatenum=datenum(peaksdate);

%% ����������ԣ�20191211��
datediff=diff(peaksdatenum);
for ii=1:length(peaks)-1
    index_1=find(flow_datenum==peaksdatenum(ii));
    index_2=find(flow_datenum==peaksdatenum(ii+1));
    span=index_1:index_2;
    midqmin=min(flow_s(span));
    Qmin66=3/4*min(peaks(ii),peaks(ii+1));
    if datediff(ii)>=interval&&midqmin<Qmin66
        event_count(ii,1)=ii;%���������ĺ����
    end
end
%% ��Գ���˫��ĺ�ˮ�����ж�
for kk=1:10
    c_num=0;
indextemp1=find(event_count==0);%
indextemp2=find(diff(indextemp1)>1);
indextemp2=[0;indextemp2;length(indextemp1)];
for nn=1:size(indextemp2,1)-1
    if indextemp1(indextemp2(nn+1))-indextemp1(indextemp2(nn)+1)+2==3
       
            datespan=(peaksdatenum(indextemp1(indextemp2(nn)+1)):peaksdatenum(indextemp1(indextemp2(nn+1))+1))';
            eventflow=flow_s(flow_datenum>=datespan(1)&flow_datenum<=datespan(end));
            minvalue=min(eventflow);
            minvalue_date=datespan(eventflow==minvalue);
            if length(datespan)>=interval&&minvalue<3/4*min(peaks(indextemp1(indextemp2(nn)+1)),peaks(indextemp1(indextemp2(nn+1))+1))
                c_num=c_num+1;
                if minvalue_date<peaksdatenum(indextemp1(indextemp2(nn)+1)+1)
                    event_count(indextemp1(indextemp2(nn)+1))=indextemp1(indextemp2(nn)+1);
                else
                    event_count(indextemp1(indextemp2(nn)+1)+1)=indextemp1(indextemp2(nn)+1)+1;
                end
                
            end
    elseif  indextemp1(indextemp2(nn+1))-indextemp1(indextemp2(nn)+1)+2>=4
            datespan=(peaksdatenum(indextemp1(indextemp2(nn)+1)):peaksdatenum(indextemp1(indextemp2(nn+1))+1))';
            eventflow=flow_s(flow_datenum>=datespan(1)&flow_datenum<=datespan(end));
            minvalue=min(eventflow);
            minvalue_date=datespan(eventflow==minvalue);
            index_min=find(peaksdatenum(indextemp1(indextemp2(nn)+1):indextemp1(indextemp2(nn+1))+1)<minvalue_date(1));
            temp=peaks(indextemp1(indextemp2(nn)+1):indextemp1(indextemp2(nn+1))+1);
            %20200108���ж�����
%             temp(temp<quantile(temp,0.5))=[];
            %how to judge
            if minvalue<3/4*max(temp)
                c_num=c_num+1;
            event_count(indextemp1(indextemp2(nn)+1)+index_min(end)-1)=indextemp1(indextemp2(nn)+1)+index_min(end)-1;       
            end
    end
            
end
if c_num==0
    break
end
end
event_count(event_count==0)=[];
event_p_info=zeros(length(event_count)+1,3);
event_p_datenuminfo=zeros(length(event_count)+1,3);

for jj=0:length(event_count)
    if jj==0
        [event_peak,loc]=max(peaks(1:event_count(jj+1)));
        e_p_star=peaks(1);
        e_p_end=peaks(event_count(jj+1));
        e_p_star_dnum=peaksdatenum(1);
        e_p_end_dnum=peaksdatenum(event_count(jj+1));
        e_p_datenum=peaksdatenum(loc);
    elseif jj==length(event_count)
        [event_peak,loc]=max(peaks(event_count(jj)+1:end));
        e_p_star=peaks(event_count(jj)+1);
        e_p_end=peaks(end);
        e_p_star_dnum=peaksdatenum(event_count(jj)+1);
        e_p_end_dnum=peaksdatenum((end));
        e_p_datenum=peaksdatenum(event_count(jj)+loc);
    else
        [event_peak,loc]=max(peaks(event_count(jj)+1:event_count(jj+1)));
        e_p_star=peaks(event_count(jj)+1);
        e_p_end=peaks(event_count(jj+1));
        e_p_star_dnum=peaksdatenum(event_count(jj)+1);
        e_p_end_dnum=peaksdatenum(event_count(jj+1));
        e_p_datenum=peaksdatenum(event_count(jj)+loc);
    end
    
    event_p_info(jj+1,:)=[e_p_star,event_peak,e_p_end];
    event_p_datenuminfo(jj+1,:)=[e_p_star_dnum,e_p_datenum,e_p_end_dnum];
end
peaks_datenum=[event_p_info,event_p_datenuminfo];
% % ��������ԣ�ȥ���������ĺ�壬����һ����С���ʱ��
% for jj=1:10
%     peaks_1=peaks;
%     peaksdatenum_1=peaksdatenum;
%     for ii=1:length(peaks)-1
%         index_1=find(flow_datenum==peaksdatenum(ii));
%         index_2=find(flow_datenum==peaksdatenum(ii+1));
%         span=index_1:index_2;
%         midq=flow_s(span);
%         Qmin50=3/4*min(peaks(ii),peaks(ii+1));
%         if peaksdatenum(ii+1)-peaksdatenum(ii)<=interval||isempty(find(midq<=Qmin50, 1))==1
%             if peaks(ii+1)>peaks(ii)
%                 peaks_1(ii)=0;peaksdatenum_1(ii)=0;
%             else
%                 peaks_1(ii+1)=0;peaksdatenum_1(ii+1)=0;%��С�ڵķ�ֵ��ֵΪ0
%             end
%         end
%     end
%     peaks=peaks_1(peaks_1>0);
%     peaksdatenum=peaksdatenum_1(peaksdatenum_1>0);
%     if length(peaks)==length(peaks_1)
%         break
%     end
% end
% peaksdatevec=datevec(peaksdatenum);
% date_peaks=[peaksdatevec(:,1:3),peaks];

end

