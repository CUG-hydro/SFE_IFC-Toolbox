function [ newflood_p,newflood_f ,baseflow,K] = newfloodevents( data,date_peaks ,s_e_date_q,dura)
%floodfeature ������ˮ����������������κ�ˮ����
%   dataΪԭʼ�������̣�date_peaksΪ��ˮ�¼���newfloodΪ������ĺ�ˮ���̣�eventsfeatureΪ������ĺ�ˮ����ֵ���������������ˮ����ʱ��
%% ������ˮ���ߵĲ���
% �㷨���������������еĺ�һʱ��ε���������ǰһ���������ȡ��ӦֵΪ0.5��0.96���Ϊ��ˮ���̵���������ϵ�������ֵ�õ�K
flow_s=data(:,4);
flow_datenum=datenum(data(:,1:3));
flow_s_sort=sort(flow_s);
%baseflow=mean(flow_s_sort(1:floor(length(flow_s)*0.05)));%�������
baseflow=quantile(s_e_date_q(:,4),0.1);%�������
baseflow(baseflow<=1)=1;
%% remove old method
% cg=flow_s(2:end)./flow_s(1:end-1);
% cg(cg<0.5)=0;cg(cg>0.95)=0;
% cg=cg(cg>0);
% cg(isempty(cg))=0.94;
% K=-1/log(mean(cg));
%% Ӧ����ˮ�����и��
%δȥ������

stardatenum=datenum(s_e_date_q(:,1:3));
enddatenum=datenum(s_e_date_q(:,5:7));
newflood_p=cell(size(date_peaks,1),1);%��ˮ����
newflood_f=zeros(size(date_peaks,1),3);%��ˮ���������зֱ�Ϊ��壬����������ʱ��
for ii=1:size(date_peaks,1)
    e_process=flow_s(find(flow_datenum==stardatenum(ii)):find(flow_datenum==enddatenum(ii)));%���κ�ˮ��������
    step=dura(ii)/2;
    step(step<10)=10;step(step>30)=30;
    T=(0:ceil(step))';
    if e_process(1)<=baseflow&&e_process(end)<=baseflow
        qguocheng_3=e_process;
    elseif e_process(1)<=baseflow&&e_process(end)>baseflow
        receflow=e_process(end)*exp(-T/K);%������ˮ��ˮ����
        qguocheng_3=[e_process;receflow(2:end)];%ԭʼ��ˮ��������ˮ�������
    elseif e_process(1)>baseflow&&e_process(end)<=baseflow
        lastreceflow=e_process(1)*exp(-T/K);%��һ����ˮ����
        lastreceflow_1=lastreceflow-baseflow;%��ˮ����������Ĳ�ֵ
        lastreceflow_1=lastreceflow_1(lastreceflow_1>0);%��ȡ����0�Ĳ���
        if length(lastreceflow_1)>length(e_process)
            lastreceflow_1=lastreceflow_1(1:length(e_process));
        end
        qguocheng_2=e_process(1:length(lastreceflow_1))-lastreceflow_1;
        qguocheng_3=[qguocheng_2;e_process(length(lastreceflow_1)+1:end)];%������ĺ�ˮ����
        
    else
        %     e_process(e_process<=baseflow)=baseflow;%ȥ��
        lastreceflow=e_process(1)*exp(-T/K);%��һ����ˮ����
        lastreceflow_1=lastreceflow-baseflow;%��ˮ����������Ĳ�ֵ
        lastreceflow_1=lastreceflow_1(lastreceflow_1>0);%��ȡ����0�Ĳ���
        receflow=e_process(end)*exp(-T/K);%������ˮ��ˮ����
        qguocheng_1=[e_process;receflow(2:end)];%ԭʼ��ˮ��������ˮ�������
        qguocheng_1=qguocheng_1(qguocheng_1>=baseflow);%�޳�β��С�ڻ����Ĳ���
        if length(lastreceflow_1)>length(qguocheng_1)
            lastreceflow_1=lastreceflow_1(1:length(qguocheng_1));
        end
        qguocheng_2=qguocheng_1(1:length(lastreceflow_1))-lastreceflow_1;
        qguocheng_3=[qguocheng_2;qguocheng_1(length(lastreceflow_1)+1:end)];%������ĺ�ˮ����
    end
    realstartime=stardatenum(ii,1);
    realendtime=stardatenum(ii,1)+length(qguocheng_3)-1;
    span=(realstartime:realendtime)';
    spandate=datevec(span);
    floodprocess=[spandate(:,1:3),qguocheng_3];%��ˮ����δ������
    newflood_f(ii,1)=max(qguocheng_3-baseflow);
    newflood_f(ii,2)=sum(qguocheng_3-baseflow)*24*3600/10^4;
    newflood_p{ii,1}=floodprocess;
    newflood_f(ii,3)=realendtime-realstartime+1;%��ˮ����ֵ�����˻���
    
end

