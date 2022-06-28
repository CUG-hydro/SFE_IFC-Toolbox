function [newflood_p,newflood_f,T2] = FloodCharacteristics(date_flow,s_e_date_q,COR_MRCdata)
%FloodCharacteristics 识别洪水特征
%   此处显示详细说明
flow_s=date_flow(:,4);
flow_datenum=datenum(date_flow(:,1:3));
stardatenum=datenum(s_e_date_q(:,1:3));
enddatenum=datenum(s_e_date_q(:,5:7));
newflood_p=cell(size(s_e_date_q,1),1);%洪水过程
newflood_f=zeros(size(s_e_date_q,1),3);%洪水特征，三列分别为洪峰，洪量，持续时间

% baseflow=quantile(s_e_date_q(:,4),0.1);%定义基流
baseflow=Deepbaseflow(date_flow);
 baseflow(baseflow<=1)=1;% 20200409改
% a=para(1);
% b=para(2);
% in_value=COR_MRCdata(1,2);
fit_t=COR_MRCdata(:,1);
fit_t1=round(fit_t*100);
fit_q=COR_MRCdata(:,2);



for ii=1:size(s_e_date_q,1)
    e_process=flow_s(find(flow_datenum==stardatenum(ii)):find(flow_datenum==enddatenum(ii)));%场次洪水径流过程
%     step=dura(ii);
%     step(step<10)=10;step(step>60)=60;
    T=(0:1:70)';
        t1=fit_t(abs(fit_q-e_process(1))==min(abs(fit_q-e_process(1))));
        t2=fit_t(abs(fit_q-e_process(end))==min(abs(fit_q-e_process(end))));
        Tt1=T+t1;
        Tt2=T+t2;
        Tt1(Tt1>fit_t(end))=[];
        Tt2(Tt2>fit_t(end))=[];
        Tt11=round(Tt1*100);
        Tt22=round(Tt2*100);
        [~,l_locb]=ismember(Tt11,fit_t1);
        [~,r_locb]=ismember(Tt22,fit_t1);
        lastreceflow=fit_q(l_locb);
        receflow=fit_q(r_locb);

lastreceflow(abs(diff(lastreceflow))<baseflow*0.01)=[];
receflow(abs(diff(receflow))<baseflow*0.01)=[];
    
    
    if e_process(1)<=baseflow&&e_process(end)<=baseflow
        qguocheng_3=e_process;
    elseif e_process(1)<=baseflow&&e_process(end)>baseflow
        receflow(receflow<baseflow)=[];
        qguocheng_3=[e_process;receflow(2:end)];%原始洪水过程与退水曲线相接
    elseif e_process(1)>baseflow&&e_process(end)<=baseflow
        lastreceflow_1=lastreceflow-baseflow;%退水曲线与基流的差值
        lastreceflow_1=lastreceflow_1(lastreceflow_1>0);%提取大于0的部分
        lastreceflow_2=lastreceflow_1+baseflow;
        if length(lastreceflow_1)>length(e_process)
%             lastreceflow_1=lastreceflow_1(1:length(e_process));
            lastreceflow_2=lastreceflow_2(1:length(e_process));
            temp_judge=e_process-lastreceflow_2;
            temp_judge(1)=0.1;
            temp_judge(temp_judge<0)=[];
            qguocheng_3=temp_judge+lastreceflow_2(length(temp_judge));
        else
        qguocheng_2=e_process(1:length(lastreceflow_1))-lastreceflow_1;
        qguocheng_3=[qguocheng_2;e_process(length(lastreceflow_1)+1:end)];%修正后的洪水过程
        end
    else
        %     e_process(e_process<=baseflow)=baseflow;%去掉
        lastreceflow_1=lastreceflow-baseflow;%退水曲线与基流的差值
        lastreceflow_1=lastreceflow_1(lastreceflow_1>0);%提取大于0的部分
        lastreceflow_2=lastreceflow_1+baseflow;
        qguocheng_1=[e_process;receflow(2:end)];%原始洪水过程与退水曲线相接
        qguocheng_1=qguocheng_1(qguocheng_1>=baseflow);%剔除尾部小于基流的部分
        if length(lastreceflow_1)>length(qguocheng_1)
%             lastreceflow_1=lastreceflow_1(1:length(qguocheng_1));
            lastreceflow_2=lastreceflow_2(1:length(qguocheng_1));
            temp_judge=qguocheng_1-lastreceflow_2;
             temp_judge(1)=0.1;
            temp_judge(temp_judge<0)=[];
            qguocheng_3=temp_judge+lastreceflow_2(length(temp_judge));
        else
        qguocheng_2=qguocheng_1(1:length(lastreceflow_1))-lastreceflow_1;
        qguocheng_3=[qguocheng_2;qguocheng_1(length(lastreceflow_1)+1:end)];%修正后的洪水过程
        end
    end
    realstartime=stardatenum(ii,1);
    realendtime=stardatenum(ii,1)+length(qguocheng_3)-1;
    span=(realstartime:realendtime)';
    spandate=datevec(span);
    realenddate1(ii,1)=realendtime;
    
    floodprocess=[spandate(:,1:3),qguocheng_3];%洪水过程未减基流
    newflood_f(ii,1)=max(qguocheng_3-min([baseflow,min(qguocheng_3),0]));
    newflood_f(ii,2)=sum(qguocheng_3-min([baseflow,min(qguocheng_3),0]))*24*3600/10^4;
    newflood_p{ii,1}=floodprocess;
    newflood_f(ii,3)=realendtime-realstartime+1;%洪水特征值均减了基流
    
    
end
Number=(1:size(newflood_f))';
Stardate=datetime(s_e_date_q(:,1:3));
Enddate=datetime(datevec(realenddate1),'Format','uuuu-MM-dd');
Peakvalue=newflood_f(:,1);
Volume=newflood_f(:,2);
Duration=newflood_f(:,3);
T2=table(Number,Stardate,Enddate,Peakvalue,Volume,Duration);
end

