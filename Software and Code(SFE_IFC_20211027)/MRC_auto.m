function [posi,sort_num,MRCdata,COR_MRCdata] = 	MRC_auto(rece_series,dura)
%MRC_auto 平移退水过程，求出模板退水曲线
%   此处显示详细说明
recessions=zeros(max(dura),length(rece_series));
for ii=1:length(rece_series)
    recess=rece_series{ii,1}(:,2);
    rece_min(ii,2)=recess(end);
    rece_min(ii,1)=ii;
    recessions(1:length(recess),ii)=recess;   
end
sort_rece_min=sortrows(rece_min,2);
sort_num=sort_rece_min(:,1);


%jj=1时情景
posi_end(1,1)=max(dura);
posi_star(1,1)=posi_end(1,1)-dura(sort_num(1))+1;


for jj=2:length(rece_series)
    s_recess1=rece_series{sort_num(jj-1),1}(:,2);
    s_r_f1=flipud(s_recess1);
    s_recess2=rece_series{sort_num(jj),1}(:,2);
    s_r_f2=flipud(s_recess2);
    s_temp1=s_r_f2(1)-s_r_f1;
     s_temp2=s_r_f2(2)-s_r_f1;
    m_s=sqrt(s_temp1(1:end-1).^2+s_temp2(2:end).^2);
    [~,index]=min(m_s);
    posi_end(jj,1)=posi_end(jj-1,1)-index+1;
    posi_star(jj,1)=posi_end(jj,1)-dura(sort_num(jj))+1;
end
posi=[posi_star,posi_end];
posi=posi-min(posi(:))+1;
dpeindex=-diff(posi_end);
posi_MRC=[];
data_MRC=[];
for mm=1:length(dpeindex)
    posi_mrc=(posi(mm,2)-dpeindex(mm):posi(mm,2))';
    data_mrc=rece_series{sort_num(mm),1}(end-dpeindex(mm):end,2);
    posi_MRC=[posi_MRC;posi_mrc];
    data_MRC=[data_MRC;data_mrc];
end
posi_MRC=[posi_MRC;(posi(length(dpeindex)+1,1):posi(length(dpeindex)+1,2))'];
data_MRC=[data_MRC;rece_series{sort_num(length(dpeindex)+1),1}(:,2)];
MRCdata=[posi_MRC,data_MRC];
for nn=min(posi_MRC):max(posi_MRC)
    cor_MRCdata(nn,1)=mean(data_MRC(posi_MRC==nn));
end
COR_MRCdata=[(min(posi_MRC):max(posi_MRC))',cor_MRCdata(min(posi_MRC):max(posi_MRC))];
end


