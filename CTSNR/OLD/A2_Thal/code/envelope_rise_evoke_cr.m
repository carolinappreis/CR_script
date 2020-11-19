% cd('/Users/Carolina/Documents/MATLAB/hayriye code for tc analysis')
% open('envelope_rise_evoke')
clear all
% cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/mat')
load ('lesion_param_2000_BZ_new.mat','evoked_all','power_all','beta_coh_all');
for i =1:size(beta_coh_all,1)
beta_coh_all{i}=beta_coh_all{i}';
end
beta_coh_all=beta_coh_all(~cellfun('isempty',beta_coh_all));
P=cell2mat(beta_coh_all);
evoked_all=evoked_all(~cellfun('isempty',evoked_all));
evoked_all=cell2mat(evoked_all);
power_all=power_all(~cellfun('isempty',power_all));
power_all=cell2mat(power_all);
output5=evoked_all;
Pxx=power_all;


for i=1:size(output5,1)
    [b,a]=butter(2,[15/(0.5*1000) 35/(0.5*1000)]);
    d(i,:)=abs(hilbert(filtfilt(b,a,output5(i,:))));
    [m,e(i)]=max(d(i,1800:2200));
end
for i=1:size(output5,1)
    dd(i,:)=(d(i,:)-mean(d(i,:)))./std(d(i,:));
    dd_max(i)=dd(i,e(i)+1799);
end
%P=sum(Pxx(:,16:36),2)./sum(Pxx(:,1:501),2);
sub=find(P>0.1 & dd_max'>1.96);
for k=1:length(sub)
    dum_var2=diff(smooth(d(sub(k),:),101));
    dum_var4=dum_var2(1800:e(sub(k))+1799);
    x0=find(dum_var4>0);
    if ~isempty(find(diff(x0)>1))
        XXXstr(k,1)=x0(max(find(diff(x0)>1))+1)-199;
    elseif isempty(find(diff(x0)>1)) && ~isempty(x0)
        XXXstr(k,1)=x0(1)-199;
        onset=1800;
        while x0(1)==1
            dum_var4=dum_var2(onset-50:e(sub(k))+1799);
            x0=find(dum_var4>0);
            onset=onset-50;
            if ~isempty(find(diff(x0)>1))
                XXXstr(k,1)=x0(max(find(diff(x0)>1))+1)-(2000-onset+1);
            elseif isempty(find(diff(x0)>1)) && ~isempty(x0)
                XXXstr(k,1)=x0(1)-(2000-onset+1);
            end
            
        end
    elseif isempty(x0)
        clear dum_var4
        dum_var4=dum_var2(1800:2200);
        x00=find(dum_var4>0);
        if ~isempty(x00)
            XXXstr(k,1)=x00(1)-199;
        elseif isempty(x00)
            XXXstr(k,1)= NaN;
        end
    end
    clear dumvar2 dumvar4 x0
end
% mean(XXXstr)