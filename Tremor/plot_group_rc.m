clear all
new=[1:5 7:9]; %%without PD patient (i.e., pt number 6)
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\A_group.mat')
% load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/A_group.mat')
a.s=S(new,:); clearvars -except a
% load ('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/F_group1.mat')
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\F_group1.mat')
f.s=S; clearvars -except a f

ref1=1;%%%% amp(=0) vs. frequency
iii=1; %%%%% amp (=0) vs. supressive effect
metric=1; %%%%%plotting 0 amp; ~=0 freq


if ref1==0
    ref=a.s; %%% max amplitude change vs. max frequecy change
else
    ref=f.s;
end

idmi=[];
idma=[];
for i =1:size(ref,1)
    if (numel(find(ref(i,:)<0)))~=0 % all -except all amplifying subjects
        idmi=[idmi i];
    end
    if (numel(find(ref(i,:)>0)))~=0 % all -except all supressive subjects
        idma=[idma i];
    end
end

if iii==0;
    sub=ref(idma,:);  % iii=0 amplifying effect;
else
    sub=ref(idmi,:);  % iii~=0 supressive effect;
end



for i=1:size(sub,1);
    if iii==0
        phase_peak(1,i)=find(sub(i,:)==max(sub(i,:)));
    else
        phase_peak(1,i)=find(sub(i,:)==min(sub(i,:)));
    end
    
    if phase_peak(i)==1;
        a_s_al(i,:)=a.s(i,:);
        f_s_al(i,:)=f.s(i,:);
        
        
    else
        a_s_al(i,:)=[a.s(i,phase_peak(i):end) a.s(i,1:phase_peak(i)-1)];
        f_s_al(i,:)=[f.s(i,phase_peak(i):end) f.s(i,1:phase_peak(i)-1)];
        
        check=[phase_peak(i):size(a_s_al,2) 1:phase_peak(i)-1];
    end
    
end


for i=1:size(a_s_al,1)
    new_as_al(i,:)=zscore(a_s_al(i,:));
    new_fs_al(i,:)=zscore(f_s_al(i,:));
end



close all

if metric==0;
z_metric=[new_as_al(:,8:12) new_as_al(:,1:7)];
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
cl=blushred;
cl1=squash;

else
z_metric=[new_fs_al(:,8:12) new_fs_al(:,1:7)];
% load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone');
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');
cl=aegean;
cl1=stone;
end

sm=[z_metric z_metric];

for ii=1:size(sm,1)
    dum=[];
    
    for i=1:size(sm,2)
        if i+2<size(sm,2)
            dum=[dum sum(sm(ii,(i:i+2)))./length(sm(ii,(i:i+2)))];
        end
    end
    sm_s(ii,:)=dum;
end

smo_s=sm_s(:,1:12);

figure(1)
plot(smo_s','Color',cl1)
xlim([1 12])
xticks([ 3 6 9 12])
xticklabels({'-90','0','90','180'})
box('off')
hold on
plot(mean(smo_s),'k','LineWidth',4,'Color', cl)
yline(0,'k', 'LineWidth',1,'LineStyle','--')


figure(2)
bar(mean(z_metric),'LineWidth',1,'FaceColor',cl,'EdgeColor',cl)
xlim([0 13])
xticks([ 3 6 9 12])
xticklabels({'-90','0','90','180'})
box('off')

if metric==0
    
    figure(1)
    ylabel('Change in amplitude')
    figure(2)
    ylabel('Change in amplitude')
    
else
    
    figure(1)
    ylabel('Change in frequency')
    figure(2)
    ylabel('Change in frequecy')
end
    


if ref1==0 && iii==0
    
    figure(1)
    xlabel('Aligment to phase with max tremor amplification')
    figure(2)
    xlabel('Aligment to phase with max tremor amplification')
    
elseif ref1==0 && iii==1
    
    figure(1)
    xlabel('Aligment to phase with max tremor supression')
    figure(2)
    xlabel('Aligment to phase with max tremor supression')
    
elseif ref1==1 && iii==0
    
    figure(1)
    xlabel('Aligment to phase with max tremor acceleration')
    figure(2)
    xlabel('Aligment to phase with max tremor acceleration')
    
else ref1==1 && iii==1
    
    figure(1)
    xlabel('Aligment to phase with max tremor deceleration')
    figure(2)
    xlabel('Aligment to phase with max tremor deceleration')
end
    
    
