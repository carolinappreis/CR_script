clear all

dist=2;

if dist==1
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\arc_ARC.mat','ttall')
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\frc_FRC.mat','ttall')
Sa=ttall; 
Sf=ttall;
elseif dist==2 %% arc for a<median
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\arc_mediansplit.mat','arc1')
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\frc_mediansplit.mat','frc1')
Sa=arc1; 
Sf=frc1;
elseif dist==3 %% arc for a>median
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\arc_mediansplit.mat','arc2')
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\frc_mediansplit.mat','frc2')
Sa=arc2; 
Sf=frc2;
end


sm=[Sa Sa Sa];
for ii=1:size(sm,1)
    for i=size(Sa,2)+1:size(Sa,2)*2
        a.s(ii,i-12)=sum(sm(ii,(i-1:i+1)))./length(sm(ii,(i-1:i+1)));
    end
end

sm=[Sf Sf Sf];
for ii=1:size(sm,1)
    for i=size(Sf,2)+1:size(Sf,2)*2
        f.s(ii,i-12)=sum(sm(ii,(i-1:i+1)))./length(sm(ii,(i-1:i+1)));
    end
end
clearvars -except a f



ref1=0;%%%% amp(=0) vs. frequency
iii=0; %%%%% amp (=0) vs. supressive effect
metric=0; %%%%%plotting 0 amp; ~=0 freq

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
    idx=idma;  % iii=0 amplifying effect;
else
    idx=idmi;  % iii~=0 supressive effect;
end

sub=ref(idx,:);

for i=1:length(idx);
    if iii==0
        phase_peak(1,i)=find(sub(i,:)==max(sub(i,:)));
    else
        phase_peak(1,i)=find(sub(i,:)==min(sub(i,:)));
    end
    
    if phase_peak(i)==1;
        a_s_al(i,:)=a.s(idx(i),:);
        f_s_al(i,:)=f.s(idx(i),:);
        
        
    else
        a_s_al(i,:)=[a.s(idx(i),phase_peak(i):end) a.s(idx(i),1:phase_peak(i)-1)];
        f_s_al(i,:)=[f.s(idx(i),phase_peak(i):end) f.s(idx(i),1:phase_peak(i)-1)];
        
        check=[phase_peak(i):size(a_s_al,2) 1:phase_peak(i)-1];
    end
    
end


% close all


aa=[a_s_al(:,8:12) a_s_al(:,1:7)];
%     load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','squash')
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','blushred','squash');
acl=blushred;
acl1=squash;


ff=[f_s_al(:,8:12) f_s_al(:,1:7)];
%     load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','aegean','stone');
load('C:\Users\creis\Documents\GitHub\CR_script\colour_pal.mat','aegean','stone');
fcl=aegean;
fcl1=stone;



figure(1)
for i=1:size(a_s_al,1)
    subplot(1,size(a_s_al,1),i)
    plot(aa(i,:),'Color',acl1)
    hold on
    plot(ff(i,:),'Color',fcl1)
    xlim([1 12])
    xticks([ 3 6 9 12])
    xticklabels({'-90','0','90','180'})
    box('off')
    
end

figure(2)
for i=1:size(a_s_al,1)
    subplot(1,size(a_s_al,1),i)
    y2=plot(aa(i,:),ff(i,:),'k+');
    y3=lsline;
    set(y3,'LineWidth',2,'Color','red')
    box('off')
    c2=corrcoef(aa(i,:)',ff(i,:)')
    legend(y3,[num2str(c2(1,2))],'box','off')
end

figure(3)
subplot(1,2,1)
plot(mean(aa),'k','LineWidth',3,'Color', acl)
hold on
plot(mean(ff),'k','LineWidth',3,'Color', fcl)
yline(0,'k', 'LineWidth',1,'LineStyle','--')
box('off')
subplot(1,2,2)
y2=plot(mean(aa),mean(ff),'k+');
y3=lsline;
set(y3,'LineWidth',2,'Color','red')
box('off')
c2=corrcoef(mean(aa)',mean(ff)')
legend(y3,[num2str(c2(1,2))],'box','off')


if ref1==0 && iii==0
    
    figure(1)
    title('Aligment to phase with max tremor amplification')
    figure(2)
    title('Aligment to phase with max tremor amplification')
    figure(3)
    title('Aligment to phase with max tremor amplification')
    
elseif ref1==0 && iii==1
    
    figure(1)
    title('Aligment to phase with max tremor supression')
    figure(2)
    title('Aligment to phase with max tremor supression')
    figure(3)
    title('Aligment to phase with max tremor supression')
    
    
elseif ref1==1 && iii==0
    
    figure(1)
    title('Aligment to phase with max tremor acceleration')
    figure(2)
    title('Aligment to phase with max tremor acceleration')
    figure(3)
    title('Aligment to phase with max tremor acceleration')
    
    
else ref1==1 && iii==1
    
    figure(1)
    title('Aligment to phase with max tremor deceleration')
    figure(2)
    title('Aligment to phase with max tremor deceleration')
    figure(3)
    title('Aligment to phase with max tremor deceleration')
    
end

