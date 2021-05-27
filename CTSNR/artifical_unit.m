% Artificial
%% CY BY CY

clear all; close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_SNR.mat')

data=cell2mat(data_region);
ecog=[];
for i=1:size(data_region,1)
    ecog = [ecog ; repmat(Ecog_region(i,:),size(data_region{i,1},1),1)];
end

clearvars -except ecog name

srn=1000;
[b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
tt=0;
gg=0;
all=[];

j=2;
ctx=ecog(j,:);
Ecogfiltered=filtfilt(b,a,ctx);
[maxvalM,maxidxM] = findpeaks(Ecogfiltered);
data(j,:)=zeros(1,length(Ecogfiltered));
data(j,maxidxM)=1;

% plot(Ecogfiltered)
% hold on
% plot(data(j,:).*100)

env=abs(hilbert(Ecogfiltered));
data_ones=find(data(j,:)==1);
hp=wrapToPi(angle(hilbert(Ecogfiltered)));
ang=hp(data_ones);

if (circ_rtest(ang))<0.05
    cy_bursts=cycles_10(env,Ecogfiltered);
    block=cell2mat(cy_bursts);
    
    for d1=1:size(block,1)
        for d2=1:size(block,2)
            if d2+1<length(block(d1,:))
                epoch=block(d1,d2):block(d1,d2+1);
                l=find(data(j,epoch)==1);
                if ~isempty(epoch(l))
                    pha_b(d1,d2)=hp(epoch(l(1))); %%% picking just the first spike in a cycle l(1)
                    all=[all hp(epoch(l(1)))];
                else
                    pha_b(d1,d2)=NaN;
                end
            end
        end
    end
    
    
    % %     n=40; (bu(1:10,1) if we want fixed number of burst with spikes per
    % cycle
    for x=1:size(pha_b,2)
        bu=pha_b(:,x); bu=bu(~isnan(bu));
        check1(1,x)=length(bu);
        vl(1,x)=circ_r(bu);
        pp(1,x)=circ_mean(bu); clear bu bu1
    end
    
    gg=gg+1;
    new=reshape(pha_b,1,size(pha_b,1)*size(pha_b,2));
    new1=new(~isnan(new));
    
    figure(1)
    subplot(1,2,1)
    polarhistogram(ang,'BinWidth',2*pi/12)
    hold on
    polarhistogram(new,'BinWidth',2*pi/12)
    subplot(1,2,2)
    polarplot([0 circ_mean(ang')], [0,circ_r(ang') ],'linewidth',6)
    hold on
    polarplot([0 circ_mean(new1')], [0,circ_r(new1') ],'linewidth',6)

end



%% -200 +200

clear all; close all
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/probe SUA_act_mat')
load('data_SUA_SNR.mat')

data=cell2mat(data_region);
ecog=[];
for i=1:size(data_region,1)
    ecog = [ecog ; repmat(Ecog_region(i,:),size(data_region{i,1},1),1)];
end

clearvars -except ecog name

srn=1000;
[b,a]=butter(2,[15/(0.5*srn) 35/(0.5*srn)],'bandpass');
tt=0;
r=0;
ep_all=[];

j=2;
clearvars -except j b a ecog data srn tt name pref_r srate in_b n_spk r in_bs in_mb in_ms vl idx ep_all
epochs_t=[];
ctx=ecog(j,:);

Ecogfiltered=filtfilt(b,a,ctx);
[maxvalM,maxidxM] = findpeaks(Ecogfiltered);
data(j,:)=zeros(1,length(Ecogfiltered));
data(j,maxidxM)=1;
env=abs(hilbert(Ecogfiltered));

[onset,offset]=bursts_aligned(env,ctx);
onset1=onset; clear onset;  onset=horzcat(onset1{:});
offset1=offset; clear offset;  offset=horzcat(offset1{:});

data_ones=find(data(j,:)==1);
hp=wrapToPi(angle(hilbert(Ecogfiltered)));
ang=hp(data_ones);

data_zeros=find(data(j,:)==0); dat_b=hp;
dat_b(data_zeros)=NaN;


if (circ_rtest(ang))<0.05
    tt=tt+1;
    pref_r(tt,:)=circ_mean(ang');
    vl(tt,:)=circ_r(ang');
    idx(tt,:)=j;
    % % %         figure(1)
    % % %         polarplot([0 circ_mean(ang')], [0, circ_r(ang')],'linewidth',2)
    % % %         hold on
    %         title(sprintf('vl %d',(j)))
    
    
    d=data(j,:);
    srn=1000;
    spkrate_1=[];
    for i =1:srn:(length(d)-srn);
        spkrate_1=[spkrate_1 numel(find(d(i:i+srn)==1))];
    end
    srate(tt,1)=mean(spkrate_1);
    
    for jj=1:length(onset)
        r=r+1;
        el=200;
        if onset(jj)>el && onset(jj)+el<length(dat_b)
            epoch=dat_b(onset(jj)-el:onset(jj)+el);
            dum=find(~isnan(epoch));
            if (numel(dum))>5
                epochs_t=[epochs_t  epoch(dum)];
                ep_all=[ep_all epoch(dum)];
            end
            n_spk(1,r)= numel(dum);
        end
    end
    in_b(j,1)=circ_mean(epochs_t');
    in_mb(j,1)=circ_r(epochs_t');
    
    figure(1)
    subplot(1,2,1)
    polarhistogram(ang,'BinWidth',2*pi/12)
    hold on
    polarhistogram(epochs_t,'BinWidth',2*pi/12)
    subplot(1,2,2)
    polarplot([0 circ_mean(ang')], [0,circ_r(ang') ],'linewidth',6)
    hold on
    polarplot([0 circ_mean(epochs_t')], [0,circ_r(epochs_t') ],'linewidth',6)
end