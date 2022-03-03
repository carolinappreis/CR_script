
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; effect=[]; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];
rng('default')
gen=(rng);
res=10;

for n=1:2
    clearvars -except  effect res n cohort cond iii clust s start ending yy out gen h_up spiral freq_bl amp_bl avg power FWHM
    if n==1
        spiral=0;
    else
        spiral=1;
    end
    for iii = 1:length(cohort)
        clearvars -except effect  res  n cohort cond iii clust s start ending yy out gen h_up spiral freq_bl amp_bl avg power FWHM
        
        co=1;
        load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
        [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
        [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
        [match_ax]=link_ax(spiral);
        dat=s.raw{iii, 1} (match_ax(co,iii,1),h_up{iii,1});
        [Pxx,F] = pwelch(dat, res*samplerate, [], res*samplerate, samplerate);
        
        cut_1=find(F==3);
        cut_2=find(F==50);
        % take power and frequency from 3-50Hz;
        power1=smooth(Pxx(cut_1:cut_2));
        F1=(F(cut_1:cut_2))';
        peak_idx=find(power1==(max(power1)));
        half_max=max(power1)./2;
        [val,idx]=sort((abs(power1-half_max)),'ascend');
        if idx(1)-idx(2)==1
            st=idx(1);
            ed=idx(3);
        else
            st=idx(1);
            ed=idx(2);
        end
        
        FWHM(n,iii)=abs((st-ed)/res);
        figure(n)
        subplot(1,4,iii)
        plot(F1,power1)
        hold on
        yline(half_max,'r','LineWidth',2)
        plot(F1(st),power1(st),'rd', 'MarkerFaceColor','green')
        plot(F1(ed),power1(ed),'rd', 'MarkerFaceColor','green')
        xlim([3 10])
        box('off')
        title(sprintf('FWHM = %d',FWHM(n,iii)),'FontSize',12)
        
        %
        %     find((diff(below))>1)
        %     plot(time, power1(1,:))
        %     hold on
        %     plot(time(above),power1(1,above),'.')
        %
        if n==1
            load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
            dns=prctile((squeeze(out.mod_amp{iii,1}(match_ax(1,iii,1),:))),0.2083);
            d1=nanmedian(squeeze(out.mod_amp{iii,2}{match_ax(2,iii,1),1}));
            effect(1,iii)=-(min(d1)-dns);
        end
        
    end
    
end

clearvars -except FWHM effect

for n=1:2
        if n==1
        spiral=0;
    else
        spiral=1;
        end
    
load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','blushred','aegean','stone','squash','sapphire','azure');
color_b1=[aegean;stone;blushred];
x=FWHM(n,:)';

f1=figure(n+2)
if spiral==0
    y=(effect)';
else
    load('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper/all_cond_spiral.mat','Ppeak')
    y=abs((Ppeak(:,3)-Ppeak(:,1))./-Ppeak(:,1));
    x=(x(2:end));
end
[fitresult] = fit( x, y, 'poly1' );
h=plot(fitresult,x,y);
legend('off')
hold on
plot(x,y,'k.','MarkerSize',10,'HandleVisibility','off');
[r,p]=corrcoef(x,y); r=vpa(round(r,2));p=vpa(round(p,2));
title(sprintf('r = %s',r),'FontSize',12)
box('off')

xlabel('FWHM')
% ylabel(sprintf( arc.label{ii,1}))
ylabel('tremor suppresion')

set(gca,'FontSize',12)
f1.Units ='centimeters';
f1.OuterPosition= [10, 10, 8, 10];
set(f1,'color','w');


end
