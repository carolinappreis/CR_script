
clear all
close all
cond={'NS';'HF';'C'};
cohort=[ 1 3 4 6];


for iii=2:4
    trial=1;
    p=[[0 0.5 0.5];[0.8 0.5 0.5];[0.8 0.5 0.8]];
    
    for co=1:3
        
        load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));
        
        [a,b]=  butter(3, [2/(0.5*samplerate2)], 'high'); %15
        filt_cs=(filtfilt(a,b,signal(3,:)));
        [Pxx,F] = pwelch(filt_cs, samplerate2, [], floor(samplerate2), samplerate2);
        
        idx=(find(F>2));
        cut=idx(1); clear idx
        idx=(find(F<20));
        cut(1,2)=idx(end); clear idx
        dum=Pxx(cut(1):cut(2));
        Peak(co,:)=find(dum==max(dum))+cut(1)-1;
        
        figure(iii-1)
        subplot(1,3,co)
        plot(tempo,filt_cs)
        
        
        ps_curves (co,:)= Pxx;
        figure(iii+3)
        plot((F),(ps_curves(co,:)./(sum(Pxx(30:50)))),'Color',p(co,:),'LineWidth',3)
        hold on
        xlim([1 10])
        xticks([2:2:10])
        % ylim([0 3.5e-4])
        box('off')
        xlabel('Frequency (Hz)')
        ylabel('Power(\muV^2)')
        legend({'NS','PLS','HFS'})
        legend('boxoff')
    end
end
