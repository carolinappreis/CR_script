

clear all
close all
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];
rng('default')
gen=(rng);
spiral=1;

if spiral==0
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
else
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
end

co=1;

for iii=1:4
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
    [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
    
    tremor_ds=d.data_ds;
    
    
    decide=0; %0= freq of filtering is task specific else is coming from prosture
    if decide==0
        if spiral==0
            NS_BE_P
        else
            NS_BE_S
        end
    else
        NS_BE_P
    end
    
    start{iii,co}=hu{iii,:};
    ending{iii,co}=hd{iii,:};
    
    
    % when patient's hand is up
    handup = [];
    for i = 1:length(start{iii,co})
        handup = [handup start{iii,co}(i):ending{iii,co}(i)]; %#ok<*AGROW>
    end
    clear i
    handup = sort(handup,'ascend');
    
    p1=figure(1)
    for aa = 1:3
        [Pxx,F] = pwelch(tremor_ds(aa,handup), samplerate, [], samplerate, samplerate);
        %          [Pxx,F] = pwelch(tremor_ds(aa,:), samplerate, [], samplerate, samplerate);
        
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;
        
        subplot(1,4,iii)
        plot(F(3:10),Pxx(3:10))
        hold on
        box('off')
        xlabel({'Frequency (Hz)'})
        ylabel({'Power (µV^2)'})
        set(gca,'FontSize',12)
        set(p1,'color','w');
        title(sprintf('patient %d',(iii)))
        clear Pxx F
    end
end
p1.Units = 'centimeters';
p1.OuterPosition= [1,100,1000,300];