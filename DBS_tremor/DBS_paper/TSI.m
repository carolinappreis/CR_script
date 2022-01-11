
cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];
rng('default')
gen=(rng);
spiral=0;

if spiral==0
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture.mat');
else
    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral.mat');
end


for iii = 1:length(cohort)
    clearvars -except  cohort cond iii clust s start ending yy out gen h_up spiral freq_bl amp_bl avg val_tsi
    
    co=1;
    
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
    [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
    [peak_ax, start, ending, yy, h_up]= srtend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up,spiral);
    [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);
    
    [match_ax]=link_ax(spiral);
    [Pxx,F] = pwelch(s.raw{iii, 1} (match_ax(co,iii,1),h_up{iii,1}), samplerate, [], samplerate, samplerate);
    power=Pxx(3:end);
    find(power==max(power));
    if sum(F(2)-F(1))==1
    f=F(find(power==max(power))+2);
    else
        error(' adjust freuqency from pxx');
    end
    dat=s.raw{iii, 1} (match_ax(co,iii,1),h_up{iii,1});
    
    
    if (~isvector(dat))
        error(' Input must be univariate vector time-series.');
    end
    dat=dat(:);
    
    % Detrend
    
    rate = 1e3;
    ft_freq = f;
    ft_bandwidth = 2;
    rejectmode = 3;
    analysis = compute_tsi( dat, ft_freq, ft_bandwidth, rate, rejectmode );
    
    % TSI
    val_tsi(1,iii)=iqr(analysis.deltaIF);
    
end
