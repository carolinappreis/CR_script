clear all

cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
load 'animal_lesion_nolesion'

state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
power_all_nc=cell(1,1);
power_all_c=cell(1,1);
check=[];
rr=1;
for nn=1:length(state_mat)
    
    cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/MATLAB/KN')
    clearvars -except nfile_i state_mat nn A power_all_c power_all_nc n rr  ts_c  ts_nc m
    name=A(state_mat(nn),1:(find(A((state_mat((nn))),:)=='.')-1))
    load(name)
    B=who;
    ctxchan=[];
    thalchan=[];
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='i')) & min(find(B{ii}=='p')) & min(find(B{ii}=='s')))
            ctxchan=ii;
        elseif ~isempty(min(find(B{ii}=='p')) & min(find(B{ii}=='r')) & min(find(B{ii}=='o')) & min(find(B{ii}=='b')))
            thalchan=[thalchan ii];
        end
    end
    
    thal_BZ=[];
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='B')) & min(find(location=='Z')));
            thal_BZ=[thal_BZ,i];
        end
    end
    
    thal_SNR=[];
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='S')) & min(find(location=='N')) &  min(find(location=='r')));
            thal_SNR=[thal_SNR,i];
        end
    end
    
    
    thal_local=thal_BZ;
    
    if ~isempty (thal_local)
        chanofinterest=thalchan(thal_local);
        
        eval(['samprateold=1/' B{ctxchan} '.interval;']);
        eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
        WaveData=double(WaveData);
        
        cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/code')
        
        ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
        ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
        WaveData_DC(1,:)=ts1.data;
        timeold=0:(1/samprateold):(length(WaveData(1,:))-1)*(1/samprateold);
        for ii=1:length(chanofinterest)
            eval(['WaveData(1+ii,:)=' B{chanofinterest(ii)} '.values;']);
            WaveData=double(WaveData);
            WaveData_DC(1+ii,:)=makemua_CR_1(WaveData(1+ii,:)',0.001,0.003,round(samprateold),1000,4);
        end
        
        clearvars muafilt_c muafilt_nc
        samprate=1000;
        
        time=0:(1/samprate):(length(WaveData_DC(1,:))-1)*(1/samprate);
        
        muafilt=[];
        [Pxx_ind,F_ind]=pwelch(WaveData_DC(1,:),samprate,[],samprate,samprate);
        Pxx_ind_beta=Pxx_ind(16:36);
        filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
        [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
        muafilt(1,1:length(time)) = filtfilt(b,a,WaveData_DC(1,:));
        
        
        power_probe_nc=[];
        power_probe_c=[];
        muafilt_c=[];
        muafilt_nc=[];
        m=2;
        for r=2:size(WaveData_DC,1)
            [power,f]=pwelch(WaveData_DC(r,:),1000,[],1000,1000);
            [Pxx_ind,F_ind]=mscohere(WaveData_DC(1,:),WaveData_DC(r,:),samprate,[],samprate,samprate);
            Pxx_ind_beta=Pxx_ind(16:36);
            filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
            [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
            P=sum(Pxx_ind(16:36))./sum(Pxx_ind(1:end));
            if P>0.1
                muafilt_c(m,1:length(time)) = filtfilt(b,a,WaveData_DC(r,:));
                power_probe_c(m,1:length(power))=(power)';
            elseif  P<0.1
                muafilt_nc(m,1:length(time)) = filtfilt(b,a,WaveData_DC(r,:));
                power_probe_nc(m,1:length(power))=(power)';              
            end
            m=m+1;
        end
        
        if ~isempty(power_probe_c)
        power_all_c{rr,:}=power_probe_c;
        end
        if ~isempty(power_probe_nc)
        power_all_nc{rr,:}=power_probe_nc;
        end
    end
    
    rr=rr+1;
end

clearvars -except power_all_c power_all_nc 

for j=1:size(power_all_c,1)
    if ~isempty(power_all_c{j})
        power_all_c{j}( ~any(power_all_c{j},2), : ) = [];
    end
end



for j=1:size(power_all_nc,1)
    if ~isempty(power_all_nc{j})
        power_all_nc{j}( ~any(power_all_nc{j},2), : ) = [];
    end
end

power_c=cell2mat(power_all_c(~cellfun('isempty',power_all_c)));
power_nc=cell2mat(power_all_nc(~cellfun('isempty',power_all_nc)));

clearvars -except power_all_c power_all_nc power_c power_nc
cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A1_Thal/mat')
%save 'BZ_power'

figure(1)
plot(power_nc(:,5:80)')
hold on
plot(median(power_nc(:,5:80)),'LineWidth',2)
figure(2)
plot(power_c(:,5:80)')
hold on
plot(median(power_c(:,5:80)),'LineWidth',2)