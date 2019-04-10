%%% ME1_LFP is not a good recording (wideband recording done with glass
%%% electrode) 

clear all
%cd ('/Users/Carolina/Documents/MATLAB/SUA/Juxta SUA:act:mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
% load 'animal_region.mat' % comes from code pre_filegen_SUA_act
load 'non_repeat_animal_region.mat'
region_mat=thal_VM

for nn=1:length(region_mat);
    
    
    clearvars -except region_mat nn A WaveData_DC filtall dif_angs meanang angindex
    name=A(region_mat(nn),1:(find(A((region_mat((nn))),:)=='.')-1))
    cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat')
    load(name)
    B=who;
    
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='i')) & min(find(B{ii}=='p')) & min(find(B{ii}=='s')) & min(find(B{ii}=='i')))
            ctxchan=ii;
        elseif ~isempty(min(find(B{ii}=='M')) & min(find(B{ii}=='E')) & min(find(B{ii}=='1')) & min(find(B{ii}=='_')))
            thalchan=ii;
        end
    end
    
    eval(['samprateold=1/' B{ctxchan} '.interval;']);
    eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
    WaveData=double(WaveData);
    ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
    ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
    WaveData_DC(nn,1,:)=ts1.data;
    
    %   cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/code')
    cd ('\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A2_Thal\code')
    Fs=1000;
    eval(['WaveThal=' B{thalchan} '.values;']);
    WaveThal=double(WaveThal);
    WaveData_DC(nn,2,:)=makemua_CR_1(WaveThal,0.001,0.003,round(samprateold),1000,4);
    
    
    filtrange=25;
    [b,a]=butter(2,[(filtrange-5)/(0.5*Fs) (filtrange+5)/(0.5*Fs)],'bandpass');
    filtall(nn,1,:)=filtfilt(b,a,WaveData_DC(nn,1,:));
    filtall(nn,2,:)=filtfilt(b,a,WaveData_DC(nn,2,:));
    
    non_norm=squeeze(angle(hilbert(filtall(nn,1,:)))-angle(hilbert(filtall(nn,2,:))))';
    for x =1:size(non_norm,2)
                if non_norm(1,x)>pi
                    non_norm(1,x)=non_norm(1,x)-(2.*pi);
                elseif non_norm(1,x)<-pi
                    non_norm(1,x)=non_norm(1,x)+(2.*pi);
                end
    end
    
    dif_angs(nn,:)=non_norm;
    angindex(nn,1)=abs(sum(exp(sqrt(-1)*(dif_angs(nn,:))))./(length(filtall)));
    meanang(nn,1)=angle(sum(exp(sqrt(-1)*(dif_angs(nn,:))))./(length(filtall)));
end

