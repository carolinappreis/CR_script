clear all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb cc cond nostim
    DBS_Fpeak
    
    
    %         load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'NS_PS.mat'));
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_NS_PS.mat'))
    
    in2=1;
    
    if in2==1
        in=3;
    elseif in2==2 % other axis 1
        in=6;
    elseif in2==3 % other axis 2
        in=7;
    end
    data=SmrData.WvData;
    samplerateold=SmrData.SR;
    %%------------------------
    tremor3=data([3 6 7],:);
    
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerateold) (Fpeak+2)/(0.5*samplerateold)],'bandpass');
    
    [m,n]=butter(3,[1.5/(0.5*samplerateold)],'low');
    
    for i=1:size(tremor3,1)
        dc_t1(i,:)=filtfilt(m,n,tremor3(i,:));
        filt_t1(i,:)=filtfilt(b,a,tremor3(i,:));
        env_t1(i,:)=abs(hilbert(filt_t1(i,:)));
    end
    data1=vertcat(filt_t1,env_t1,dc_t1);
    
    Fs=(Fpeak+2)*2+5;
    time=0:1/samplerateold:(size(data1,2)-1)/samplerateold;
    ts=timeseries(data1,0:(1/samplerateold):((size(data1,2)-1)/samplerateold));
    ts1=resample(ts,0:1/Fs:((size(data1,2)-1)/samplerateold),'linear');
    data2(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    
    
    filt_t3=data2(1:3,:);
    env_t3=data2(4:6,:);
    dc_t3=data2(7:9,:);
    
    
    
    % all conditions
    s=round(([8646 88131 174501 234200 315901 393800].*Fs)./1000,0);
    e=round(([82701 147501 226501 297200 377500 454600].*Fs)./1000,0);
    
    %just posture
    
    s=round(([8646  174501  315901].*Fs)./1000,0);
    e=round(([82701  226501  377500].*Fs)./1000,0);
    
    
    % jj=3;%% axis
    %
    % s=s(pca_idx{1,1}(:)==jj);
    
    
    % %just spiral
    % s=round(([ 88131 234200  393800].*Fs)./1000,0);
    % e=round(([ 147501  297200  454600].*Fs)./1000,0);
    t=5*Fs;
    for th=1:size(tremor3,1)
        RS_raw1=[];
        RS_t1=[];
        RS_e1=[];
        RS_dc1=[];
        for tr=1:length(s)
            RS_raw1=[RS_raw1 tremor3(th,s(tr)-Fs:e(tr))];
            RS_t1=[RS_t1 filt_t3(th,s(tr)-Fs:e(tr))];
            RS_e1=[RS_e1 env_t3(th,s(tr)-Fs:e(tr))];
            RS_dc1=[RS_dc1 dc_t3(th,s(tr)-Fs:e(tr))];
        end
        
        
        RS_r(th,:)=RS_raw1;
        RS_t(th,:)=RS_t1;
        RS_e(th,:)=RS_e1;
        RS_dc(th,:)=RS_dc1;
        
    end
end

% clearvars -except Fs rs ns fs

