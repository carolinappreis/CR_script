clear all
iiii=[ 2 3 4 5 8 10 11 13 16 17];

for numb=1:length(iiii)
    clearvars -except   iiii numb  WaveData_DC all_pxx
    %     load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iiii(numb)),'_RS.mat'))
    load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iiii(numb)),'_RS.mat'))
    
    data=SmrData.WvData;
    
    samplerateold=SmrData.SR;
    ts=timeseries(data,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
    ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
    ds_data(1:size(ts1.data,1),1:size(ts1.data,3))=ts1.data;
    samplerate=1000;
    tt=0:1/samplerate:(size(ds_data,2)-1)/samplerate;
    tre_3=ds_data([3 5 6],:);
    emg=[7,8];
    
    index=[];
    for i=2:size(data,2)-1
        if data(2,i-1)<2.5 && data(2,i)>2.5
            index=[index i];
        end
    end
    clear i
    
    indexes4=[index(1) index(find(diff(index)./samplerateold > 0.95)+1)];
    indexes3=[index(find(diff(index)./samplerateold > 0.95)) index(end)];
    
    dd2=round(data(4,:)*100)./100;
    for i=1:length(indexes4)
        xx(i)=round(dd2(indexes4(i))./0.1); %#ok<*SAGROW>
    end
    clear i
    
    start=floor((indexes4./samplerateold)*samplerate);
    ending=floor((indexes3./samplerateold)*samplerate);%floor(5*samplerate);
    
    %%% when patient's hand is up
    handup=[];
    for i=1:length(start)
        handup=[handup start(i):ending(i)]; %#ok<*AGROW>
    end
    clear i
    handup=sort(handup,'ascend');
    
    
    for ii=1:length(emg)
        
%         zdata(ii,:)=zscore(data(emg(ii),:));
        
        WaveData_DC(numb,ii,:)=makemua_emg(data(emg(ii),:)',0.008,0.012,round(samplerateold),1000,4);
        
%         [b,a]=butter(2,[3/(0.5*samplerate) 80/(0.5*samplerate)],'bandpass'); %15
%         emg_zf(ii,:)=filtfilt(b,a,zdata(ii,:));
        [Pxx,F]=pwelch(squeeze(WaveData_DC(numb,ii,handup)),samplerate,[],samplerate,samplerate);
        all_pxx(numb,ii,:)=Pxx;
        
    end
    
end