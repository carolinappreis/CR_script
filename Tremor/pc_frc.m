clear all
iii=[1 2 3 4 5 8 10 11 12 13];
for in2=1:3;
    for numb=1:length(iii);
        clearvars -except iii numb frc in2
        %     load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
        load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))
        
        % analysing the "main tremor axis"
        start_cleaner;
        
        clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency
        
        handup=[];
        for i=1:length(start)
            handup=[handup start(i):ending(i)]; %#ok<*AGROW>
        end
        handup=sort(handup,'ascend');
        
        [Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);
        
        frange=F(3:10);
        Pxxrange=Pxx(3:10);
        
        Fpeak=frange(find(Pxxrange==max(Pxxrange)));
        
        if (Fpeak-2)>=1
            [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
        else
            [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
        end
        tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
        
        dummy=hilbert(tremor_or);
        envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
        phase=angle(dummy);
        frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
        
        tremor=(data(3,:));% %score(:,1)';%
        ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
        ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
        tremorx(1:size(ts1.data,3))=ts1.data;
        tremor=(data(5,:));% %score(:,1)';%
        ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
        ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
        tremory(1:size(ts1.data,3))=ts1.data;
        tremor=(data(6,:));% %score(:,1)';%
        ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
        ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
        tremorz(1:size(ts1.data,3))=ts1.data;
        if (Fpeak-2)>=1
            [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
        else
            [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
        end
        tremorxf=filtfilt(b,a,tremorx);
        tremoryf=filtfilt(b,a,tremory);
        tremorzf=filtfilt(b,a,tremorz);
        
        
        for j=1:length(start)
            x=[tremorxf(start(j):ending(j));tremoryf(start(j):ending(j));tremorzf(start(j):ending(j))];
            [pc,score,latent,tsquare] = pca(x');
            xxx(j,1:3)=pc(1:3,1);
            ma(j)=(find(abs(xxx(j,1:3))==max(abs(xxx(j,1:3)))));
        end
        
        tremor_or2=NaN(20,5001);
        tremor_or22=NaN(20,5001);
        
        if in2==1
            for i=1:length(start)
                if (~isnan(start(i))&& ma(i)==1)
                    tremor_or2(i,1:(ending(i)-start(i)+1))=unwrap(phase(start(i):ending(i)));
                    tremor_or22(i,1:(ending(i)-start(i)+1))=(phase(start(i))+(0:1:(ending(i)-start(i)))*2*pi/(1000./mean(frequency(start(i)-1000:start(i)))));
                    tremor_k(i,1)= (tremor_or2(i,(ending(i)-start(i)+1))-tremor_or22(i,(ending(i)-start(i)+1)))/(2*pi*0.001*(ending(i)-start(i))); %mean(frequency(ending(i)-1000:ending(i)));%
                else
                    tremor_or22(i,1:5001)=NaN;
                    tremor_or2(i,1:5001)=NaN;
                    tremor_k(i,1)=NaN;
                end
            end
        else
            for i=1:length(start)
                if ~isnan(start(i))
                    tremor_or2(i,1:(ending(i)-start(i)+1))=unwrap(phase(start(i):ending(i)));
                    tremor_or22(i,1:(ending(i)-start(i)+1))=(phase(start(i))+(0:1:(ending(i)-start(i)))*2*pi/(1000./mean(frequency(start(i)-1000:start(i)))));
                    tremor_k(i,1)= (tremor_or2(i,(ending(i)-start(i)+1))-tremor_or22(i,(ending(i)-start(i)+1)))/(2*pi*0.001*(ending(i)-start(i))); %mean(frequency(ending(i)-1000:ending(i)));%
                else
                    tremor_or22(i,1:5001)=NaN;
                    tremor_or2(i,1:5001)=NaN;
                    tremor_k(i,1)=NaN;
                end
            end
        end
        
        
        %     idx_outl=find(tremor_k>(mean(tremor_k+2*(std(tremor_k))))|tremor_k<(mean(tremor_k-2*(std(tremor_k)))));
        %     tremor_k(idx_outl,1)=NaN;
        %     xx(1,idx_outl)=NaN;
        
        clear tt
        tt=[];
        k=1;
        
        tt=NaN(20,12);
        
        for i=1:12
            tt(1:sum(xx==i),i)=tremor_k(find(xx==i));
        end
        frc(numb,:)=nanmedian(tt);
        
        cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
        if in2==1
            fm_ax=frc;
            save('fm_ax.mat','fm_ax')
            
        elseif in2==2
            s_frc=frc;
            save('s_frc.mat','s_frc')
            
        elseif in2==3
            t_frc=frc;
            save('t_frc.mat','t_frc')
        end
    end
end
