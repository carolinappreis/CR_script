%%% check start and end points and cf. with notes.

clear all
close all
iii=[1];

for numb=1;
    %     :length(iii);
    clearvars -except iii numb ttall amp_1b amp_1l ph_stim LS tt1
    
    in2=3; % analysing the "main tremor axis"
%     load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));
      load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))

    
    DBS_find_cond;
    clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency
    
    
    for hh=1:2;
        
        handup=[];
        for i=1:length(start{hh,1})
            handup=[handup start{hh,1}(i):ending{hh,1}(i)]; %#ok<*AGROW>
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
        % tremor_or=zscore(tremor_or);
        dummy=hilbert(tremor_or);
        envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
        phase=angle(dummy);
        frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
        
        
 
        
        
        
        for i=1:length(start{hh,1})
            if (~isnan(start{hh,1}(i)))
                amp_bef(i,:)=(envelope(start{hh,1}(i)-1000:start{hh,1}(i)+5000-1));
                xx{hh,1}(i)= xx{hh,1}(i);
            else
                amp_bef(i,:)=NaN;
                xx{hh,1}(i)= NaN;
            end
        end
        if hh==1 
        dr(1,:)=find(xx{hh,1}==4);
        else
        dr(1,:)=find(xx{hh,1}==7);
        end

%         if hh==1 
%         dr(1,:)=find(xx{hh,1}==2);
%         else
%         dr(1,:)=find(xx{hh,1}==8);
%         end
%         
%         
        for i=1:length(dr)
        env_ep{hh,1}=amp_bef(dr,:);
        end
        clear dr
     load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\HFS_median.mat','hfs');

   
        
f=figure;
plot(median(env_ep{hh,1},1),'LineWidth',2,'Color','k')
hold on
yline(hfs(in2,hh),'LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5]);
box('off')
ylabel('Instantaneous tremor severity (m/s^2)')
xlabel('Time (seconds)')   
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 12, 12];
set(gca,'FontSize',14)
set(f,'color','w');     
    end
end
