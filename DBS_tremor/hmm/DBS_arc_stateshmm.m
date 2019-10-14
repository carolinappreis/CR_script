%%% check start and end points and cf. with notes.

clear all
close all
iii=[1];
cc=1
numb=1;
%     :length(iii);
clearvars -except iii numb ttall ampall ph_stim LS tt1 cc
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));


in2=3; % analysing the "main tremor axis"

if in2==1
    in=3;
elseif in2==2 % other axis 1
    in=6;
elseif in2==3 % other axis 2
    in=7;
end

data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(in,:));
addon=92; addon_end=35;

time=0:1/samplerateold:(size(data,2)-1)/samplerateold;

%%% downsample
samplerate=1000;
ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
ts1=resample(ts,0:1/samplerate:((size(data,2)-1)/samplerateold),'linear');
tremor2(1:size(ts1.data,3))=ts1.data;

Fpeak=4;

if (Fpeak-2)>=1
    [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
else
    [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
end

tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;

dummy=hilbert(tremor_or);
envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));

%     load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\start_end_RS.mat')
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/start_end_RS.mat')

if cc==1
    s=round((start1*samplerate)./samplerateold,0);
    e=round((ending1*samplerate)./samplerateold,0);
    dum1= xx1;
else
    s=round((start2*samplerate)./samplerateold,0);
    e=round((ending2*samplerate)./samplerateold,0);
    dum1= xx2;
end


%     load('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\states_spiral_3ch.mat')
load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/Tnew_states2_motor',num2str(cc),'_3ch.mat'),'states')


for cr=1:3
    if ~isempty (states{cr,1})
        bg_states{cr,1}(1,:)=s(1,states{cr,1}(1,:));
        end_states{cr,1}(1,:)=e(1,states{cr,1}(1,:));
        dum{cr,1}(1,:)=dum1(1,states{cr,1}(1,:));
    end
end

for hh=1:size(bg_states,1)
    if ~isempty (bg_states{hh,1})
        tremor_or2=NaN(length(bg_states{hh,1}),1);
        tremor_or3=NaN(length(bg_states{hh,1}),1);
        
        
        for i=1:length(bg_states{hh,1})
            if (~isnan(bg_states{hh,1}(i)))
                tremor_or3(i,1)=mean(envelope(bg_states{hh,1}(i)-1000:bg_states{hh,1}(i)));
                tremor_or2(i,1)=(mean(envelope(end_states{hh,1}(i)-1000:end_states{hh,1}(i)))-mean(envelope(bg_states{hh,1}(i)-1000:bg_states{hh,1}(i))))/mean(envelope(bg_states{hh,1}(i)-1000:bg_states{hh,1}(i)));
                xx_st{hh,1}(i)= dum{hh,1}(i);
            else
                tremor_or2(i,1)=NaN;
                tremor_or3(i,1)=NaN;
                xx_st{hh,1}(i)= NaN;
            end
        end
        
        
        tt=NaN(25,12);
        amp=NaN(25,12);
        yy=xx_st{hh,1}(:);
        
        for i=1:12
            tt(1:sum(yy==i),i)=tremor_or2(find(yy==i));
            amp(1:sum(yy==i),i)=tremor_or3(find(yy==i));
            ph_stim{hh,1}(numb,i)=sum(yy==i);
        end
        
        clear yy xx_st
        
        tt1{hh,1}{numb,1}=tt;
        
        ttall {hh,1}(numb,:)=nanmedian(tt);
        ampall {hh,1}(numb,:)=nanmedian(amp);
        
    end
end
% %     clearvars -except ttall iii numb ampall ph_stim LS tt1 hh


% clearvars -except ttall ampall ph_stim LS tt1
f=figure()
for i=1:3
    subplot(1,3,i)
    if ~isempty (ttall{i,1})
        bar(ttall{i,1})
        hold on
        plot((cell2mat(tt1{i,1}))','.')
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        set(gca,'XTickLabelRotation',45)
        set(gca,'FontSize',12)
        set(f,'color','w');
        box('off')
        if i==3
            title('State1+State2','FontSize',12)
        else
            title(['State',num2str(i)],'FontSize',12)
        end
    end
end
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 30, 10];
% f=figure()
% for i=1:3
%     subplot(1,3,i)
%     if ~isempty (ttall{i,1})
%         bar(0:30:330,ttall{i,1},'LineWidth',0.5,'FaceColor','k','EdgeColor','k')
%         ylabel ('Power spectral Density')
%         xlabel ('Frequency (Hz)')
%         yline(0,'LineWidth',1)
%         % yticks([ -0.5:0.25:0.5])
%         box('off')
%         ylabel('Change in tremor severity')
%         xlabel('Stimulation phase (degrees)')
%         set(gca,'XTickLabelRotation',45)
%         set(gca,'FontSize',12)
%         set(f,'color','w');
%         if i==3
%             title('State1+State2','FontSize',12)
%         else
%             title(['State',num2str(i)],'FontSize',12)
%         end
%
%     end
% end
% f.Units = 'centimeters';
% f.OuterPosition= [10, 10, 20, 10];



for pp=1:3
    if ~isempty (ttall{pp,1})
        SO=repmat(ttall{pp,1},1,3);
        
        for ii=1:size(SO,1)
            for i=12+1:12*2
                smooth_r{pp,1}(ii,i-12)=sum(SO(ii,(i-1:i+1)))./length(SO(ii,(i-1:i+1)));
                
            end
        end
    end
end


f=figure()
for i=1:3
    subplot(1,3,i)
    if ~isempty (smooth_r{i,1})
        bar(0:30:330,smooth_r{i,1},'LineWidth',0.5,'FaceColor','k','EdgeColor','k')
        yline(0,'LineWidth',1)
        % yticks([ -0.5:0.25:0.5])
        box('off')
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        set(gca,'XTickLabelRotation',45)
        set(gca,'FontSize',12)
        set(f,'color','w');
        if i==3
            title('State1+State2','FontSize',12)
        else
            title(['State',num2str(i)],'FontSize',12)
        end
    end
end
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 30, 10];

load('/Users/Carolina/Documents/GitHub/CR_script/colour_pal.mat','m1_1','m1_2','m2_1','m2_2','blushred','stone','aegean');
col(1,:)=blushred;
col(2,:)=aegean;
f=figure()
for i=2
        bar(0:30:330,smooth_r{i,1},'LineWidth',0.5,'FaceColor',col(cc,:),'EdgeColor',col(cc,:))
        yline(0,'LineWidth',1)
        % yticks([ -0.5:0.25:0.5])
        box('off')
        ylabel('Change in tremor severity')
        xlabel('Stimulation phase (degrees)')
        set(gca,'XTickLabelRotation',45)
        set(gca,'FontSize',14)
        set(f,'color','w');
end
f.Units = 'centimeters';
f.OuterPosition= [10, 10, 12, 12];