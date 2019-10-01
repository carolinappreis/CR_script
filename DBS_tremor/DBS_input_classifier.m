% %%% check start and end points and cf. with notes.
% 
% clear all
% close all
% iii=[1];
% 
% for numb=1;
%     %     :length(iii);
%     clearvars -except iii numb ttall ampall ph_stim LS tt1
%            load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
% % load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));
% 
%     in2=1; % analysing the "main tremor axis"
%     
%     DBS_find_cond;
%     clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency
% 
% 
%     for hh=2;
% 
%         handup=[];
%         for i=1:length(start{hh,1})
%             handup=[handup start{hh,1}(i):ending{hh,1}(i)]; %#ok<*AGROW>
%         end
%         handup=sort(handup,'ascend');
% 
%         [Pxx,F]=pwelch(tremor2(handup),samplerate,[],samplerate,samplerate);
% 
%         frange=F(3:10);
%         Pxxrange=Pxx(3:10);
% 
%         Fpeak=frange(find(Pxxrange==max(Pxxrange)));
% 
%         if (Fpeak-2)>=1
%             [b,a]=butter(2,[(Fpeak-2)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
%         else
%             [b,a]=butter(2,[(1)/(0.5*samplerate) (Fpeak+2)/(0.5*samplerate)],'bandpass'); %15
%         end
%         tremor_or=filtfilt(b,a,tremor2)*10*9.81/0.5;
%         % tremor_or=zscore(tremor_or);
%         dummy=hilbert(tremor_or);
%         envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
%         phase=angle(dummy);
%         frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
% 
%         tremor=(data(3,:));% %score(:,1)';%
%         ts=timeseries(tremor,0:(1/samplerateold):((size(data,2)-1)/samplerateold));
%         ts1=resample(ts,0:0.001:((size(data,2)-1)/samplerateold),'linear');
%         [b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
%         % tremor_or=zscore(tremor_or);
%         dummy=hilbert(tremor_or);
%         envelope=sqrt((real(dummy).^2)+(imag(dummy).^2));
%         phase=angle(dummy);
%         frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
%         tremor=(data(3,:));% %score(:,1)';%
% %         timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
% 
% 
%         dt1_s=floor((dt1_s./samplerateold)*samplerate)+addon;
%         dt1_e=floor((dt1_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
%         dt2_s=floor((dt2_s./samplerateold)*samplerate)+addon;
%         dt2_e=floor((dt2_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
% 
%         tremor_or2=NaN(length(start{hh,1}),1);
%         tremor_or3=NaN(length(start{hh,1}),1);
% 
%         % amplitude
%         tremor_or2=NaN(length(start{hh,1}),1);
%         for i=1:length(start{hh,1})
%             if (~isnan(start{hh,1}(i)))
%                 tremor_or3(i,1)=mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
%                 tremor_or2(i,1)=(mean(envelope(ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
%                 pha_idx(i,1)=xx{hh,1}(i);
%             else
%                 tremor_or2(i,1)=NaN;
%                 tremor_or3(i,1)=NaN;
%                 pha_idx(i,1)=NaN;
%             end
%         end
% 
%         input_2=[tremor_or2 tremor_or3 pha_idx];
% 
%     end
% 
% end

clear all
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA')
load('input.mat')
input_2=input_2(1:length(input_1),:);
motor(1:size(input_1,1),1)=0;
motor(size(input_1,1)+1:size(input_1,1)*2,1)=1;
input_all=vertcat(input_1,input_2);
input_new=[input_all motor];

amp_idx=find(input_new(:,1)>0);
sup_idx=find(input_new(:,1)<0);

amp=input_new(amp_idx,2:end);
sup=input_new(sup_idx,2:end);


inp=[amp;sup(1:size(amp,1),:)];
effect(1:size(amp,1),1)={'supr'};
effect(size(amp,1)+1:(2*size(amp,1)),1)={'amp'};

clearvars -except inp effect



% motor=cell(2*size(input_1,1),1);
% motor(1:size(input_1,1),1)={'posture'};
% motor(size(input_1,1)+1:end,1)={'spiral'};
% input_all=vertcat(input_1,input_2);


rng(1); % For reproducibility
Mdl = TreeBagger(size(amp,1),inp,effect,'OOBPrediction','On','Method','classification')

view(Mdl.Trees{1},'Mode','graph')

