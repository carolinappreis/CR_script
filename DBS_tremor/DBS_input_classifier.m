% %%% check start and end points and cf. with notes.
% 
% clear all
% close all
% iii=[1];
% 
% for numb=1;
%     %     :length(iii);
%     clearvars -except iii numb ttall ampall ph_stim LS tt1
%     load(strcat('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA\0',num2str(iii(numb)),'_RS_PS.mat'))
%     % load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/0',num2str(iii(numb)),'_RS_PS.mat'));
%     
%     in2=1; % analysing the "main tremor axis"
%     
%     DBS_find_cond;
%     clear handup Pxx F frange Pxxrange Fpeak tremor_or dummy envelope phase frequency
%     
%     
%     for hh=1:2;
%         
%         clear state_idx ef1 ef2 ef3 ef4 tbl phase_s
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
%         z_env=zscore(envelope);
%         diff_zenv=diff(z_env);
%         phase=angle(dummy);
%         frequency=(smooth((1000/(2*pi))*diff(unwrap(angle(dummy))),500))';
%         tremor=(data(3,:));% %score(:,1)';%
%         %         timeor=0:1/samplerate:(size(tremorx,2)-1)/samplerate;
%         
%         
%         dt1_s=floor((dt1_s./samplerateold)*samplerate)+addon;
%         dt1_e=floor((dt1_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
%         dt2_s=floor((dt2_s./samplerateold)*samplerate)+addon;
%         dt2_e=floor((dt2_e./samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);
%         
%         amp_1_m=NaN(length(start{hh,1}),1);
%         change=NaN(length(start{hh,1}),1);
%         amp_1_dir=NaN(length(start{hh,1}),1);
%         amp_1_var=NaN(length(start{hh,1}),1);
%         amp_stim_dir=NaN(length(start{hh,1}),1);
%         amp_stim_var=NaN(length(start{hh,1}),1);
%         
%         %         freq_1sb=NaN(length(start{hh,1}),1);
%         %         freq_dsb=NaN(length(start{hh,1}),1);
%         %         freq_ddur=NaN(length(start{hh,1}),1);
%         %         freq_adur=NaN(length(start{hh,1}),1);
%         
%         tremor_or2=NaN(length(start{hh,1}),1);
%         tremor_or3=NaN(length(start{hh,1}),1);
%         tremor_k=NaN(length(start{hh,1}),1);
%         
%         for i=1:length(start{hh,1})
%             if (~isnan(start{hh,1}(i)))
%                 change(i,1)=(mean(envelope(ending{hh,1}(i)-1000:ending{hh,1}(i)))-mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i))))/mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
%                 amp_1_m(i,1)=mean(envelope(start{hh,1}(i)-1000:start{hh,1}(i)));
%                 amp_1_dir(i,1)=z_env(start{hh,1}(i))-z_env(start{hh,1}(i)-1000);
%                 amp1_var(i,1)=std(diff(z_env(start{hh,1}(i)-1000:start{hh,1}(i))));
%                 amp_stim_dir(i,1)=z_env(ending{hh,1}(i))-z_env(start{hh,1}(i));
%                 amp_stim_var(i,1)=std(diff(z_env(start{hh,1}(i):ending{hh,1}(i))));
%                 
%                 tremor_or2(i,1:(ending{hh,1}(i)-start{hh,1}(i)+1))=unwrap(phase(start{hh,1}(i):ending{hh,1}(i)));
%                 tremor_or22(i,1:(ending{hh,1}(i)-start{hh,1}(i)+1))=(phase(start{hh,1}(i))+(0:1:(ending{hh,1}(i)-start{hh,1}(i)))*2*pi/(1000./mean(frequency(start{hh,1}(i)-1000:start{hh,1}(i)))));
%                 tremor_k(i,1)= (tremor_or2(i,(ending{hh,1}(i)-start{hh,1}(i)+1))-tremor_or22(i,(ending{hh,1}(i)-start{hh,1}(i)+1)))/(2*pi*0.001*(ending{hh,1}(i)-start{hh,1}(i))); %mean(frequency(ending{hh,1}(i)-1000:ending{hh,1}(i)));%
%                 motor(i,1)=hh;
%                 
%                 pha_idx(i,1)=xx{hh,1}(i);
%                 
%                 
%             else
%                 change(i,1)=NaN;
%                 amp_1_m(i,1)=NaN;
%                 amp_1_dir(i,1)=NaN;
%                 amp1_var(i,1)=NaN;
%                 amp_stim_dir(i,1)=NaN;
%                 amp_stim_var(i,1)=NaN;
%                 tremor_k(i,1)=NaN;
%                 pha_idx(i,1)=NaN;
%                 motor(i,1)=NaN;
%             end
%         end
%         
% %         for i=1:size(pha_idx,1)
% %             if pha_idx(i,1)>=1 & pha_idx(i,1)<4
% %                 phase_s(i,1)=1;
% %             elseif pha_idx(i,1)>=4 & pha_idx(i,1)<=6
% %                 phase_s(i,1)=2;
% %             elseif pha_idx(i,1)>=7 & pha_idx(i,1)<=9
% %                 phase_s(i,1)=3;
% %             elseif pha_idx(i,1)>=10 & pha_idx(i,1)<=12
% %                 phase_s(i,1)=4;
% %             end
% %         end
% %         
%         phase_s=pha_idx;
%         input=[amp_1_m phase_s motor];
%         
%         th=[prctile(change,25) prctile(change,75)];
%         state_idx{1,1}=input(find(change<th(1)),:);
%         state_idx{2,1}=input(find(change>th(1) & change<th(2)),:);
%         state_idx{3,1}=input(find(change>th(2)),:);
%         
%         
%         ef1=ones(numel(find(change<th(1))),1);
%         ef2(1:numel(find(change>th(1) & change<th(2))),1)=2;
%         ef3(1:numel(find(change>th(2))),1)=3;
%         
%         if hh==1;
%             ef_1=[ef1; ef2; ef3];
%             tbl1=cell2mat(state_idx);
%             
%         else
%             ef_2=[ef1; ef2; ef3];
%             tbl2=cell2mat(state_idx);
%             
%         end
%         
%     end
%     obs=[tbl1;tbl2];
%     effects=[ef_1;ef_2];
% end
% clearvars -except obs effects ef_1 ef_2
% 
% stim_phase=categorical(obs(:,2));
% motor_state=categorical(obs(:,3));
% amp_1sec=obs(:,1);
% 
% 
% X=table(amp_1sec, motor_state, stim_phase);
% 
% 
% 
% close all
% 
% rng(1); % For reproducibility
% Mdl = TreeBagger(90,X,effects,'OOBPrediction','on','Method','classification','OOBPredictorImportance','on')
% 
% view(Mdl.Trees{1},'Mode','graph')
% 
% figure;
% subplot(2,1,1)
% oobErrorBaggedEnsemble = oobError(Mdl);
% plot(oobErrorBaggedEnsemble)
% xlabel 'Number of grown trees';
% ylabel 'Out-of-bag classification error';
% 
% imp = Mdl.OOBPermutedPredictorDeltaError;
% 
% subplot(2,1,2)
% bar(imp);
% title('Curvature Test');
% ylabel('Predictor importance estimates');
% xlabel('Predictors');
% h = gca;
% 
% 
% 
% 
% 
% 
% 



clear all
cd('C:\Users\creis\OneDrive - Nexus365\Phasic_DBS\patient data\DBS_DATA')
load('input.mat')
input_2=input_2(1:length(input_1),:);
motor(1:size(input_1,1),1)=0;
motor(size(input_1,1)+1:size(input_1,1)*2,1)=1;
input_all=vertcat(input_1,input_2);
input_new=[input_all motor];

eff1=find(input_new(:,1)<prctile(input_new(:,1),25));
eff2=find(input_new(:,1)>prctile(input_new(:,1),25) & input_new(:,1)<prctile(input_new(:,1),75));
eff3=find(input_new(:,1)>prctile(input_new(:,1),75));


ef1=input_new(eff1,2:end);
ef2=input_new(eff2,2:end);
ef3=input_new(eff3,2:end);

obs=vertcat(ef1, ef2 ,ef3);

res(1:size(ef1,1),1)={'supr'};
res2(1:size(ef2,1),1)={'NA'};
res3(1:size(ef3,1),1)={'amp'};

R=vertcat(res,res2,res3);

R=categorical(R);
amp_1s=obs(:,1);
phase=categorical(obs(:,2));
motor=categorical(obs(:,3));


X=table(amp_1s,phase,motor);



% motor=cell(2*size(input_1,1),1);
% motor(1:size(input_1,1),1)={'posture'};
% motor(size(input_1,1)+1:end,1)={'spiral'};
% input_all=vertcat(input_1,input_2);


rng(1); % For reproducibility
Mdl = TreeBagger(70,X,R,'OOBPrediction','On','Method','classification','OOBPredictorImportance','On')

view(Mdl.Trees{1},'Mode','graph')

figure;
oobErrorBaggedEnsemble = oobError(Mdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';


imp = Mdl.OOBPermutedPredictorDeltaError;
bar(imp);
title('Curvature Test');
ylabel('Predictor importance estimates');
xlabel('Predictors');
h = gca;
