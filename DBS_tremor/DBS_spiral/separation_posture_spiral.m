%% NS
%%% patient 1 - joingning two files with NS

% clear
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P01/01_nostim_baseline.mat')
% dum=SmrData; clear SmrData;
% load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P01/01_nostim_pohuureandspiral.mat')
% dum2=SmrData; clear SmrData.WvData;
% SmrData.WvData=[dum.WvData dum2.WvData];
% if(length(dum.WvData) + length(dum2.WvData)) ==length(SmrData.WvData)
% clearvars -except SmrData
% cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA')
% save ('P01_NS')
% else
% error('error')
% end

%%% comments:
% pateint 1 ---- baseline attacjed to beginning very different from baseline pohuure spiral - posture 2 and 3 are dodgy 


clear
cohort =[ 1 3 4 6];
for iii=4
%     1:size(cohort,2)
    clearvars -except cohort iii
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA')
load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_NS.mat'))
[d]=dbs_preprocess(SmrData); 
tremor3=d.data_ds;
[p,q]=butter(3,[5/(0.5*d.samplerateold)],'low');
[m,n]=butter(3,[1.5/(0.5*d.samplerateold)],'low');
for i=1:size(tremor3,1)
    t2(i,:)=filtfilt(p,q,tremor3(i,:));
    t3(i,:)=filtfilt(m,n,tremor3(i,:));
end

%% plot to asjuhu NS_BE_P NS_BE_S
% % subplot(3,1,1)
% % plot(tremor3(1,:))
% % subplot(3,1,2)
% % plot(t2(1,:))
% % subplot(3,1,3)
% % plot(t2(3,:))



% NS_BE_P %code with indicies for beggining end points of NS
NS_BE_S

%%%check NS_BE_P
close all

plot(tremor3(1,:))
hold on
for i=1:size(hu{iii,1}(1,:),2)
xline(hu{iii,1}(i),'r')
xline(hd{iii,1}(i),'r')
end

end
plot3(1:length(t3(1,hu{iii,1}(i):hd{iii,1}(i))),t3(2,hu{iii,1}(i):hd{iii,1}(i)),t3(3,hu{iii,1}(i):hd{iii,1}(i)))


%%%% -------------

%% RS
clear
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P01/01_random_postureandspiral.mat')
dum=SmrData.WvData; 
%%% clear SmrData.WvData; MANUALLY
dum2=dum(:,2998011:end);
SmrData.WvData=dum2;
clearvars -except SmrData
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA')
% % % save ('P01_RS')

%%%%%-------

%% PLS posture

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P03/P03_PLS.mat')
cut=6535190;
data=SmrData.WvData([3 6 7],:);
plot3(1:length(data(2,1:cut)),data(2,1:cut),data(3,1:cut))
plot3(1:length(data(2,cut:end)),data(2,cut:end),data(3,cut:end))
dum=SmrData.WvData; 
%%% clear SmrData.WvData; MANUALLY
dum2=dum(:,1:cut);
SmrData.WvData=dum2;
clearvars -except SmrData
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd')
% % % save ('P03_PLS_P')

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P06/P06_120POSTURE3TIMES_240posture2times.mat')
cut=4982100;
% data=SmrData.WvData([3 6 7],:);
% plot3(1:length(data(2,1:cut)),data(2,1:cut),data(3,1:cut))
% plot3(1:length(data(2,cut:end)),data(2,cut:end),data(3,cut:end))
dum=SmrData.WvData; 
%%% clear SmrData.WvData; MANUALLY
dum2=dum(:,1:cut);
SmrData.WvData=dum2;
clearvars -except SmrData
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd')
% % % save ('P06_PLS_P1')


load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P06/P06_120POSTURE3TIMES_240posture2times.mat')
cut=4982100;
% data=SmrData.WvData([3 6 7],:);
% plot3(1:length(data(2,1:cut)),data(2,1:cut),data(3,1:cut))
% plot3(1:length(data(2,cut:end)),data(2,cut:end),data(3,cut:end))
dum=SmrData.WvData; 
%%% clear SmrData.WvData; MANUALLY
dum2=dum(:,cut:end);
SmrData.WvData=dum2;
clearvars -except SmrData
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd')
% % % save ('P06_PLS_P2')


%%%PLS spiral

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P03/P03_PLS.mat')
cut=7418625;
data=SmrData.WvData([3 6 7],:);
plot3(1:length(data(2,1:cut)),data(2,1:cut),data(3,1:cut))
plot3(1:length(data(2,cut:end)),data(2,cut:end),data(3,cut:end))
dum=SmrData.WvData; 
%%% clear SmrData.WvData; MANUALLY
dum2=dum(:,cut:end);
SmrData.WvData=dum2;
clearvars -except SmrData
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_spiral')
% % % save ('P03_PLS_S')

load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P04/p04_pls_90 spiral_3posture300_2spiral.mat')
% data=SmrData.WvData([3 6 7],:);
% plot3(1:length(data(2,1:cut)),data(2,1:cut),data(3,1:cut))
% plot3(1:length(data(2,cut:end)),data(2,cut:end),data(3,cut:end))
dum=SmrData.WvData; 
%%% clear SmrData.WvData; MANUALLY
dum2=dum(:,[1:1519204 7309372:end]);
SmrData.WvData=dum2;
clearvars -except SmrData
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_spiral')
% % %  save ('P04_PLS_S')


load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/raw_patient_data/P06/p06_60spiral_2.mat')
clearvars -except SmrData
cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_spiral')
% % % save ('P06_PLS_S')


%%% HFS



load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf/P03_HFS.mat')
% dum=SmrData.WvData; 
% % dum2=dum(:,[8929:73906]);
% dum2=dum(:,[95161:157547]);
% SmrData.WvData=dum2;
% clearvars -except SmrData
% cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf')
% save ('P03_HFS_S')

clear
load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf/P04_HFS.mat')

% dum=SmrData.WvData; 
% % dum2=dum(:,[42667:93145]);
% dum2=dum(:,[164174: 231805]);
% SmrData.WvData=dum2;
% clearvars -except SmrData
% cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf')
% save ('P04_HFS_S')


clear
load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf/P06_HFS.mat')
% dum=SmrData.WvData; 
% %  dum2=dum(:,[1:98244]);
% dum2=dum(:,[142181:259995]);
% SmrData.WvData=dum2;
% clearvars -except SmrData
% cd('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_hf')
% save ('P06_HFS_S')

[d]=dbs_preprocess(SmrData); 
tremor3=d.data_ds;
[p,q]=butter(3,[5/(0.5*d.samplerateold)],'low');
[m,n]=butter(3,[1.5/(0.5*d.samplerateold)],'low');
for i=1:size(tremor3,1)
    t2(i,:)=filtfilt(p,q,tremor3(i,:));
    t3(i,:)=filtfilt(m,n,tremor3(i,:));
end








