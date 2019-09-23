% % %%% Create smooth data and save
% %
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------AMP
% % %%% loading all ARC axis
clear all
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\amp_ARC')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\am_ax')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\s_arc')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\t_arc')

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/amp_ARC')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/am_ax')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/s_arc')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/t_arc')

SO=repmat(ttall,1,3);
S1=repmat(am_ax,1,3);
S2=repmat(s_arc,1,3);
S3=repmat(t_arc,1,3);
 for ii=1:size(SO,1)
    for i=size(ttall,2)+1:size(ttall,2)*2
        smo_s(ii,i-12)=sum(SO(ii,(i-1:i+1)))./length(SO(ii,(i-1:i+1)));
        smo_s1(ii,i-12)=sum(S1(ii,(i-1:i+1)))./length(S1(ii,(i-1:i+1)));
        smo_s2(ii,i-12)=sum(S2(ii,(i-1:i+1)))./length(S2(ii,(i-1:i+1)));
        smo_s3(ii,i-12)= sum(S3(ii,(i-1:i+1)))./length(S3(ii,(i-1:i+1)));
    end
end

clearvars -except smo_s smo_s1 smo_s2 smo_s3
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')

save('smooth_arc3.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------FREQ

% %%% loading all FRC axis
clear all
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\freq_FRC')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\fm_ax')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\s_frc')
% load('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\t_frc')

load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/freq_FRC')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/fm_ax')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/s_frc')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/t_frc')


SO=repmat(ttall,1,3);
S1=repmat(fm_ax,1,3);
S2=repmat(s_frc,1,3);
S3=repmat(t_frc,1,3);


for ii=1:size(SO,1)
    for i=size(ttall,2)+1:size(ttall,2)*2
        smo_s(ii,i-12)=sum(SO(ii,(i-1:i+1)))./length(SO(ii,(i-1:i+1)));
        smo_s1(ii,i-12)=sum(S1(ii,(i-1:i+1)))./length(S1(ii,(i-1:i+1)));
        smo_s2(ii,i-12)=sum(S2(ii,(i-1:i+1)))./length(S2(ii,(i-1:i+1)));
        smo_s3(ii,i-12)= sum(S3(ii,(i-1:i+1)))./length(S3(ii,(i-1:i+1)));
    end
end

clearvars -except smo_s smo_s1 smo_s2 smo_s3
cd('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data')
% cd('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data')
save('smooth_frc3.mat')

%
%

