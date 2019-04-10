clear all
% cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
load list_juxta; A=SUA_juxta;
newfile=[];
nolesion=[];
lesion=[];
newfile=size(A,1);
for nfile_i=1:(newfile)
    name=A(nfile_i,1:(find(A(nfile_i,:)=='.')-1));
    load(name)
    if IpsiEEG.dopamine=='6-OHDA'
        lesion=[lesion nfile_i];

    else nolesion=[nolesion nfile_i];

    end
end

clearvars -except nolesion lesion A
% save 'SUAjuxta_lesion_nolesion'
% 
% clear all
% %cd ('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat')
% 
% load ('SUAjuxta_lesion_nolesion.mat','lesion','A')

state_mat=lesion;% lesion vs nolesioned
newfile=size(A,1);
thal_VA=[];
thal_VM=[];
thal_VL=[];

for nn=1:length(state_mat);  %%% chose state (lesion OR nolesion)
   
    name=A(state_mat(nn),1:(find(A((state_mat((nn))),:)=='.')-1))
    load(name)
    B=who;
    
    
    ctxchan=[];
    thalchan=[];
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='i')) & min(find(B{ii}=='p')) & min(find(B{ii}=='s')))
            ctxchan=ii;
        elseif ~isempty(min(find(B{ii}=='u')) & min(find(B{ii}=='n')) & min(find(B{ii}=='i')) & min(find(B{ii}=='t')))
            thalchan=[thalchan ii];
        end
    end

        eval(['location=' B{thalchan} '.location']);
        if ~isempty(min(find(location=='V')) & min(find(location=='A')))
            thal_VA=[thal_VA,nn];
        elseif ~isempty(min(find(location=='V')) & min(find(location=='M')))
            thal_VM=[thal_VM,nn];
        elseif ~isempty(min(find(location=='V')) & min(find(location=='L')))
            thal_VL=[thal_VL,nn];         
     
    end
end

clearvars -except thal_VA thal_VM thal_VL A lesion
thal_VA=lesion(thal_VA);
thal_VM=lesion(thal_VM);
thal_VL=lesion(thal_VL);
% save 'animal_region'


clear all
%cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act

thal_VA=thal_VA(1,2:3);
thal_VL=thal_VL(1,[4 5 7]);
thal_VM=thal_VM(1,6:7);
% save 'non_repeat_animal_region.mat'

clear all
%cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act
thal_VL=thal_VL(1,1:4);
thal_VM=thal_VM(1,1:6);
% save 'single_subj_VM_VL.mat'



