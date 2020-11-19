clear all
cd ('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A2_Thal\mat')
load 'animal_lesion_nolesion.mat'

state_mat=lesion;% lesion vs nolesioned
nfile_i=3:(newfile);
med_thal_cell=cell(1,1);
med_ctx=[];
x=1;
check=[];
thal_all=cell(1,1);
subjects=[];
for nn=1:length(state_mat);
    
    cd C:\Users\wolf5173\Documents\MATLAB\KN
    clearvars -except  subjects x med_thal_cell med_ctx nfile_i state_mat nn A check thal_all
    name=A(state_mat(nn),1:(find(A((state_mat((nn))),:)=='.')-1))
    load(name)
    B=who;
    
    ctxchan=[];
    thalchan=[];
    for ii=1:size(B,1)
        if ~isempty(min(find(B{ii}=='i')) & min(find(B{ii}=='p')) & min(find(B{ii}=='s')))
            ctxchan=ii;
        elseif ~isempty(min(find(B{ii}=='p')) & min(find(B{ii}=='r')) & min(find(B{ii}=='o')) & min(find(B{ii}=='b')))
            thalchan=[thalchan ii];
        end
    end
    
    thal_BZ=[];%basal-ganglia recipient thalamus
    thal_CZ=[];%cerebellar recipient thalamus
    thal_ZI=[];%Zona incerta
    thal_AV=[];%Amteroventral thalamus
    thal_SNr=[];%Substantia Nigra reticulata
    thal_NW=[];%No where
    thal_PC=[];%Para-central thalamus
    thal_AM=[];%Antero medial thalamus
    thal_Rt=[];%Retucular thalamic nucleus
    for i=1:length(thalchan)
        eval(['location=' B{thalchan(i)} '.location']);
        if ~isempty(min(find(location=='B')) & min(find(location=='Z')));
            thal_BZ=[thal_BZ,i];
        elseif ~isempty(min(find(location=='C')) & min(find(location=='Z')));
            thal_CZ=[thal_CZ,i];
        elseif ~isempty(min(find(location=='A')) & min(find(location=='V')));
            thal_AV=[thal_AV,i];
        elseif ~isempty(min(find(location=='Z')) & min(find(location=='I')));
            thal_ZI=[thal_ZI,i];
        elseif ~isempty(min(find(location=='S')) & min(find(location=='N')) & min(find(location=='r')));
            thal_SNr=[thal_SNr,i];
        elseif ~isempty(min(find(location=='n')) & min(find(location=='o')) & min(find(location=='w')));
            thal_NW=[thal_NW,i];
        elseif ~isempty(min(find(location=='P')) & min(find(location=='C')));
            thal_PC=[thal_PC,i];
        elseif ~isempty(min(find(location=='A')) & min(find(location=='M')));
            thal_AM=[thal_AM,i];
        elseif ~isempty(min(find(location=='R')) & min(find(location=='t')));
            thal_Rt=[thal_Rt,i];
        end
    end
    
 
    
    if (length ([thal_BZ thal_CZ thal_ZI thal_AV thal_Rt thal_AM thal_SNr thal_NW thal_PC]))== length (thalchan)
        thal_local=thalchan(thal_BZ); %%% Chose thal location
    else thal_local=NaN;
        
    end
    
    if ~isnan (thal_local);
        subjects=[subjects nn];
    end
        
    if ~isnan (thal_local)
        thal_all{x,1}= thal_local;
        x=x+1;
    end
    
end

clearvars -except thal_all subjects
cd ('C:\Users\wolf5173\Documents\Carolina_code\codes_thal\A1_Thal\mat')
save thal_BZ_all