clear all
%cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')

load('BZ.mat');
% for ii=1:length(BZ.idrat)
% data{ii,1}=BZ.filt_thal{BZ.idrat(ii),1}
% end
% data=vertcat(data{:});
bins=[100:10:300];

BZ.idrat=[1:size(BZ.env_ctx,1)];
for ik=1:length(BZ.idrat)
    ref1=BZ.onset_raw{1,(BZ.idrat(ik))}{2,1};
    ref1_1=BZ.onset_phase_al{1,(BZ.idrat(ik))}{2,1};
    ref2=BZ.offset_raw{1,BZ.idrat(ik)}{2,1};
    if length(ref1) ~= length(ref1_1)
        ref1=ref1(1:length(ref1_1));
        ref2=ref2(1:length(ref1_1));
    end
    %         [dur,dur_idx]=sort(ref2-ref1,'ascend');
    dur=ref2-ref1;
    for i=1:size(bins,2)-1
        for nn=1:size(dur,2);
            if ~isempty (dur(nn)>bins(i) && dur(nn)<=bins(i+1))
                lala(ik,i)=nn;
            end
        end
    end
end

