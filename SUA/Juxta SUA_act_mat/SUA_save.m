clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
%cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all' , 'freq')

% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act

all_regions={[thal_VL];[thal_VA thal_VM]};
for ii=1
%     :size(all_regions,1)
    for  j=1:size(all_regions{ii,:},2)
        
        clearvars -except data_all data_ones WaveData_DC rec_pa rec_npa A all_regions ii j region_pl region_npl region_spl region_snpl output_pa output_npa ISI_all spikerate_all
        name=A(all_regions{ii,:}(j),1:(find(A((all_regions{ii,:}(j)),:)=='.')-1));
        cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
        %         cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
        
        load(name);
        B=who;
        ctxchan=[];
        for rr=1:size(B,1)
            if ~isempty(min(find(B{rr}=='E')) & min(find(B{rr}=='E')) & min(find(B{rr}=='G')))
                ctxchan=rr;
            end
        end
        
        %--------------------------------------------------------------
        eval(['samprateold=1/' B{ctxchan} '.interval;']);
        eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
        WaveData=double(WaveData);
        ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
        ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
        WaveData_DC(j,:)=ts1.data;
        
        %--------------------------------------------------------------
        sr=1/unite.interval;
        srn=1000;
        dataold=unite.values';
        dataold=full(dataold);
        data=zeros(1,100000);
        timeold=0:1/sr:(size(dataold,2)-1)/sr;
        time=0:1/srn:(size(data,2)-1)/srn;
        spk_t=timeold(find(dataold==1));
        spk_tround=round(spk_t,3);
        nn=[];
        for i=1:length(spk_t)
            [ d, ix ] = min( abs( time-spk_tround(i) ) );
            nn=[nn ix];
        end
        data(nn)=1;
        data_all(j,:)=data;
        data_ones{j,:}=data(data==1);
    end
end

 clearvars -except data_ones data_all WaveData_DC time srn
 save 'SUA_CZ'