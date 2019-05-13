clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
%  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all' , 'freq')

%  cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/SUA/Juxta SUA_act_mat/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load 'animal_region.mat' % comes from code pre_filegen_SUA_act
% load 'non_repeat_animal_region.mat'
% load 'single_subj_VM_VL.mat'
% all_regions={thal_VL; thal_VM};
thal_CZ= thal_VL;
thal_BZ=[thal_VM thal_VA];
% all_regions={thal_CZ; thal_BZ;}
all_regions={thal_BZ;};
n=[17.5:7.5:32.5];
freq=cell(1,1);
for i=1:length(n)
    freq{1,i}=[n(i)-5 n(i)+5];
end

for t=1:size(freq,2)
    for i=1:size(all_regions,1)
        for  j=1:size(all_regions{i,:},2)
            
            clearvars -except A all_regions i j ang ang_cons regional_ang m data_all euler euler1 ang_mean ang_cons t freq ang1 vect_length ray_test ecog_units
            name=A(all_regions{i,:}(j),1:(find(A((all_regions{i,:}(j)),:)=='.')-1));
            load(name);
            B=who;
            ctxchan=[];
            for ii=1:size(B,1)
                if ~isempty(min(find(B{ii}=='E')) & min(find(B{ii}=='E')) & min(find(B{ii}=='G')))
                    ctxchan=ii;
                end
            end
            
            
            eval(['samprateold=1/' B{ctxchan} '.interval;']);
            eval(['WaveData(1,:)=' B{ctxchan} '.values;']);
            WaveData=double(WaveData);
            ts=timeseries(WaveData(1,:),0:(1/samprateold):((size(WaveData,2)-1)/samprateold));
            ts1=resample(ts,0:0.001:((size(WaveData,2)-1)/samprateold),'linear');
            WaveData_DC(1,:)=ts1.data;
            
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
            for ii=1:length(spk_t)
                [ d, ix ] = min( abs( time-spk_tround(ii) ) );
                nn=[nn ix];
            end
            data(nn)=1;
            data_all{i,j}=data;
            data_ones=find(data==1);
            [b,a]=butter(2,[freq{1,t}(1)/(0.5*srn) freq{1,t}(end)/(0.5*srn)],'bandpass');
            Ecogfiltered=filtfilt(b,a,WaveData_DC);
            ecog_units(j,t,:)=Ecogfiltered;
            ang{1,j}(:,t)=wrapTo2Pi(angle(hilbert(Ecogfiltered(data_ones))));
            vect_length(j,t)=circ_r(ang{1,j}(:,t));
            ray_test(j,t)=circ_rtest(ang{1,j}(:,t));
        end
    end
end

r=ray_test<0.05;
stat_d=[];
for i=1:10
    if sum(r(i,:))>0
        stat_d=[stat_d;i];
    end
end

for i =1:length(stat_d)
    if find(vect_length(stat_d(i),1:2)==max(vect_length(stat_d(i),1:2)))==1
        ecogbf_match(i,:)=ecog_units(i,1,:);
    else find(vect_length(stat_d(i),1:2)==max(vect_length(stat_d(i),1:2)))==2
        ecogbf_match(i,:)=ecog_units(i,2,:);
    end
end


units_match=cell2mat(data_all(1,stat_d(1:end))');

clearvars -except units_match ecogbf_match stat_d
save 'BZ_cycle'


