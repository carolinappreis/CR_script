cd('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA')
cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')
clear; close
cohort = [ 1 3 4 6];
cond={'NS';'RS'};
clust=struct; out=struct; start=cell(10,1); ending=cell(10,1); yy=cell(10,3); h_up=cell(10,3); s=struct; freq_bl=[]; amp_bl=[];cf=struct;
rng('default')
gen=(rng);
spiral=0; %%% posture=0; spiral=1;
opt_out=1; %%% with cluster=0; without cluster =1;

% % if spiral==0
% %    load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_posture_unclust.mat');
% % else
% %     load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/DBS_final_spiral_unclust.mat');
% % end


for iii = 1:length(cohort)
    clearvars -except  opt_out cohort cond iii clust s start ending yy out gen h_up spiral freq_bl amp_bl avg power dat cf Freq_peak

    for co=1:size(cond,1)

        load(strcat('/Users/Carolina/Desktop/oxford/data_code_thesis/DBS/paper_analyses/paper_matfiles/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'));
        %       load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA/P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))
        [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;




        data=d.data_raw;  tremor_ds=d.ds_dall(1,:);

        decide=0; %0= freq of filtering is task specific else is coming from posture



        if co==1

            if decide==0
                if spiral==0
                    NS_BE_P
                else
                    NS_BE_S
                end
            else
                NS_BE_P
            end

            start{iii,co}=hu{iii,:};
            ending{iii,co}=hd{iii,:};


            % when patient's hand is up
            handup = [];
            for i = 1:length(start{iii,co})
                handup = [handup start{iii,co}(i):ending{iii,co}(i)]; %#ok<*AGROW>
            end
            clear i
            handup = sort(handup,'ascend');

            h_up{iii,co}=handup;

            [Pxx,F] = pwelch(tremor_ds(handup), samplerate, [], 2*samplerate, samplerate);
            Ppeak1=max(Pxx(5:100));
            Ppeak=F(find(Pxx==Ppeak1));
            clear Pxx F
        else

            %----------------------------------------------------------------------------------------------

            ts=[];
            te=[]
            addon=92; addon_end=35;

            % determine stimulation time points
            index = [];
            for i = 2:size(data, 2)-1
                if data(2,i-1)<2.5 && data(2,i)>2.5
                    index = [index i];
                end
            end
            clear i

            % Find trigger
            indexes4 = [index(1) index(find(diff(index) ./ samplerateold > 0.95)+1)];
            indexes3 = [index(find(diff(index) ./ samplerateold > 0.95)) index(end)];


            dd2 = round(data(4, :) *100) ./ 100;
            for i = 1:length(indexes4)
                xx(i) = round(dd2(indexes4(i)) ./ 0.1); %#ok<*SAGROW>
            end
            clear i


            ss = floor((indexes4 ./ samplerateold)*samplerate)+addon;
            ee = floor((indexes3 ./ samplerateold)*samplerate)+addon+addon_end;%floor(5*samplerate);

            % when patient's hand is up
            handup = [];
            for i = 1:length(ss)
                handup = [handup ss(i):ee(i)]; %#ok<*AGROW>
            end
            clear i
            handup = sort(handup,'ascend');
            h_up{iii,co}=handup;


            [Pxx,F] = pwelch(tremor_ds(handup), samplerate, [], 2*samplerate, samplerate);
             Ppeak1=max(Pxx(5:100));
            Ppeak=F(find(Pxx==Ppeak1));
            clear Pxx F


        end
        Freq_peak(iii,co)=Ppeak;
        
    end
end



