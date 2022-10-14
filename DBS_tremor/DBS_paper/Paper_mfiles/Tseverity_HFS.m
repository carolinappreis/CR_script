clear all; close all
cohort = [1 3 4 6];
cond={'HFS'};
folder={'/DATA_hf/'};
spiral=0;
% start=cell(10,1); ending=cell(10,1); yy=cell(10,3); s=struct;
yy=cell(10,3); s=struct;
HFS_BE_P;
cohort = [ 1 3 4 6];
spiral=0; [match_ax]=link_ax(spiral);


for iii=1:length(cohort)
    clearvars -except hd hu iii folder cohort cond yy s spiral match_ax sev

    co=1;
    load(strcat('/Users/Carolina/OneDrive - Nexus365/DBS-STIM',num2str(folder{co,1}),'P0',num2str(cohort(iii)),'_',num2str(cond{co,1}),'.mat'))

    [d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;

    start=hu{iii,:}; ending=hd{iii,:};

    handup = [];
    for i = 1:size(start,2)
        handup = [handup start(i):ending(i)];
    end
    h_up = sort(handup,'ascend');

    tremor_ds=d.data_ds;
    for aa = 1:3
        [Pxx,F] = pwelch(tremor_ds(aa,handup), samplerate, [], samplerate, samplerate);
        frange = F(3:10);
        Pxxrange = Pxx(3:10);
        Freqpeak(aa,:) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(aa,:) = max(Pxxrange);
        ps_curves(aa,:) = Pxx;

        clear Pxx F
    end

    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    [s]=zfiltenv_simple(d,peak_ax,co,iii,s,samplerate);


    sev(iii,:)=[mean(s.env_acc{iii,1}(match_ax(2,iii,1),h_up)) std(s.env_acc{iii,1}(match_ax(2,iii,1),h_up))];
end
