
clear
% load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P03_PLS_P')

%  load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P1')
%
 load('/Users/Carolina/OneDrive - Nexus365/DBS-STIM/DATA_pls_sd/P06_PLS_P2')

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;

co=3;
iii=1;
clust=struct; out=struct; start=cell(1,3); ending=cell(1,3); yy=cell(1,3); h_up=cell(1,3); s=struct;

[peak_ax, start, ending, yy, h_up]= dbs_startend(d,samplerateold,samplerate,iii,co,start,ending,yy,h_up);

[s]=dbs_zfiltenv(d,peak_ax,co,iii,s,samplerate);

[b,a]=butter(2,[0.8/(0.5*samplerate) ],'low'); %15
for i=1:3
s.l_filt{iii,co}(i,:)=filtfilt(b,a,s.raw{iii,co}(i,:));
end


for p=1
%     1:length(start{iii,co})
subplot(3,1,1)
plot(d.ds_dall(2,start{iii,co}(p)-5000:ending{iii,co}(p)))
subplot(3,1,2)
plot(s.l_filt{iii,co}(1,start{iii,co}(p)-5000:ending{iii,co}(p)))
subplot(3,1,3)
plot(s.env{iii,co}(1,start{iii,co}(p)-5000:ending{iii,co}(p)))
end

%%% find mhd mhu points
% take mean evelope second before hand up before mhd during handdown second
% after mhup

