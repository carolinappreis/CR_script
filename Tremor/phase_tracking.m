clear all
iii=[1 2 3 4 5 8 10 11 13];
cor=[];
for numb=1:length(iii);
clearvars -except iii numb psi opsi qnt_psi psi cor psi_w
% load(strcat('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/Random_Stim/RS/P0',num2str(iii(numb)),'_RS.mat'))
load(strcat('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\Random_Stim\RS\P0',num2str(iii(numb)),'_RS.mat'))

data=SmrData.WvData;
samplerateold=SmrData.SR;
tremor=(data(1,:));
addon=92; addon_end=35;

time=0:1/samplerateold:(size(data,2)-1)/samplerateold;

%%% determine stimulation time points
index=[];
for i=2:size(data,2)-1
    if data(2,i-1)<2.5 && data(2,i)>2.5
        index=[index i];
    end
end

clear i
indexes4=[index(1) index(find(diff(index)./samplerateold > 0.95)+1)];
indexes3=[index(find(diff(index)./samplerateold > 0.95)) index(end)];

ph_tremor=angle(hilbert(tremor));

dd2=round(data(4,:)*100)./100;
for l=1:length(indexes4)
    xx(l)=round(dd2(indexes4(l))./0.1); %#ok<*SAGROW>
end


for i=1:12
    idx_1(1:sum(xx==i),i)=indexes4(find(xx==i));
    idx_2(1:sum(xx==i),i)=indexes3(find(xx==i));
end


for k=1:size(idx_1,2)
    run=[];
    for mk=1:size(idx_1,1)
        run=[run index(index>=idx_1(mk,k) & index<=idx_2(mk,k))];
    end
    sav_run{1,k}(:)=run;
    psi(numb,k)=circ_r(ph_tremor(run)'); clear run
end

r=[];
for i=1:12
r=[r length(sav_run{1,i})];
end


if sum(r)==length(index)
    cor(numb)=numb;
else
    cor(numb)=NaN;
end
    
clear run
for n=1:length(indexes4)
    run=index(index>=indexes4(n) & index<=indexes3(n));
    psi_wi{numb,1}(1,n)=circ_r(ph_tremor(run)'); clear run
end

% plot(time,data(2,:))
% hold on 
% plot(time,tremor)
% plot(time,ph_tremor)
% plot(time(run),tremor(run),'r*')
% plot(time(run),ph_tremor(run),'ko')

qnt_psi(numb,1)=numel(find(psi_wi{numb,1}(:)>=0.5));
qnt_psi(numb,2)=(numel(find(psi_wi{numb,1}(:)>=0.5)).*100)./(length(psi_wi{numb,1}));
end


numel(find(psi<0.5))
numel(find(psi<0.8))
