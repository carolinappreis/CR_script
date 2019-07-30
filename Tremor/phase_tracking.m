clear all
iii=[1 2 3 4 5 6 8 10 11];
for numb= 1:length(iii);
clearvars -except iii numb psi opsi qnt_psi
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

for n=1:length(indexes4)
    run=index(index>=indexes4(n) & index<=indexes3(n));
    psi{numb,1}(1,n)=circ_r(ph_tremor(run)'); clear run
end

% plot(time,data(2,:))
% hold on 
% plot(time,tremor)
% plot(time,ph_tremor)
% plot(time(run),tremor(run),'r*')
% plot(time(run),ph_tremor(run),'ko')

qnt_psi(numb,1)=numel(find(psi{numb,1}(:)>=0.5));
qnt_psi(numb,2)=(numel(find(psi{numb,1}(:)>=0.5)).*100)./(length(psi{numb,1}));
end



