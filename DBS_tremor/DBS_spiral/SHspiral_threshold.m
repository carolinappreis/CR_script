clear
load('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/P01_NS1_SH.mat')
bins=find(NS1(:,1)==-1);
for y=1:length(bins)
    if(y+1)<length(bins)
        data1{1,y}=NS1(bins(y)+1:bins(y+1)-1,:);
    end
end
data1{1,length(bins)-1}=NS1(bins(y)+1:end,:);


NS1(bins,:)=[];

% find(NS1(:,1)<10) %%% to find turn new whole spiral

%%% Data by spiral
%%% for nr=1:size(data1,2)
nr=1;
data=data1{1,nr};

%%%OR all data
% data=NS1;

samplerate=floor(1000/median(diff(data(:,3))));
[a,b]=  butter(2, [2/(0.5*samplerate) 6/(0.5*samplerate)], 'bandpass'); %15

% dum=factor(length(data));
% dum1=abs(dum-samplerate);
% dum2=(find(dum1==min(dum1)));
% L=dum(dum2);
% %     L=13;
% n = floor(length(data) / L); % Number of segments

L=samplerate;
n = floor(length(data) / L); % Number of segments

for o=1:2
    filt(:,o)=filtfilt(a,b,data(:,o));
    idx_st = NaN(n, 1);
    idx_st(1, 1) = 1;
    idx_end = NaN(n, 1);
    idx_end(1, 1) = L;
    f_new = NaN(n, L);
    for i = 2:n
        idx_st(i) = idx_st(i-1) + L;
        idx_end(i) = idx_end(i-1) + L;
    end
    for i = 1:n
        f_new(i, :) = filt( idx_st(i, 1) : idx_end(i, 1),o);
    end
    resamp{1,o}=f_new; clear f_new
end



new_f=[];
for j = 1:size(resamp{1,1},1)
    x = [resamp{1,1}(j,:); resamp{1,2}(j,:)];
    [pc, score, latent, tsquare, explained] = pca(x');
    comb_pc(j, 1:2) = pc(1:2, 1);
    maxx=find(pc(:, 1)==max(pc(:,1)));
    new_f=[new_f resamp{1,maxx}(j,:)];
    clear x
end

figure()
plot(comb_pc(:,1), 'r.')
hold on
plot(comb_pc(:,2), 'b.')
xlim([0,size(comb_pc, 1)])

new_e=abs(hilbert(new_f));

th=prctile(new_e,75);
id_sp1=find(new_e>th);
id_sp2=find(new_e<th);

data_s=data(1:length(new_e),:);
time=1:length(data_s);

figure
plot(time,new_f)
hold on 
plot(time,new_e)
yline(th)
plot(time(id_sp1),new_e(id_sp1),'r.')
plot(time(id_sp2),new_e(id_sp2),'b.')


%     figure
%     plot3(time,data_s(:,1),data_s(:,2),'Color',[0.5 0.5 0.5])
%     hold on
%     plot3(time(id_sp1),data_s(id_sp1,1),data_s(id_sp1,2),'r.','MarkerSize',10)
%     plot3(time(id_sp2),data_s(id_sp2,1),data_s(id_sp2,2),'b.','MarkerSize',10)

figure
plot(data_s(:,1),data_s(:,2),'Color',[0.5 0.5 0.5])
hold on
plot(data_s(id_sp1,1),data_s(id_sp1,2),'r.','MarkerSize',10)
plot(data_s(id_sp2,1),data_s(id_sp2,2),'b.','MarkerSize',10)


close all

