% clear all
% cd ('/Users/Carolina/Documents/MATLAB/SUA/Juxta SUA:act:mat')
% % cd ('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A2_Thal/code')

sr=1/unite.interval;
srn=1000;
dataold=unite.values';
dataold=full(dataold);
data=zeros(1,100000);
timeold=0:1/sr:(size(dataold,2)-1)/sr;
time=0:1/srn:(size(data,2)-1)/srn;
spk_t=timeold(find(dataold==1));
spk_tround=round(spk_t,3);
n=[];
for i=1:length(spk_t)
    [ d, ix ] = min( abs( time-spk_tround(i) ) );
    n=[n ix];
end
data(n)=1;
data_ones=find(data==1);

% plot(timeold,dataold)
% hold on
% plot(time,data)

cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load 'data_all.mat'
filtrange=betal_mean;
[b,a]=butter(2,[(filtrange-5)/(0.5*srn) (filtrange+5)/(0.5*srn)],'bandpass');
Ecogfiltered=filtfilt(b,a,WaveData_DC);
env=abs(hilbert(Ecogfiltered));
ang=angle(hilbert(Ecogfiltered));
threshold=prctile(env,75);
tt(size(env,1):size(env,2))=threshold;

indexexceed=find(env>threshold);
diffindex=diff(indexexceed);
pnts=find(diffindex>1);
begin=indexexceed(pnts+1);
ending=indexexceed(pnts);
begin2=[indexexceed(1) begin];
ending2=[ending indexexceed(end)];

ind_b=[];
for i=1:(length(begin2))
    if (ending2(i)-begin2(i))>=50
        ind_b=[ind_b i];
    end
end

begin3=begin2(ind_b);
ending3=ending2(ind_b);

% duration=ending3-begin3;
% median_b=median(duration);
% SD_b=std(duration);

thr=150;
ind_b1=[];
for i=1:(length(begin3)-2)
    if (begin3(i+1)-ending3(i))>=thr && (begin3(i+2)-ending3(i+1))>=thr
        ind_b1=[ind_b1 i+1];
    end
end

if (begin3(2)-ending(1))>=thr
    ind_b1= [1 ind_b1];
end

if (begin3(length(begin3))-ending3(length(begin3)-1))>=thr
    ind_b1=[ind_b1 length(begin3)];
end
    
onset1=begin3(ind_b1);
offset1=ending3(ind_b1);

% plot(time,env)
% hold on
% plot(time,tt)
% plot(time(begin2),env(begin2),'r.')
% plot(time(begin3),env(begin3),'y*')
% plot(time(onset1),env(onset1),'bo')

ang_ones=ang(data_ones);
