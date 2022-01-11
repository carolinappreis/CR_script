clear all
close all
cond={'NS';'HF';'C'};
cohort=[ 1 3 4 6];

iii=3;
trial=1;
co=2;

load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/P0',num2str(cohort(iii)),'_clean_',num2str(cond{co,1}),num2str(trial),'_SH.mat'));


t_spi=cell(size(cohort,2),2,size(cond,1));
t_spi{3,1,1}=[[1 2279];[2280 3815];[3816 5742];[5743 7404];[7405 length(tempo)]];
t_spi{3,2,1}=[[1 2102];[2103 3898];[3899 5780];[5781 length(tempo)]];
t_spi{3,1,2}=[[1 1657];[1658 2815];[2816 5040];[5041 6940];[6941 length(tempo)]];
t_spi{3,1,3}=[[1 1936];[1937 3719];[3720 5449];[5450 6968];[6969 8714];[8715 10420];[10421 12130];[12131 13890];[13891 15650];[15651 length(tempo)]];



%%% filt & var
com_sig=sqrt((signal(1,:).^2)+(signal(2,:).^2));
[a,b]=  butter(2, [2/(0.5*samplerate2) 8/(0.5*samplerate2)], 'bandpass'); %15
filt_cs=(filtfilt(a,b,com_sig));
shift=abs(min(filt_cs)-1);
fdata=filt_cs+shift;

env=abs(hilbert(filt_cs));

%%% spatial
% signal=data;
centre=[535 361];
% centre=[0 0];
d1=signal(1,:)-centre(1);
d2=signal(2,:)-centre(2);
[px,py]=cart2pol(d1,d2);


[xx,yy]=cart2pol(centre(1),centre(2));

for mm=1:size(t_spi{iii,trial,co}(:,1),1)
    
    dumpxx=(px(t_spi{iii,trial,co}(mm,1): t_spi{iii,trial,co}(mm,2)))-xx;
    
    pxx=rad2deg(dumpxx);
    
    idx{1,1}=find(pxx>0 & pxx< 90);
    idx{2,1}=find(pxx<0 & pxx>-90);
    idx{3,1}=find(pxx>-180 & pxx<-90);
    idx{4,1}=find(pxx>90 & pxx<180);
    
    
    for i=1:4
        
        ddum=diff(idx{i,1});
        nid=idx{i,1}(find(ddum>1));
        
        figure(1)
        polarplot(px(idx{i,1}),py(idx{i,1}),'.')
        hold on
        figure(2)
        plot(tempo(idx{i,1}),fdata(idx{i,1}),'.')
        hold on
        
        
        mpipi(mm,i) = peak2peak(fdata(idx{i,1}));
        mvar(mm,i) = var(fdata(idx{i,1}));
        menv(mm,i)= mean(env(idx{i,1}));
        
        clear ddum nid
    end
    clear idx pxx dumpxx
end
