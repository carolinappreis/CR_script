[b,a]=butter(2,[0.1/(0.5*samplerate) ],'low'); %15
C=(filtfilt(b,a,envelope));
A=mean(C);
ind_s=[];
for i=11:length(C)-11
    if C(i-1)<A && C(i+1)>A
        ind_s=[ind_s i]; %#ok<*AGROW>
    end
end
for i=1:(length(ind_s)-1)
    if ind_s(i+1)-ind_s(i)==1
        ind_s(i+1)=NaN;
    end
end
ind_s2=ind_s(~isnan(ind_s));
ind_e=[];
for i=11:length(C)-11
    if C(i-1)>A && C(i+1)<A
        ind_e=[ind_e i];
    end
end
for i=1:(length(ind_e)-1)
    if ind_e(i+1)-ind_e(i)==1
        ind_e(i+1)=NaN;
    end
end
ind_e2=ind_e(~isnan(ind_e));
AA=ind_s2(find((ind_e2-ind_s2)>10000));
BB=ind_e2(find((ind_e2-ind_s2)>10000));

close all
plot(C)
hold on
plot(AA,C(AA),'r.')
plot(BB,C(BB),'r.')