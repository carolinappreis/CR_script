function [start, ending, tremor_ds, samplerate]=match(SmrData,iii,co)

[d]=dbs_preprocess(SmrData); samplerateold=d.samplerateold; samplerate=d.samplerate;
data=d.data_raw;  tremor_ds=d.data_ds;
addon=92; addon_end=35;

if co==1
    NS_Spiral
    start=hu{iii,:};
    ending=hd{iii,:};
    
elseif co==2
    HFS_BE_S; start=hu{iii,:}; ending=hd{iii,:};
else
    new = find(data(2,:) > 4);
    difp = find((diff(new)) > 100000); % are you trying to threshold at 9.6 secos?
    ep_1 = [new(difp) new(end)];
    sp_1 = [new(1) new(difp+1)];
    
    dum=find((ep_1-sp_1)<(60000./samplerate)*samplerateold);
    if ~isempty(dum)
        sp_1(1,dum)=NaN;
        ep_1(1,dum)=NaN;
    end
    sp=sp_1(~isnan(sp_1));
    ep=ep_1(~isnan(ep_1));
    
%                 time=1:length(data(1,:));
%                 plot(time,data(4,:))
%                 hold on
%                 plot(time,data(3,:))
%                 plot(time(sp),data(4,sp),'r.')
%                 plot(time(ep),data(4,ep),'k.')
    start= floor((sp./samplerateold)*samplerate)+addon; ending = floor((ep./samplerateold)*samplerate)+addon+addon_end;
  
    if iii==2
        start(1)=79820;
        ending(1)=183800;
       
   end
    
    if iii==4
        ending=[150900];
    end
%     
%     if iii==3
%         ending(1)=109000;
%         start(2)=180700;
%         ending(2)=256000;
%         start(3)=529000;
%         ending(3)=595900;
%     end
end
end