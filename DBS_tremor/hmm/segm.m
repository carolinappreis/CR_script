
function [epch]=segm(tremor_or,cc,samplerate)
if cc==1
    segmentb=round(([17320].*samplerate)./1000,0);
    segmente=round(([82510].*samplerate)./1000,0);
    
elseif cc==2
    segmentb=round(([8646 88131 174501 234200 315901 393800].*samplerate)./1000,0);
    segmente=round(([82701 147501 226501 297200 377500 454600].*samplerate)./1000,0);
    
elseif cc==3
    segmentb=round(([1621 112001 205601 287501 366701 450401].*samplerate)./1000,0);
    segmente=round(([76121 179301 263501 350901 422701 517501].*samplerate)./1000,0);
end


t=(15000.*samplerate)./1000;
if cc==1
    bin=round((segmente-segmentb)./(2*t),0);
    epch=tremor_or(segmentb:segmentb+((bin*(2*t))-1)); 
else
    epch=[];
    for tr=1:length(segmentb)
        run=round((length(segmentb(tr):segmente(tr)))./2,0);    
        epch=[epch tremor_or((segmentb(tr)+run-t):(segmentb(tr)+run+t-1))];
    end
    
end
