
function [epch]=segm(envelope,cc);
if cc==1
    segmentb=[17320];
    segmente=[82510];
    
elseif cc==2
    segmentb=[8646 88131 174501 234200 315901 393800];
    segmente=[82701 147501 226501 297200 377500 454600];
    
elseif cc==3
    segmentb=[1621 112001 205601 287501 366701 450401];
    segmente=[76121 179301 263501 350901 422701 517501];
end


t=15000;
if cc==1
    epch=envelope(1:2*t);
    bin=1:(2*t):(segmente-segmentb);
    for i=1:length(bin)
        if i+1<length(bin)
        epch=[epch envelope( bin(i):bin(i+1)-1)];
        end
    end
    
    
else
    epch=[];
    for tr=1:length(segmentb)
        run=envelope(segmentb(tr):segmente(tr));
        epch=[epch (envelope(((numel(run)/2)-t):((numel(run)/2)+t)-1))];
    end
    
end
