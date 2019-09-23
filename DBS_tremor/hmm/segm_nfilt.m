
function [epch]=segm_nfilt(tremor2,cc);
if cc==2
    segmentb=[8646 88131 174501 234200 315901 393800];
    segmente=[82701 147501 226501 297200 377500 454600];
    
elseif cc==3
    segmentb=[1621 112001 205601 287501 366701 450401];
    segmente=[76121 179301 263501 350901 422701 517501];
end

epch=[];
for tr=1:length(segmentb)
t=15000;
run=tremor2(segmentb(tr):segmente(tr));
epch=[epch (tremor2(((numel(run)/2)-t):((numel(run)/2)+t)-1))];
end
dum=1;