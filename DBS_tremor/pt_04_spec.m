% PLS trials: 1st trial spiral at 90o ; 2nd-4th posture at 300o; 5-6th
% spiral at 90o

if cond==1;
    start=start(2);
    ending=ending(2);
    bef_s=floor([149.3]*samplerate);
    tap_in(1,:)=floor([253.7]*samplerate);
    tap_in(2,:)=floor([257]*samplerate);
    tap_out(1,:)=floor([266.2]*samplerate);
    tap_out(2,:)=floor([268.3]*samplerate);
    
else
    start=start([1 3 4]);
    ending=ending([1 3 4]);
    bef_s=floor([4.092 676.9  1067]*samplerate);
    tap_in(1,:)=NaN;
    tap_in(2,:)=NaN;
    tap_out(1,:)=NaN;
    tap_out(2,:)=NaN;
end