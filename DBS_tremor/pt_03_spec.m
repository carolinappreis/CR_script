%%%both conditons had phase 0 as a preferred phase

if cond==1;
    start=start(1:3);
    ending=ending(1:3);
    bef_s=floor([56.36 275.8 491.7]*samplerate);
    tap_in(1,:)=floor([112.2 339.7 557.7]*samplerate);
    tap_in(2,:)=floor([143.6 369.4 585.6]*samplerate);
    tap_out(1,:)=floor([116.1 343 560.7]*samplerate);
    tap_out(2,:)=floor([147.4 372 589.1]*samplerate);
    
else
    start=start(4:5);
    ending=ending(4:5);
    bef_s=floor([757.9 1004]*samplerate);
    tap_in(1,:)=NaN;
    tap_in(2,:)=NaN;
    tap_out(1,:)=NaN;
    tap_out(2,:)=NaN;
end