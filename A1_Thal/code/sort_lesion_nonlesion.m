clear all
cd ('\Users\creis\Documents\MATLAB\KN')
A=ls;
newfile=[];
nolesion=[];
lesion=[];
newfile=size(A,1);
for nfile_i=3:(newfile)
    name=A(nfile_i,1:(find(A(nfile_i,:)=='.')-1));
    load(name)
    if IpsiEEG.dopamine=='6-OHDA'
        lesion=[lesion nfile_i];

    else nolesion=[nolesion nfile_i];

    end
end

clearvars -except nolesion lesion A newfile
% save 'animal_lesion_nolesion'
