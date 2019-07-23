function [Ra, meanRa, keep] = getRa(data,injstart,injend,samplerate)
% Get Access Resistane Function, created 7-26-2016, modified 08-15-2016
% takes sweeps and outputs access resistance
% also returns the mean Ra
% keep is a logical variable, if max(Ra) <= 25 and does not vary > +/-20% from
% meanRa, keep == 1; if Ra violates these conditions, keep == 0
% assumes that dV is 10mV -- since this sometimes is recorded incorrectly
% all my protocols use a 10mV step for checking Ra. (EMG 08-15-2016)

%init vars
if length(size(data)) == 3
    noSweeps=size(data,3);
    noFiles=1;
elseif length(size(data)) == 4
    noSweeps=size(data,3);
    noFiles=size(data,4);
elseif length(size(data)) < 3
    error('Not enough input arguments in pClamp Data file -- it may be filtered')
end

if noFiles == 1 %calc Ra if only one .abf file
    dV=10.*ones(noSweeps,1);
    dI=zeros(noSweeps,1);
    Ra=zeros(noSweeps,1);
    for ii = 1:noSweeps
        dI(ii) = mean(data((injstart-.051*samplerate):(injstart-.01*samplerate),1,ii)) - min(data(injstart:injend,1,ii)); %difference bw baseline I and Vstep evoked peak I (pA)
        Ra(ii) = 1000*(dV(ii)/dI(ii)); %in MOhms
    end    
elseif noFiles > 1 %calc Ra if multiple .abf files
    dV=10.*ones(noSweeps*noFiles,1);
    dI=zeros(noSweeps*noFiles,1);
    Ra=zeros(noSweeps*noFiles,1);
    for jj = 1:noFiles
        for ii = 1:noSweeps
            dI(ii+jj*noSweeps-noSweeps) = mean(data((injstart-.51*samplerate):(injstart-.01*samplerate),1,ii,jj)) - min(data(injstart:injend,1,ii,jj)); %difference bw baseline I and Vstep evoked peak I (pA)
            Ra(ii+jj*noSweeps-noSweeps) = 1000*(dV(ii+jj*noSweeps-noSweeps)/dI(ii+jj*noSweeps-noSweeps)); %in MOhms
        end
    end
end

meanRa=mean(Ra);
%exclude if it doesn't meet the above criteria
if max(Ra) > 25
    keep = false;
elseif (max(Ra) > meanRa * 1.2) || (min(Ra) < meanRa * .8)
    keep = false;
else
    keep = true;
end

%plot Ra to view
Rafig=figure;
plot(Ra,'-o')
line([0 noSweeps*noFiles+1],[meanRa * 1.2, meanRa * 1.2],'LineStyle','--','Color','k')
line([0 noSweeps*noFiles+1],[meanRa * .8, meanRa * .8],'LineStyle','--','Color','k')
ylim([min(Ra)-5 max(Ra)+5])
xlim([0 noSweeps*noFiles+1])
close(Rafig)

end