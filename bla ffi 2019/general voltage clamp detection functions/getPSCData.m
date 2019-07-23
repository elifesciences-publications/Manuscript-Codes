function [amplitude, latency, rise, decay, tPeaks] = getPSCData(data,stimstart,stimdur,samplerate,inputCurrDir)
% Function to get data on post-synaptic currents, created 04-27-2018,
% modified 04-26-18
% takes sweeps and outputs data on the evoked post-synaptic currents
% returns PSC amplitude, latency, 10-90 rise time, weighted decay tau

%%inits 
%input inits
if length(size(data)) == 2
    noSweeps=size(data,2);
    noFiles=1;
elseif length(size(data)) == 3
    noSweeps=size(data,2);
    noFiles=size(data,3);
elseif length(size(data)) > 3
    error('Too many input arguments in data file')
end
%baseline inits -- 500ms before 1ms before stimulus start
blstart=stimstart-0.5*samplerate-1*(samplerate/1000);
blend=blstart+0.5*samplerate;
%stim inits
stimend=stimstart+stimdur*samplerate;
%create vectors for output
amplitude=zeros(noSweeps*noFiles,1);
latency=zeros(noSweeps*noFiles,1);
rawEventStart=zeros(noSweeps*noFiles,1);
rise=zeros(noSweeps*noFiles,1);
%create other vectors
peakRawCurrent=zeros(noSweeps*noFiles,1);
tPeaks=zeros(noSweeps*noFiles,1);

%% find data
%sample random 25% of sweeps to find mean peak time (or 5, if less than 20 sweeps)
if noSweeps*noFiles < 20
    noRSamples=5;
else
    noRSamples=round(0.25*noSweeps*noFiles);
end
sampleSweeps=randperm(noSweeps*noFiles,noRSamples); %sweeps to sample

%pick stimulus artefact window
thisTime = 1000.*(1/samplerate:1/samplerate:size(data,1)/samplerate);
stimArtefactFig = figure;
plot(thisTime,data)
line([stimend*(1000/samplerate) stimend*(1000/samplerate)],[-50 50],'color','k','linewidth',2)
xlim([stimstart*(1000/samplerate)-5 stimend*(1000/samplerate)+10]) %15ms window around stimulus
artDelay = inputdlg('stim artefact buffer time (ms) from end of stimulus to start window?');
artDelay = (samplerate/1000)*str2double(artDelay);
close(stimArtefactFig)

%get times 
tCurrent=zeros(1,length(noRSamples));
currDir=zeros(1,length(noRSamples));
for ii = 1:noRSamples
    if rem(sampleSweeps(ii),noSweeps) == 0
        idFile=sampleSweeps(ii)/noSweeps;
        idSweep=noSweeps;
    else
        idFile=ceil(sampleSweeps(ii)/noSweeps);
        idSweep=rem(sampleSweeps(ii),noSweeps);
    end
    if inputCurrDir == -1 %20ms window for epsc
        [~,~,tCurrent(ii),currDir(ii)]=findPSCpeak(data,stimend+artDelay,.02*samplerate,blstart,blend,idFile,idSweep,inputCurrDir);
    elseif inputCurrDir == 1 %20ms window for ipsc
        [~,~,tCurrent(ii),currDir(ii)]=findPSCpeak(data,stimend+artDelay,.02*samplerate,blstart,blend,idFile,idSweep,inputCurrDir);
    end
end
tCurrent(isnan(tCurrent))=[];

if isempty(tCurrent) ~=1
    tPeakMedian=uint32(median(tCurrent)); %median time of event to center window for EPSC detection
    
    %get data on the currents
    for jj = 1:noFiles
        for ii = 1:noSweeps
            
            %amplitude
            if inputCurrDir == -1 %if EPSC have a window for 10ms
                if (tPeakMedian-.005*samplerate) > (stimend+artDelay)
                    [amplitude(ii+jj*noSweeps-noSweeps),peakRawCurrent(ii+jj*noSweeps-noSweeps),tPeaks(ii+jj*noSweeps-noSweeps)]=...
                        findPSCpeak(data,tPeakMedian-.005*samplerate,.01*samplerate,blstart,blend,jj,ii,inputCurrDir);
                else
                    windowShorten = (stimend+artDelay)-(tPeakMedian-.005*samplerate);
                    window = .01*samplerate - windowShorten;
                    [amplitude(ii+jj*noSweeps-noSweeps),peakRawCurrent(ii+jj*noSweeps-noSweeps),tPeaks(ii+jj*noSweeps-noSweeps)]=...
                        findPSCpeak(data,stimend+artDelay,window,blstart,blend,jj,ii,inputCurrDir);
                end
            elseif inputCurrDir == 1 %if IPSC have a window for 15ms
                if (tPeakMedian-.01*samplerate) > (stimend+artDelay)
                    [amplitude(ii+jj*noSweeps-noSweeps),peakRawCurrent(ii+jj*noSweeps-noSweeps),tPeaks(ii+jj*noSweeps-noSweeps)]=...
                        findPSCpeak(data,tPeakMedian-.0075*samplerate,.015*samplerate,blstart,blend,jj,ii,inputCurrDir);
                else
                    windowShorten = (stimend+artDelay)-(tPeakMedian-.0075*samplerate);
                    window = .015*samplerate - windowShorten;
                    [amplitude(ii+jj*noSweeps-noSweeps),peakRawCurrent(ii+jj*noSweeps-noSweeps),tPeaks(ii+jj*noSweeps-noSweeps)]=...
                        findPSCpeak(data,stimend+artDelay,window,blstart,blend,jj,ii,inputCurrDir);
                end
            end
            if tPeaks(ii+jj*noSweeps-noSweeps) < stimend+artDelay %if "event" occurred before window start, likely part of stim artefact so exclude
                tPeaks(ii+jj*noSweeps-noSweeps) = NaN;
                amplitude(ii+jj*noSweeps-noSweeps) = NaN;
                peakRawCurrent(ii+jj*noSweeps-noSweeps) = NaN;
            else
                amplitude=abs(amplitude); %gets absolute value of current amplitude in pA
            end
            
            %get rise time and latency
            %rise time: 20-80% rise time
            %latency: time from stimulus end to 20% rise
            if isnan(amplitude(ii+jj*noSweeps-noSweeps)) ~= 1 %if event occurred
                if noFiles == 1
                    dataForRiseLat = data(stimend+artDelay:tPeaks(ii+jj*noSweeps-noSweeps),ii); %get window of data from stim end to peak
                    dataForRiseLat = dataForRiseLat - mean(data(blstart:blend,ii)); %baseline subtract data
                elseif noFiles > 1
                    dataForRiseLat = data(stimend+artDelay:tPeaks(ii+jj*noSweeps-noSweeps),ii,jj); %get window of data from stim end to peak
                    dataForRiseLat = dataForRiseLat - mean(data(blstart:blend,ii,jj)); %baseline subtract data
                end
                [rise(ii+jj*noSweeps-noSweeps), ~, latency(ii+jj*noSweeps-noSweeps)] = ...
                    findRise(dataForRiseLat,1,length(dataForRiseLat),samplerate);
                latency(ii+jj*noSweeps-noSweeps) = latency(ii+jj*noSweeps-noSweeps) +artDelay - 1; %scales back to end of stimulus
                latency(ii+jj*noSweeps-noSweeps) = latency(ii+jj*noSweeps-noSweeps)*(1000/samplerate); %converts to ms
            elseif isnan(amplitude(ii+jj*noSweeps-noSweeps)) == 1   %if no event occurred
                latency(ii+jj*noSweeps-noSweeps) = NaN;
                rise(ii+jj*noSweeps-noSweeps) = NaN;
            end
            

            %check that correct peak is found
            viewsweepsFig=figure('Position',[2445 200 900 600]);
            hold on
            plot(data(:,ii,jj))
            scatter(tPeaks(ii+jj*noSweeps-noSweeps),peakRawCurrent(ii+jj*noSweeps-noSweeps))
            if isnan(latency(ii+jj*noSweeps-noSweeps)) ~= 1
                scatter(uint32(latency(ii+jj*noSweeps-noSweeps)*(samplerate/1000)+stimend-1),data(uint32(latency(ii+jj*noSweeps-noSweeps)*(samplerate/1000)+stimend-1),ii,jj))
            end
            if noFiles > 1
                line([stimstart-.002*samplerate stimend+.05*samplerate],...
                    [mean(data(blstart:blend,ii,jj))-6*getMAD(data(blstart:blend,ii,jj)) ...
                    mean(data(blstart:blend,ii,jj))-6*getMAD(data(blstart:blend,ii,jj))],'color','k')
                line([stimstart-.002*samplerate stimend+.05*samplerate],...
                    [mean(data(blstart:blend,ii,jj))+6*getMAD(data(blstart:blend,ii,jj)) ...
                    mean(data(blstart:blend,ii,jj))+6*getMAD(data(blstart:blend,ii,jj))],'color','k')
            elseif noFiles == 1
                line([stimstart-.002*samplerate stimend+.05*samplerate],...
                    [mean(data(blstart:blend,ii))-6*getMAD(data(blstart:blend,ii)) ...
                    mean(data(blstart:blend,ii))-6*getMAD(data(blstart:blend,ii))],'color','k')
                line([stimstart-.002*samplerate stimend+.05*samplerate],...
                    [mean(data(blstart:blend,ii))+6*getMAD(data(blstart:blend,ii)) ...
                    mean(data(blstart:blend,ii))+6*getMAD(data(blstart:blend,ii))],'color','k')
            end
            xlim([stimstart-.002*samplerate stimend+.05*samplerate])
            if isnan(amplitude(ii+jj*noSweeps-noSweeps)) ~= 1 %if event occurred
                ylim([-1.25*amplitude(ii+jj*noSweeps-noSweeps) 1.25*amplitude(ii+jj*noSweeps-noSweeps)])
            else
                ylim([-10*getMAD(data(blstart:blend,ii)) 10*getMAD(data(blstart:blend,ii))])
            end
            if inputCurrDir == -1
                pfft = 1; %pause here to examine every sweep by eye for EPSCs 
            elseif inputCurrDir == 1
                pfft = 1; %pause here to examine every sweep by eye for IPSCs 
            end
            close(viewsweepsFig)
        end
    end
    
    %get mean trace aligned max rise
    if sum(~isnan(amplitude)) > 0 %if any successes get mean data and tau
        if inputCurrDir == -1
            [~, meanTrace, ~, tMaxSlope] = meanTraceMaxRise(data(:,~isnan(amplitude)),tPeaks(~isnan(amplitude)),samplerate,.125,inputCurrDir); 
        elseif inputCurrDir == 1
            [~, meanTrace, ~, tMaxSlope] = meanTraceMaxRise(data(:,~isnan(amplitude)),tPeaks(~isnan(amplitude)),samplerate,.25,inputCurrDir); 
        end
        
        %weighted decay tau
        decay=findTau(meanTrace,tMaxSlope,samplerate);
    else
        decay = NaN;
    end
    
elseif isempty(tCurrent)
    amplitude=NaN;
    latency=NaN;
    rise=NaN;
    decay=NaN;
    tAllPeaks=NaN;
end
end

function [amp, peakRawCurrent, time, direction] = findPSCpeak(thisData,tStart,tWindow,blstart,blend,fileNo,sweepNo,currentDir)

%get shortened vector of data for analysis
if length(size(thisData)) == 2
    blCurrent=mean(thisData(blstart:blend,sweepNo));
    madBaseline=getMAD(thisData(blstart:blend,sweepNo)); %median absolute deviation instead of standard deviation as it isn't affected by large deflections from mean as much as sd is
    data2analyze=thisData(uint32(tStart):uint32(tStart+tWindow),sweepNo)-blCurrent;
elseif length(size(thisData)) == 3
    blCurrent=mean(thisData(blstart:blend,sweepNo,fileNo));
    madBaseline=getMAD(thisData(blstart:blend,sweepNo,fileNo));
    data2analyze=thisData(uint32(tStart):uint32(tStart+tWindow),sweepNo,fileNo)-blCurrent;
end

%find peaks
[valleys(:,1),valleys(:,2)]=findvalleys(data2analyze);
[peaks(:,1),peaks(:,2)]=findpeaks(data2analyze);

if nargin == 7 %if direction of current is not given
    if isempty(peaks) %if inward is larger
        direction=-1;
    elseif isempty(valleys) %if outward is larger
        direction =1;
    elseif max(abs(valleys(:,1))) >= max(abs(peaks(:,1)))%if inward is larger
        direction=-1;
    elseif max(abs(valleys(:,1))) < max(abs(peaks(:,1))) || isempty(valleys) %if outward is larger
        direction=1;
    end
elseif nargin == 8 %if direction of current is given
    direction = currentDir;
end

if direction == -1 %if inward is larger
    abovethresh=valleys(:,1)<-6*madBaseline; % find those valleys from above that are above thresh
    if sum(abovethresh) > 0
        peaksAboveThresh(:,1)=valleys(abovethresh,1);
        peaksAboveThresh(:,2)=valleys(abovethresh,2);
        [amp,peakInd]=min(peaksAboveThresh(:,1));
        peakRawCurrent=amp+blCurrent;
        time=peaksAboveThresh(peakInd,2)+tStart-1;
        
        %check that correct peak is found
        viewsweepsFig=figure;
        hold on
        plot(data2analyze)
        scatter(peaksAboveThresh(peakInd,2),amp)
        line([1 length(data2analyze)],[-6*madBaseline -6*madBaseline],'color','k')
        line([1 length(data2analyze)],[6*madBaseline 6*madBaseline],'color','k')
        close(viewsweepsFig)
    else
        amp = NaN;
        peakRawCurrent=NaN;
        time = NaN;
    end
elseif direction == 1 %if outward is larger
    abovethresh=peaks(:,1)>6*madBaseline; % find those peaks from above that are above thresh
    if sum(abovethresh) > 0
        peaksAboveThresh(:,1)=peaks(abovethresh,1);
        peaksAboveThresh(:,2)=peaks(abovethresh,2);
        [amp,peakInd]=max(peaksAboveThresh(:,1));
        peakRawCurrent=amp+blCurrent;
        time=peaksAboveThresh(peakInd,2)+tStart-1;
        
        %check that correct peak is found
        viewsweepsFig=figure;
        hold on
        plot(data2analyze)
        scatter(peaksAboveThresh(peakInd,2),amp)
        line([1 length(data2analyze)],[-6*madBaseline -6*madBaseline],'color','k')
        line([1 length(data2analyze)],[6*madBaseline 6*madBaseline],'color','k')
        close(viewsweepsFig)
    else
        amp = NaN;
        peakRawCurrent=NaN;
        time = NaN;
    end
end
end
