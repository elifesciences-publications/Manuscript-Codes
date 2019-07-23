function [amplitude, latency, rise, decay, tPeaks] = getUnitaryData(data,tSpks,samplerate,inputCurrDir)
% Function to get data on post-synaptic currents, created 02-05-2019,
% modified 02-05-19
% takes sweeps and outputs data on the evoked unitary post-synaptic currents
% returns PSC amplitude, latency, 20-80 rise time, weighted decay tau

%%inits
%input inits
if length(size(data)) == 2
    noSweeps=size(data,2);
elseif length(size(data)) > 2
    error('Too many input arguments in data file')
end
%baseline inits -- 500ms before 1ms before stimulus start
blstart=tSpks-0.5*samplerate-1*(samplerate/1000);
blend=blstart+0.5*samplerate;
%create vectors for output
amplitude=zeros(noSweeps,1);
latency=zeros(noSweeps,1);
rawEventStart=zeros(noSweeps,1);
rise=zeros(noSweeps,1);
%create other vectors
peakRawCurrent=zeros(noSweeps,1);
tPeaks=zeros(noSweeps,1);

%% find data
%get data on the currents
for ii = 1:noSweeps
    
    %amplitude, unitary event window up to 3ms after spike
    if inputCurrDir == -1
        [amplitude(ii),peakRawCurrent(ii),tPeaks(ii)]=...
            findPSCpeak(data,tSpks(ii),.003*samplerate,blstart(ii),blend(ii),ii,inputCurrDir);
    elseif inputCurrDir == 1 %if IPSC have a window for 15ms
        if (tPeakMedian-.01*samplerate) > (stimend+artDelay)
            [amplitude(ii),peakRawCurrent(ii),tPeaks(ii)]=...
                findPSCpeak(data,tPeakMedian-.0075*samplerate,.015*samplerate,blstart,blend,jj,ii,inputCurrDir);
        else
            windowShorten = (stimend+artDelay)-(tPeakMedian-.0075*samplerate);
            window = .015*samplerate - windowShorten;
            [amplitude(ii),peakRawCurrent(ii),tPeaks(ii)]=...
                findPSCpeak(data,stimend+artDelay,window,blstart,blend,jj,ii,inputCurrDir);
        end
    end
%     if tPeaks(ii) < stimend+artDelay %if "event" occurred before window start, likely part of stim artefact so exclude
%         tPeaks(ii) = NaN;
%         amplitude(ii) = NaN;
%         peakRawCurrent(ii) = NaN;
%     else
%         amplitude=abs(amplitude); %gets absolute value of current amplitude in pA
%     end
    
    %get rise time and latency
    %rise time: 20-80% rise time
    %latency: time from stimulus end to 20% rise
    if isnan(amplitude(ii)) ~= 1 %if event occurred
        dataForRiseLat = data(tSpks(ii):tPeaks(ii),ii); %get window of data from stim end to peak
        [rise(ii), ~, latency(ii)] = ...
            findRise(dataForRiseLat,1,length(dataForRiseLat),samplerate);
        latency(ii) = latency(ii) - 1; %scales back to end of stimulus
        latency(ii) = latency(ii)*(1000/samplerate); %converts to ms
    elseif isnan(amplitude(ii)) == 1   %if no event occurred
        latency(ii) = NaN;
        rise(ii) = NaN;
    end
    
    %check to verify real event vs part of decay (if non-existant rise but "peak" detected)
    if isnan(amplitude(ii)) ~= 1 %if event occurred
        if isnan(rise(ii)) %if nonexistant rise
            amplitude(ii) = NaN;
            latency(ii) = NaN;
        end
    end
    
    %check that correct peak is found
    viewsweepsFig=figure('Position',[2445 200 900 600]);
    hold on
    plot(data(:,ii))
    if isnan(amplitude(ii)) ~= 1
        scatter(tPeaks(ii),peakRawCurrent(ii))
    end
    line([tSpks(ii)-.002*samplerate tSpks(ii)+.05*samplerate],...
        [mean(data(blstart:blend,ii))-6*getMAD(data(blstart:blend,ii)) ...
        mean(data(blstart:blend,ii))-6*getMAD(data(blstart:blend,ii))],'color','k')
    line([tSpks(ii)-.002*samplerate tSpks(ii)+.05*samplerate],...
        [mean(data(blstart:blend,ii))+6*getMAD(data(blstart:blend,ii)) ...
        mean(data(blstart:blend,ii))+6*getMAD(data(blstart:blend,ii))],'color','k')
    xlim([tSpks(ii)-.002*samplerate tSpks(ii)+.05*samplerate])
    if inputCurrDir == -1
        pfft = 1; %pause here to examine every sweep by eye for EPSCs
    elseif inputCurrDir == 1
        pfft = 1; %pause here to examine every sweep by eye for IPSCs
    end
    close(viewsweepsFig)
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
end

function [amp, peakRawCurrent, time] = findPSCpeak(thisData,tStart,tWindow,blstart,blend,sweepNo,currentDir)

%get shortened vector of data for analysis
madBaseline=getMAD(thisData(blstart:blend,sweepNo)); %median absolute deviation instead of standard deviation as it isn't affected by large deflections from mean as much as sd is
data2analyze=thisData(uint32(tStart):uint32(tStart+tWindow),sweepNo);

%find peaks
[valleys(:,1),valleys(:,2)]=findvalleys(data2analyze);
[peaks(:,1),peaks(:,2)]=findpeaks(data2analyze);

%check that correct peak is found
        viewsweepsFig=figure;
        hold on
        plot(data2analyze)
        
        line([1 length(data2analyze)],[-6*madBaseline -6*madBaseline],'color','k')
        line([1 length(data2analyze)],[6*madBaseline 6*madBaseline],'color','k')

if currentDir == -1 %if inward is larger
    abovethresh=valleys(:,1)<-6*madBaseline; % find those valleys from above that are above thresh
    if sum(abovethresh) > 0
        peaksAboveThresh(:,1)=valleys(abovethresh,1);
        peaksAboveThresh(:,2)=valleys(abovethresh,2);
        [amp,peakInd]=min(peaksAboveThresh(:,1));
        peakRawCurrent=amp;
        time=peaksAboveThresh(peakInd,2)+tStart-1;
        
        %check that correct peak is found
        scatter(peaksAboveThresh(peakInd,2),amp)
    else
        amp = NaN;
        peakRawCurrent=NaN;
        time = NaN;
    end
elseif currentDir == 1 %if outward is larger
    abovethresh=peaks(:,1)>6*madBaseline; % find those peaks from above that are above thresh
    if sum(abovethresh) > 0
        peaksAboveThresh(:,1)=peaks(abovethresh,1);
        peaksAboveThresh(:,2)=peaks(abovethresh,2);
        [amp,peakInd]=max(peaksAboveThresh(:,1));
        peakRawCurrent=amp;
        time=peaksAboveThresh(peakInd,2)+tStart-1;
        
        %check that correct peak is found
        scatter(peaksAboveThresh(peakInd,2),amp)
    else
        amp = NaN;
        peakRawCurrent=NaN;
        time = NaN;
    end
end
close(viewsweepsFig)
end
