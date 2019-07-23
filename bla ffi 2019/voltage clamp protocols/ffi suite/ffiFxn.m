%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FFI analysis For Suite %%%%%%%%%
%%%%%%%%%%% Created: 04-13-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 05-06-2018 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outputdata] = ffiFxn(inputs)
%% INIT VARS
ffiDirec=inputs.ffiDir;
cd(ffiDirec);
samplerate = inputs.SamplerateHzEditField.Value;
stimstart = inputs.StimulusstartmsEditField.Value*samplerate/1000;
stimdur = textread('stimdur.txt'); %#ok<DTXTRD>
stimdur = stimdur/(10^6); %gives stimulation duration in seconds
stimend = stimstart+stimdur*samplerate;
blstart=stimstart-0.5*samplerate-1*(samplerate/1000); %same as getPSCData fxn
blend=blstart+0.5*samplerate;

%% LOAD DATA
%-70 mV
cd([ffiDirec,'/m70']);
contents = dir('*.abf');
filenames = {contents.name}';
ffiFiles.epsc = fullfile(cd,filenames);
rawData.epsc = abfload(ffiFiles.epsc{1},'sweeps','a');
filtData.epsc(:,:) = sgolayfilt(rawData.epsc(:,1,:),3,11,[],1); %savitsky-golay filter, third order, 11 frame (+/- 0.5ms)

%0 mV
cd([ffiDirec,'/zero']);
contents = dir('*.abf');
filenames = {contents.name}';
ffiFiles.ipsc = fullfile(cd,filenames);
rawData.ipsc = abfload(ffiFiles.ipsc{1},'sweeps','a');
filtData.ipsc(:,:) = sgolayfilt(rawData.ipsc(:,1,:),3,11,[],1); %savitsky-golay filter, third order, 11 frame (+/- 0.5ms)
noSweeps = size(filtData.ipsc,2);

%% CHECK RA
%injection for Ra
injstart = inputs.RainjectionstartmsEditField.Value*(samplerate/1000); %voltage injection start
injend = inputs.RainjectionendmsEditField.Value*(samplerate/1000); %voltage injection end

%collate sweeps
allSweeps = cat(3,rawData.epsc,rawData.ipsc);

%Ra analysis
[Ra, meanRa, keep] = getRa(allSweeps,injstart,injend,samplerate) %#ok<*NOPRT>

if keep == 1
    %% DO ANALYSIS
    %baseline subtract traces
    blData.epsc = -1*ones(size(filtData.epsc));
    for ii = 1:size(filtData.epsc,2)
        blData.epsc(:,ii) = filtData.epsc(:,ii) - mean(filtData.epsc(blstart:blend,ii));
    end
    blData.ipsc = -1*ones(size(filtData.ipsc));
    for ii = 1:size(filtData.ipsc,2)
        blData.ipsc(:,ii) = filtData.ipsc(:,ii) - mean(filtData.ipsc(blstart:blend,ii));
    end
    
    %-70 mV
    [amplitude.epsc, latency.epsc, outputdata.rise.epsc, outputdata.decay.epsc, tPeaks.epsc] = getPSCData(blData.epsc,stimstart,stimdur,samplerate,-1);
    
    %0 mV
    [amplitude.ipsc, latency.ipsc, outputdata.rise.ipsc, outputdata.decay.ipsc, tPeaks.ipsc] = getPSCData(blData.ipsc,stimstart,stimdur,samplerate,1);
    
    %get percent failures
    perFail.epsc = 100*(sum(isnan(amplitude.epsc))/length(amplitude.epsc));
    perFail.ipsc = 100*(sum(isnan(amplitude.ipsc))/length(amplitude.ipsc));
    
    %EPSC
    %get baseline subtraced traces
    blcurrent.epsc=zeros(1,noSweeps);
    blData.epsc=zeros(size(filtData.epsc,1),noSweeps);
    for ii = 1:noSweeps
        blcurrent.epsc(ii)=mean(filtData.epsc(blstart:blend,ii));
        blsData.epsc(:,ii)=filtData.epsc(:,ii)-blcurrent.epsc(ii);
    end
    
    %if 100% failures, calculate amplitude as max response at window for PSC detection (5-25ms post stimulus)
    if perFail.epsc == 100
        %mean trace and amplitude
        outputdata.mTrace.epsc = mean(blsData.epsc,2);
        [outputdata.mAmp.epsc, tAmp.epsc] = min(outputdata.mTrace.epsc(stimend+.005*samplerate:stimend+.025*samplerate)); %amplitude
        tAmp.epsc = tAmp.epsc + stimend + .005*samplerate - 1;
        
    else %if any number of successes, calc amp as mean of success amplitudes.
        %amplitude
        outputdata.mAmp.epsc = mean(amplitude.epsc(~isnan(amplitude.epsc)));
        %latency
        outputdata.mLat.epsc = mean(latency.epsc(~isnan(latency.epsc)));
        
        %jitter
        outputdata.jtr.epsc = std(latency.epsc(~isnan(latency.epsc)));
        
        %get mean trace
        [alignedTraces.epsc, outputdata.mTrace.epsc, ~, alignPos.epsc]= meanTraceMaxRise(blsData.epsc(:,~isnan(amplitude.epsc)),tPeaks.epsc(~isnan(amplitude.epsc)),samplerate,.3,-1);
        
    end
    
    %IPSC
    %get baseline subtraced traces
    blcurrent.ipsc=zeros(1,noSweeps);
    blData.ipsc=zeros(size(filtData.ipsc,1),noSweeps);
    for ii = 1:noSweeps
        blcurrent.ipsc(ii)=mean(filtData.ipsc(blstart:blend,ii));
        blsData.ipsc(:,ii)=filtData.ipsc(:,ii)-blcurrent.ipsc(ii);
    end
    %if 100% failures, calculate amplitude as max response at window for PSC detection (5-25ms post stimulus)
    if perFail.ipsc == 100
        %mean trace and amplitude
        outputdata.mTrace.ipsc = mean(blsData.ipsc,2);
        [outputdata.mAmp.ipsc, tAmp.ipsc] = max(outputdata.mTrace.ipsc(stimend+.005*samplerate:stimend+.025*samplerate)); %amplitude
        tAmp.ipsc = tAmp.ipsc + stimend + .005*samplerate - 1;
        
    else %if any number of successes, calc amp as mean of success amplitudes.
        %amplitude
        outputdata.mAmp.ipsc = mean(amplitude.ipsc(~isnan(amplitude.ipsc)));
        %latency
        outputdata.mLat.ipsc = mean(latency.ipsc(~isnan(latency.ipsc)));
        
        %jitter
        outputdata.jtr.ipsc = std(latency.ipsc(~isnan(latency.ipsc)));
        
        %get mean trace
        [alignedTraces.ipsc, outputdata.mTrace.ipsc, ~, alignPos.ipsc]= meanTraceMaxRise(blsData.ipsc(:,~isnan(amplitude.ipsc)),tPeaks.ipsc(~isnan(amplitude.ipsc)),samplerate,.3,1);
        
    end
    %% plot data
    t.epsc = 1/samplerate:1/samplerate:length(outputdata.mTrace.epsc)/samplerate;
    t.epsc = 1000.*t.epsc; %convert to ms
    %epsc
    EPSCfig=figure(1);
    EPSCfig.Position = [420 600 400 175];
    hold on
    if perFail.epsc == 100
        plot(t.epsc,blsData.epsc,'color',[.67 .67 .67],'linewidth',1)
        plot(t.epsc,outputdata.mTrace.epsc,'color',inputs.plotColor,'linewidth',2)
        xlim((1000/samplerate)*[stimend-.015*samplerate stimend+.285*samplerate])
        ylim([-50 50])
    else
        plot(t.epsc,alignedTraces.epsc,'color',[.67 .67 .67],'linewidth',1)
        plot(t.epsc,outputdata.mTrace.epsc,'color',inputs.plotColor,'linewidth',2)
        xlim((1000/samplerate)*[alignPos.epsc-.025*samplerate length(outputdata.mTrace.epsc)])
        ylim([-1.15*max(amplitude.epsc) .15*max(amplitude.epsc)])
    end
    title('EPSC, grey = ind sweeps, color = mean trace')
    xlabel('time (ms)')
    ylabel('pA')
    epscAx = gca;
    epscAx.Box = 'off'; epscAx.YColor = 'k'; epscAx.XColor = 'k'; epscAx.LineWidth = 1; epscAx.TickDir ='out';
    
    t.ipsc = 1/samplerate:1/samplerate:length(outputdata.mTrace.ipsc)/samplerate;
    t.ipsc = 1000.*t.ipsc; %convert to ms
    %ipsc
    ipscfig=figure(2);
    ipscfig.Position = [420 350 400 175];
    hold on
    if perFail.ipsc == 100
        plot(t.ipsc,blsData.ipsc,'color',[.67 .67 .67],'linewidth',1)
        plot(t.ipsc,outputdata.mTrace.ipsc,'color',inputs.plotColor,'linewidth',2)
        xlim((1000/samplerate)*[stimend-.015*samplerate stimend+.285*samplerate])
        ylim([-50 50])
    else
        plot(t.ipsc,alignedTraces.ipsc,'color',[.67 .67 .67],'linewidth',1)
        plot(t.ipsc,outputdata.mTrace.ipsc,'color',inputs.plotColor,'linewidth',2)
        xlim((1000/samplerate)*[alignPos.ipsc-.025*samplerate length(outputdata.mTrace.ipsc)])
        ylim([-.15*max(amplitude.ipsc) 1.15*max(amplitude.ipsc)])
    end
    title('IPSC, grey = ind sweeps, color = mean trace')
    xlabel('time (ms)')
    ylabel('pA')
    ipscAx = gca;
    ipscAx.Box = 'off'; ipscAx.YColor = 'k'; ipscAx.XColor = 'k'; ipscAx.LineWidth = 1; ipscAx.TickDir ='out';

end
end