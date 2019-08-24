function output = upscFxn(inputs)
%for comparing unitary currents
%created 02-05-18 %modified 02-05-19

%% INIT VARS
close all

samplerate=inputs.SamplerateHzEditField.Value;
stimCh = 1;
if str2double(inputs.RecordingChannelDropDown.Value) == 1 
    stimCh = 2;
end
noSweeps = inputs.NumberofsweepsEditField.Value;

%% LOAD DATA
cd(inputs.upscDir);

%data
contents = dir('*.abf');
filenames = {contents.name}';
minstimFiles = fullfile(cd,filenames);
rawData = abfload(minstimFiles{1},'sweeps','a');
recordingCellData.raw = rawData(:,(2*str2double(inputs.RecordingChannelDropDown.Value)-1):2*str2double(inputs.RecordingChannelDropDown.Value),:);
recordingCellData.filt(:,:) = sgolayfilt(recordingCellData.raw(:,1,:),3,11,[],1); %savitsky-golay filter, third order, 11 frame -- creates 3D matrix (MxNxP) where M is time, N is input/command, P is sweep number
stimulusCellData = rawData(:,(2*stimCh-1):2*stimCh,:);

%% ANALYSIS
%CHECK RA
%injection for Ra
injstart = inputs.RainjectionstartmsEditField.Value*(samplerate/1000); %voltage injection start
injend = inputs.RainjectionendmsEditField.Value*(samplerate/1000); %voltage injection end

%Ra analysis
[Ra, meanRa, keep] = getRa(recordingCellData.raw,injstart,injend,samplerate) %#ok<*NOPRT>


if keep == 1
    %find Spikes and define unitary event latency as time from peak
    %find by getting max point above 0 mV, report in data point number
    for ii = 1:noSweeps
        tOvershoot = find(stimulusCellData(:,1,ii)>0);
        [~, tSpk(ii)] = max(stimulusCellData(tOvershoot,1,ii));
        tSpk(ii) = tSpk(ii) + tOvershoot(1) - 1;
    end
    
    %baseline subtract traces, set baseline -501ms to -1ms before spike peak
    blData = -1*ones(size(recordingCellData.filt));
    for ii = 1:noSweeps
        blData(:,ii) = recordingCellData.filt(:,ii) - mean(recordingCellData.filt(tSpk(ii)-.501*samplerate:tSpk(ii)-.001*samplerate,ii));
    end
    
    %get properties of events
    [output.amplitude, output.latency, output.rise, output.decay, output.tAllPeaks] = getUnitaryData(blData,tSpk,samplerate,-1);
    
    if sum(isnan(output.amplitude)) ~= length(output.amplitude)
        %mean amplitude
        amp=output.amplitude;
        amp(isnan(amp))=[];
        output.mAmp=mean(amp);
        
        %latency info
        lat=output.latency;
        lat(isnan(lat))=[];
        output.mLat=mean(lat);
        output.jitter=std(lat);
        
        %kinetics
        rt = output.rise;
        rt(isnan(rt)) = [];
        output.mRT = mean(rt);
        output.decay;
        
        %success rate
        noFail=sum(isnan(output.amplitude));
        output.sRate=100*((length(output.amplitude)-noFail)/length(output.amplitude));
        
        %mean trace
        [output.alignedTraces,output.mTrace]=meanTraceMaxRise(blData(:,~isnan(output.amplitude)),output.tAllPeaks(~isnan(output.amplitude)),samplerate,.15,-1);
        
        %get values of failures
        count = 1;
        for ii = 1:length(output.amplitude)
            if isnan(output.amplitude(ii))
                ampFailures(count)=min(blData(tSpk(ii):tSpk(ii)+.003*samplerate,ii));
                count=count+1;
            end
        end
        
    else
        %success rate
        output.sRate=0;
    end
    
    %% DISPLAY DATA
    %traces
    t=1000.*(1/samplerate:1/samplerate:size(blData,1)/samplerate); %time in milliseconds
    unitaryFig.Traces = figure(1);
    unitaryFig.Traces.Position=[390 620 450 450];
    
    subplot(2,1,1)
    hold on
    if (output.sRate ~= 0) & (output.sRate ~= 100) %if mixed successes and failures
        plot(t,blData(:,isnan(output.amplitude)),'color','k','linewidth',1.5)
        plot(t,blData(:,~isnan(output.amplitude)),'color',inputs.plotColor,'linewidth',1.5)
    elseif output.sRate == 0 %if no successes
        plot(t,blData,'color','k')
        ylim([-25 25])
    elseif output.sRate == 100 %if only successes
        plot(t,blData,'color',inputs.plotColor)
    end
    xlim((1000/samplerate).*[mean(tSpk)-.01*samplerate mean(tSpk)+0.075*samplerate])
    tAx = gca;
    setAx(tAx);
    title('Unitary Responses')
    xlabel('time (ms)')
    ylabel('current (pA)')
    output.wholeTrace = blData;
    
    subplot(2,1,2)
    vTrace(:,:) = stimulusCellData(:,1,:);
    output.meanSpkTrace = mean(vTrace,2);
    plot(t,output.meanSpkTrace,'color',[.35 .35 .35],'linewidth',2);
    xlim((1000/samplerate).*[mean(tSpk)-.01*samplerate mean(tSpk)+0.075*samplerate])
    tAx = gca;
    setAx(tAx);
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    
    if output.sRate ~= 0 %if successes
        %success vs failures
        unitaryFig.Vals = figure(2)
        unitaryFig.Vals.Position = [845 625 450 150];
        hold on
        scatter(find(~isnan(output.amplitude)),output.amplitude(~isnan(output.amplitude)),50,inputs.plotColor,'filled')
        if output.sRate ~= 100
            scatter(find(isnan(output.amplitude)),ampFailures,50,'k','filled')
        end
        valAx = gca;
        setAx(valAx);
        title('Successes vs Failures')
        xlabel('sweep number')
        ylabel('peak amplitude (pA)')
    end
    
end

output
end