function [output] = minstimFxn(inputs)
%for analysis of biophysics differences bw responding and nonresponding neurons in minstim data
%created 05-04-18 %modified 05-04-18

%% INIT VARS
%sample rate
samplerate=inputs.SamplerateHzEditField.Value;
%stimulus start in ms
stimstart = inputs.StimulusstartmsEditField.Value*(samplerate/1000);
%set baseline -501ms to -1ms before stim
blstart=stimstart-inputs.BaselinedurationmsEditField.Value*(samplerate/1000)-1*(samplerate/1000);
blend=blstart+0.5*samplerate;

%% LOAD DATA
cd(inputs.minStimDir);

%stim information
stiminfo = textread('stiminfo.txt'); %#ok<DTXTRD>
output.stimdur = stiminfo(2)/(10^6); %gives stimulation duration in seconds
output.stimIn = stiminfo(3); % gives sitmulation intensity in uA
stimend = stimstart+output.stimdur*samplerate;

%data
contents = dir('*.abf');
filenames = {contents.name}';
minstimFiles = fullfile(cd,filenames);
rawData = abfload(minstimFiles{1},'sweeps','a');
filtData(:,:) = sgolayfilt(rawData(:,1,:),3,11,[],1); %savitsky-golay filter, third order, 11 frame -- creates 3D matrix (MxNxP) where M is time, N is input/command, P is sweep number

%% ANALYSIS
%CHECK RA
%injection for Ra
injstart = inputs.RainjectionstartmsEditField.Value*(samplerate/1000); %voltage injection start
injend = inputs.RainjectionendmsEditField.Value*(samplerate/1000); %voltage injection end

%Ra analysis
[Ra, meanRa, keep] = getRa(rawData,injstart,injend,samplerate) %#ok<*NOPRT>


if keep == 1
    %baseline subtract traces
    blData = -1*ones(size(filtData));
    for ii = 1:size(filtData,2)
        blData(:,ii) = filtData(:,ii) - mean(filtData(blstart:blend,ii));
    end
    
    %get properties of events
    [output.amplitude, output.latency, output.rise, output.decay, output.tAllPeaks] = getPSCData(blData,stimstart,output.stimdur,samplerate,-1);
    
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
        windowCenter = mean(output.tAllPeaks(~isnan(output.tAllPeaks)));
        count = 1;
        for ii = 1:length(output.amplitude)
            if isnan(output.amplitude(ii))
                ampFailures(count)=min(blData(windowCenter-.005*samplerate:windowCenter+.005*samplerate,ii));
                count=count+1;
            end
        end
        
    else
        %success rate
        output.sRate=0;
    end
    
    %stimulus value
    output.stimVal=output.stimIn*output.stimdur*1000; %stimulus value in uA*msec
    
    %% DISPLAY DATA
    %traces
    t=1000.*(1/samplerate:1/samplerate:size(blData,1)/samplerate); %time in milliseconds
    minStimFig.Traces = figure(1);
    minStimFig.Traces.Position=[390 620 450 150];
    hold on
    if (output.sRate ~= 0) & (output.sRate ~= 100) %if mixed successes and failures
        plot(t,blData(:,isnan(output.amplitude)),'color','k')
        plot(t,blData(:,~isnan(output.amplitude)),'color',inputs.plotColor)
        ylim([-1.15*max(output.amplitude) .25*max(output.amplitude)])
    elseif output.sRate == 0 %if no successes
        plot(t,blData,'color','k')
        ylim([-25 25])
    elseif output.sRate == 100 %if only successes
        plot(t,blData,'color',inputs.plotColor)
        ylim([-1.15*max(output.amplitude) .25*max(output.amplitude)])
    end
    xlim((1000/samplerate).*[stimend-.01*samplerate stimend+0.075*samplerate])
    tAx = gca;
    setAx(tAx);
    title('Minimal Stim Responses')
    xlabel('time (ms)')
    ylabel('current (pA)')
    
    if output.sRate ~= 0 %if successes
        %success vs failures
        minStimFig.Vals = figure(2)
        minStimFig.Vals.Position = [390 385 450 150];
        hold on
        scatter(find(~isnan(output.amplitude)),-1.*output.amplitude(~isnan(output.amplitude)),50,inputs.plotColor,'filled')
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