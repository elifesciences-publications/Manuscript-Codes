%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% stim. induced spiking for Suite %%%%%%
%%%%%%%%%%% Created: 05-11-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 05-12-2018 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = stimSpkFxn(inputs)
%based on ic_excite.m code
%% INIT VARS
stimSpkDirec=inputs.stimSpkDir;
cd(stimSpkDirec);

%sample rate
samplerate=inputs.SamplerateHzEditField.Value;
%stimulus start
output.stimstart = -1.*ones(1,inputs.NostimuliEditField.Value);
for ii = 1:inputs.NostimuliEditField.Value
    if ii == 1
        output.stimstart(ii) = inputs.StimulusstartmsEditField.Value*(samplerate/1000);
    else
        output.stimstart(ii) = output.stimstart(ii-1) + samplerate*(1/inputs.StimulusFreqHzEditField.Value);
    end
end
%set baseline right before (1ms) first stimulus start
blstart=output.stimstart(1)-inputs.BaselinedurationmsEditField.Value*(samplerate/1000)-1*(samplerate/1000);
blend=blstart+inputs.BaselinedurationmsEditField.Value*(samplerate/1000);

%% LOAD DATA
stimdur = textread('stimdur.txt'); %#ok<DTXTRD>
stimdur = stimdur/(10^6); %gives stimulation duration in seconds
stimend = output.stimstart+stimdur*samplerate;

%data
contents = dir('*.abf');
filenames = {contents.name}';
stimSpkFiles = fullfile(cd,filenames);
rawData = abfload(stimSpkFiles{1},'sweeps','a');
output.vData(:,:) = rawData(:,1,:); %only uses voltage data, reformats to time x trial
noSweeps = size(rawData,3);

%% ANALYSIS
%baseline membrane voltage
output.baseVm=zeros(1,noSweeps);
for ii = 1:noSweeps
    output.baseVm(ii)=mean(output.vData(blstart:blend,ii));
end
output.blVm=mean(output.baseVm);


%spiking
spkPeaks=-1.*ones(noSweeps,inputs.NostimuliEditField.Value,5);
spkTimes=-1.*ones(noSweeps,inputs.NostimuliEditField.Value,5);
spksPerStim=-1.*ones(noSweeps,inputs.NostimuliEditField.Value);
rasterSpks=-1.*ones(noSweeps,inputs.NostimuliEditField.Value,5);
rasterTimes=-1.*ones(noSweeps,inputs.NostimuliEditField.Value,5);
stCount=0;
output.stTrace=[];
output.stEPSP=[];
for ii = 1:noSweeps
    for stim = 1:inputs.NostimuliEditField.Value
        if stim < inputs.NostimuliEditField.Value
            [peaks,locs,~,proms] = findpeaks(output.vData(stimend(stim):output.stimstart(stim+1),ii));
        elseif stim == inputs.NostimuliEditField.Value
            [peaks,locs,~,proms] = findpeaks(output.vData(stimend(stim):output.stimstart(stim)+samplerate*(1/inputs.StimulusFreqHzEditField.Value),ii));
        end
        if sum(proms>35) == 0 %if no APs
            peaks=[]; locs=[]; proms=[];
            stCount=stCount+1;
            [stEPSPval,~]=max(output.vData(stimend(stim)+.0025*samplerate:output.stimstart(stim)+(1/inputs.StimulusFreqHzEditField.Value-.01)*samplerate,ii));
            output.stEPSP(stCount)=stEPSPval-output.baseVm(ii);
            output.stTrace(:,stCount)=output.vData(stimend(stim)+.0025*samplerate:output.stimstart(stim)+(1/inputs.StimulusFreqHzEditField.Value-.01)*samplerate,ii)-output.baseVm(ii);
            stEPSPval=[];
        else %if APs
            peaks(proms<35)=[]; 
            if isempty(peaks) ~= 1
                spkPeaks(ii,stim,1:length(peaks))=peaks;
                rasterSpks(ii,stim,1:length(peaks))=peaks; peaks=[];
            end
            locs(proms<35)=[]; 
            if isempty(locs) ~= 1
                spkTimes(ii,stim,1:length(locs))=locs; 
                rasterTimes(ii,stim,1:length(locs))=locs; 
                spkTimes(ii,stim,1:length(locs))=spkTimes(ii,stim,1:length(locs))+stimend(stim)-1; %correct time point for AP peak within trace
                rasterTimes(ii,stim,1:length(locs))=rasterTimes(ii,stim,1:length(locs))+stimend(stim)-1; locs=[];
            end
            proms=[];
        end
        spksPerStim(ii,stim) = sum(spkTimes(ii,stim,:)~=-1);
    end
    verifyAP=figure(1);
    verifyAP.Position=[400 350 850 375];
    hold on
    plot(output.vData(:,ii))
    for stim = 1:inputs.NostimuliEditField.Value
        scatter(spkTimes(ii,stim,spkTimes(ii,stim,:)~=-1),spkPeaks(ii,stim,spkPeaks(ii,stim,:)~=-1))
    end
    xlim([output.stimstart(1)-.1*samplerate output.stimstart(end)+.15*samplerate])
    close(verifyAP)
end

%remove subthresh EPSPs that are <0 and thus IPSPs
output.stTrace(:,output.stEPSP<0)=[];
output.stEPSP(output.stEPSP<0)=[];

%get means
output.mSpkPerStim.byStim = mean(spksPerStim);
output.mSpkPerStim.bySweep = mean(spksPerStim(:));
output.mEPSPtrace = mean(output.stTrace,2);

%% DISPLAY DATA

%time vector
t.vData = 1000.*(1/samplerate:1/samplerate:size(output.vData,1)/samplerate); %time in ms
t.stEPSP = 1000.*(1/samplerate:1/samplerate:size(output.stTrace,1)/samplerate); %time in ms

%lighter color
lightColor = (1+(.85-max(inputs.plotColor))./max(inputs.plotColor)).*inputs.plotColor; %linearly scales RGB values to a point where max value is .85 for one of R,G,B

%plot raster
rasterFig=figure(1);
rasterFig.Position=[400 575 535 175];
hold on
for ii = 1:noSweeps
    for stim = 1:inputs.NostimuliEditField.Value
        if noSweeps == noSweeps
            line((1000/samplerate).*[output.stimstart(stim) output.stimstart(stim)],...
                [noSweeps+.5 noSweeps+1.5],'linewidth',4,'color','k')
        end
        for spks = 1:5
            if rasterSpks(ii,stim,spks) ~= -1
                line((1000/samplerate).*[rasterTimes(ii,stim,spks) rasterTimes(ii,stim,spks)],...
                    [ii-.5 ii+.5],'linewidth',2,'color',inputs.plotColor)
            end
        end
    end
end
xlabel('time (ms)')
xlim((1000/samplerate).*[output.stimstart(1)-.05*samplerate output.stimstart(end)+.1*samplerate])
ylabel('trial number')
ylim([0 noSweeps+2.5])
title('spiking raster')
rasterAx = gca;
setAx(rasterAx);

%plot traces
tracesFig=figure(2);
tracesFig.Position=[400 100 535 400];
hold on
for ii = 1:noSweeps
    plot(t.vData,output.vData(:,ii)-output.baseVm(ii)+(ii-1)*round(100/noSweeps),'color',inputs.plotColor,'linewidth',1.5)
end
xlim((1000/samplerate).*[output.stimstart(1)-.05*samplerate output.stimstart(end)+.1*samplerate])
ylabel('Vm (mV // trials separated)')
xlabel('time (ms)')
tracesAx = gca;
setAx(tracesAx);
ylim([.5*tracesAx.YLim(1) .95*tracesAx.YLim(2)])
if tracesAx.YLim(2) < 110
    tracesAx.YLim(2) = 110;
end

stEPSPfig=figure(3);
stEPSPfig.Position=[945 475 225 275];
hold on
plot(t.stEPSP,output.stTrace,'linewidth',.5,'color',lightColor)
plot(t.stEPSP,output.mEPSPtrace,'linewidth',2,'color',inputs.plotColor)
stAx = gca;
setAx(stAx);
ylabel('Vm (mV)')
ylim([-3 12])
xlabel('time (ms)')
if ~isempty(output.mEPSPtrace)
    xlim((1000/samplerate).*[0 length(output.mEPSPtrace)])
end
title('subthreshold EPSPs')
legend('ind EPSPs','mean EPSP')
legend('boxoff')

spkPerStim = figure(4);
spkPerStim.Position = [945 115 275 275];
hold on
plot(output.mSpkPerStim.byStim,'linewidth',3,'color',inputs.plotColor)
scatter(inputs.NostimuliEditField.Value+1,output.mSpkPerStim.bySweep,125,inputs.plotColor,'filled')
spsAx = gca;
setAx(spsAx);
xlim([.5 inputs.NostimuliEditField.Value+1.5])
xlabel('Stimulus Number')
spsAx.XTick = 1:1:inputs.NostimuliEditField.Value+1;
spsAx.XTickLabel{end} = 'mean';
if spsAx.YLim(2) < 1
    spsAx.YLim(2) = 1;
end
spsAx.YLim(1) = 0;
ylabel('Spikes per stimulus')


%% output the data into command window
output %#ok<*NOPRT>
end
