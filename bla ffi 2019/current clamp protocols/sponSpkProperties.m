%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Biophys Properties from sAPs %%%%%%%%
%%%%%%%%%%% Created: 06-13-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 06-15-2018 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INIT VARS
close all
clearvars

filt60 = 1; %remove 60Hz noise
window = 100; %in ms

%% Load files
spondir=uigetdir;
cd(spondir);

%data
contents = dir('*.abf');
filenames = {contents.name}';
sponFiles = fullfile(cd,filenames);

%load
[data, sampleInt] = abfload(sponFiles{1},'a','sweeps');
samplerate = (sampleInt/1e6)^-1;
if filt60 == 1
    for ii = 1:size(data,2)
        filtData(:,ii) = templatefilterm(data(:,ii),samplerate/60, 200/samplerate,20); %median based 60Hz removal filter
    end
else
    for ii = 1:size(data,2)
        filtData(:,ii) = data(:,ii);
    end
end
t = (1:1:length(data))/samplerate; %time in seconds
window = window*(samplerate/1000); %converts to dp

%% detect spikes
%use rolling baseline + MAD threshold

%find spontaneous events
for jj = 1:size(filtData,2)
    noSpks(jj) = 0;
    bufferEnd = 0;
    for ii = .1*samplerate+1:length(filtData(:,jj))-samplerate %get from first window/2 ms to window/2 ms before end of trace (want only whole currents)
        if ii > bufferEnd
            noiseMAD = getMAD(filtData(ii-0.025*samplerate:ii-1,jj)); %25ms rolling baseline
            spikeThresh = 10*noiseMAD+median(filtData(ii-0.025*samplerate:ii-1,jj));
            if filtData(ii,jj) > spikeThresh
                [~,loc] = max(filtData(ii-.5*window:ii+.5*window,jj)); %give 150ms buffer for spc isi
                if loc ~= 1
                    if loc ~= 0.05*samplerate+1
                        noSpks(jj) = noSpks(jj) +1
                        tEvents(noSpks(jj),jj) = ii+loc-.5*window-1;
                        spkTrace(:,noSpks(jj),jj) = filtData(tEvents(noSpks(jj),jj)-0.025*samplerate:tEvents(noSpks(jj),jj)+samplerate,jj);
                        bufferEnd = ii + loc + 0.001*samplerate;
                    end
                end
            end
        end
    end
end
% spkTrace = spkTrace';
%baseline subtract
for jj = 1:size(filtData,2)
    for ii = 1:noSpks(jj)
        baseV(ii,jj) = mean(spkTrace(1:.005*samplerate,ii,jj));
        blSpk(:,ii,jj) = spkTrace(:,ii,jj) - baseV(ii,jj);
    end
end


%% measure halfwidth and remove noise (do diff for mouse cells)
%inits
hw=zeros(noSpks(jj),2);
halfamp1=-1.*ones(noSpks(jj),2);
halfamp2=zeros(noSpks(jj),2);
tHalf1=zeros(noSpks(jj),2);
tHalf2=zeros(noSpks(jj),2);

%get amplitude
for jj = 1:size(filtData,2)
    for ii = 1:noSpks(jj)
        peakAmp(ii,jj) = blSpk(.01*samplerate+1,ii,jj);
    end

%get halfwidth
for ii = 1:noSpks(jj)
    dHalfamp1=[];
    dHalfamp2=[]; %#ok<*NASGU>
    halfamp=peakAmp(ii,jj)/2;
    dHalfamp1=abs(blSpk(1:.025*samplerate,ii)-halfamp); %closest to half amp from t = -10ms pre spk to 1dp pre peak
    [~,tminDist1]=min(dHalfamp1);
    dHalfamp2=abs(blSpk(.025*samplerate+2:.075*samplerate,ii)-halfamp); %closest to half amp from t = 1dp post peak to 50ms post peak
    [~,tminDist2]=min(dHalfamp2);
    tHalf1(ii,jj)=tminDist1;
    halfamp1(ii,jj)=blSpk(tHalf1(ii,jj),ii,jj);
    tHalf2(ii,jj)=.01*samplerate+2+tminDist2-1;
    halfamp2(ii,jj)=blSpk(tHalf2(ii,jj),ii,jj);
    hw(ii,jj)=(tHalf2(ii,jj)-tHalf1(ii,jj))/(samplerate/1000);
end

%remove noise (aka spks with hw less than 2ms)
deleteCells = [];
deleteCells = find(hw(:,jj)<2);
if ~isempty(deleteCells)
    noSpks(jj) = noSpks(jj) - length(deleteCells(jj))
    spkTrace(:,deleteCells(jj),jj) = [];
    blSpk(:,deleteCells(jj),jj) = [];
    peakAmp(deleteCells(jj),jj) = [];
    hw(deleteCells(jj),jj) = [];
end
end