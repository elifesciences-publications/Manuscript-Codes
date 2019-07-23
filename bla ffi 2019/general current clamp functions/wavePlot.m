%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% interneuron waveform plot %%%%%%
%%%%%%%%% Created: 04-24-2017 %%%%%%%%%%
%%%%%%%%%% Edited: 04-23-2019 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INIT VARS
close all
clearvars

%% LOAD DATA
%PV cells
pvdir=uigetdir;
cd(pvdir);
contents = dir('*.mat');
filenames = {contents.name}';
bp_files.pv = fullfile(cd,filenames);

%SST cells
sstdir=uigetdir;
cd(sstdir);
cd([sstdir,'/FS-Sst']);
contents = dir('*.mat');
filenames = {contents.name}';
bp_files.fs = fullfile(cd,filenames);
cd([sstdir,'/nFS-Sst']);
contents = dir('*.mat');
filenames = {contents.name}';
bp_files.nfs = fullfile(cd,filenames);
%pv
for ii = 1:length(bp_files.pv)
    load(bp_files.pv{ii},'app')
    if size(app.ciData.rWaves,1) > 1
        wForm.pv(ii,:) = mean(app.ciData.rWaves); clear app;
    else
        wForm.pv(ii,:) = app.ciData.rWaves; clear app;
    end
end

%fs
for ii = 1:length(bp_files.fs)
    load(bp_files.fs{ii},'app')
    if size(app.ciData.rWaves,1) > 1
        wForm.fs(ii,:) = mean(app.ciData.rWaves); clear app;
    else
        wForm.fs(ii,:) = app.ciData.rWaves; clear app;
    end
end

%nfs
for ii = 1:length(bp_files.nfs)
    load(bp_files.nfs{ii},'app')
    if size(app.ciData.rWaves,1) > 1
        wForm.nfs(ii,:) = mean(app.ciData.rWaves); clear app;
    else
        wForm.nfs(ii,:) = app.ciData.rWaves; clear app;
    end
end


%% Get Data
%voltage -- x axis
meanV.pv = mean(wForm.pv(:,2:end));
semV.pv = sem(wForm.pv(:,2:end),1);
meanV.fs = mean(wForm.fs(:,2:end));
semV.fs = sem(wForm.fs(:,2:end),1);
meanV.nfs = mean(wForm.nfs(:,2:end));
semV.nfs = sem(wForm.nfs(:,2:end),1);

%dV/dt
for ii = 1:size(wForm.pv,1)
    dVdt.pv(ii,:)=diff(wForm.pv(ii,:)).*10; %gives dVdt in mV/ms
end
for ii = 1:size(wForm.fs,1)
    dVdt.fs(ii,:)=diff(wForm.fs(ii,:)).*10; %gives dVdt in mV/ms
end
for ii = 1:size(wForm.nfs,1)
    dVdt.nfs(ii,:)=diff(wForm.nfs(ii,:)).*10; %gives dVdt in mV/ms
end
meandVdt.pv = mean(dVdt.pv);
semdVdt.pv = sem(dVdt.pv,1);
meandVdt.fs = mean(dVdt.fs);
semdVdt.fs = sem(dVdt.fs,1);
meandVdt.nfs = mean(dVdt.nfs);
semdVdt.nfs = sem(dVdt.nfs,1);

%% Plot
%colors
green = [.35 .6 .25]; %green
orange = [.85 .35 .1]; %orange
purple = [.6 .4 .67]; %purple
pink = [.9 .2 .4]; %pink

wPlot = figure(1);
wPlot.Position = [1035 375 350 500];
hold on
%pv
plot(meanV.pv,meandVdt.pv,'color',green)
for ii = 1:length(meanV.pv)
    line([meanV.pv(ii)-semV.pv(ii) meanV.pv(ii)+semV.pv(ii)],[meandVdt.pv(ii) meandVdt.pv(ii)],'color',green)
    line([meanV.pv(ii) meanV.pv(ii)],[meandVdt.pv(ii)-semdVdt.pv(ii) meandVdt.pv(ii)+semdVdt.pv(ii)],'color',green)
end
xlim([-65 25])
ylim([-125 205])
wAx = gca;
wAx.Box = 'off'; wAx.LineWidth = 1; wAx.XColor = 'k'; wAx.YColor = 'k';
wAx.TickDir = 'out'; wAx.XAxisLocation = 'origin'; wAx.YAxisLocation = 'origin';
wAx.XTick = [-60 -40 -20 20]; wAx.YTick = [-100 -50 50 100 150 200];

w2Plot = figure(2);
w2Plot.Position = [1035 375 350 500];
hold on
%fs
plot(meanV.fs,meandVdt.fs,'color',purple)
for ii = 1:length(meanV.fs)
    line([meanV.fs(ii)-semV.fs(ii) meanV.fs(ii)+semV.fs(ii)],[meandVdt.fs(ii) meandVdt.fs(ii)],'color',purple)
    line([meanV.fs(ii) meanV.fs(ii)],[meandVdt.fs(ii)-semdVdt.fs(ii) meandVdt.fs(ii)+semdVdt.fs(ii)],'color',purple)
end
xlim([-65 25])
ylim([-125 205])
wAx = gca;
wAx.Box = 'off'; wAx.LineWidth = 1; wAx.XColor = 'k'; wAx.YColor = 'k';
wAx.TickDir = 'out'; wAx.XAxisLocation = 'origin'; wAx.YAxisLocation = 'origin';
wAx.XTick = [-60 -40 -20 20]; wAx.YTick = [-100 -50 50 100 150 200];

w3Plot = figure(3);
w3Plot.Position = [1035 375 350 500];
hold on
%nfs
plot(meanV.nfs,meandVdt.nfs,'color',pink)
for ii = 1:length(meanV.nfs)
    line([meanV.nfs(ii)-semV.nfs(ii) meanV.nfs(ii)+semV.nfs(ii)],[meandVdt.nfs(ii) meandVdt.nfs(ii)],'color',pink)
    line([meanV.nfs(ii) meanV.nfs(ii)],[meandVdt.nfs(ii)-semdVdt.nfs(ii) meandVdt.nfs(ii)+semdVdt.nfs(ii)],'color',pink)
end
xlim([-65 25])
ylim([-125 205])
wAx = gca;
wAx.Box = 'off'; wAx.LineWidth = 1; wAx.XColor = 'k'; wAx.YColor = 'k';
wAx.TickDir = 'out'; wAx.XAxisLocation = 'origin'; wAx.YAxisLocation = 'origin';
wAx.XTick = [-60 -40 -20 20]; wAx.YTick = [-100 -50 50 100 150 200];
