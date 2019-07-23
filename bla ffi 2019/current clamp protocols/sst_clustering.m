%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% sst in clustering %%%%%%%%%%%
%%%%%%%%% Created: 06-10-2017 %%%%%%%%%%
%%%%%%%%%% Edited: 05-12-2019 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INIT VARS
close all
clearvars

noSst = 0;

%% LOAD DATA

%SST cells
sstdir=uigetdir;
cd(sstdir);
%data
contents = dir('*.mat');
filenames = {contents.name}';
sstFiles = fullfile(cd,filenames);

%load
sstVars = -1.*ones(length(sstFiles),15);
halfwidths = -1.*ones(length(sstFiles),1);
for ii = 1:length(sstFiles)
    
    %find correlations to be used
    load(sstFiles{ii})
    sstVars(ii,1)=app.ciData.rm;
    sstVars(ii,2)=app.ciData.maxFR;
    maxFR(ii) = app.ciData.maxFR;
    sstVars(ii,3)=app.ciData.dAHP;
    sstVars(ii,4)=app.ciData.vsag; 
    sags(ii) = app.ciData.vsag;
    sstVars(ii,5)=mean(app.ciData.halfwidth);
    halfwidths(ii) = mean(app.ciData.halfwidth);
    sstVars(ii,6)=mean(app.ciData.rheoAHPval); 
    sstVars(ii,7)=mean(app.ciData.rheoAHPlat);
    ahpLat(ii) = mean(app.ciData.rheoAHPlat);
    sstVars(ii,8)=mean(app.ciData.rheoThreshold); 
    sstVars(ii,9)=app.ciData.frRatio; 
    sstVars(ii,10)=app.ciData.ampRatio;
    sstVars(ii,11)=app.ciData.spkLat; 
    sstVars(ii,12)=app.ciData.broadeningRatio; 
    sstVars(ii,13)=mean((app.ciData.valRheoPeaks-app.ciData.rheoThreshold')); %spike amplitude; top left to bottom right diagonal of the rheoAmp matrix
    sstVars(ii,14)=app.ciData.noRheoSpks; 
    sstVars(ii,15)=app.ciData.mtau;

    noSst = noSst +1;
    clear app
end

%% FIND PARAMS TO USE -- only worry about if new data added.
count = 0;
for ii = 1:size(sstVars,2)
    for jj = 1:(size(sstVars,2)-ii)
        count = count +1;
        [rho(count),p(count)]=corr(sstVars(:,ii),sstVars(:,ii+jj));
        comp(count,:) = [ii,ii+jj];
    end
end

alpha = drsFDRpval(p);
sigCorrs = find(p<alpha);
theseCorrs = comp(sigCorrs,:);

%% DO CLUSTERS
%zscore transformation
sstzScores = -1.*ones(length(sstFiles),15);
for ii = 1:size(sstVars,2) %loop through variables
    sstzScores(:,ii) = zscore(sstVars(:,ii));
end
sstClus = linkage(sstzScores,'ward','euclidean');

%Scree Plot
screeFig = figure(1);
screeFig.Position = [1200 650 250 150];
hold on
for ii = 1:10 %for first 15 centroids
    scatter(ii,sstClus(end-(ii-1),3),50,'k','filled')
end
sAx = gca;
sAx.Box = 'off'; sAx.YColor= 'k'; sAx.XColor= 'k'; sAx.LineWidth = 1; sAx.TickDir = 'out';

%Dendrogram
dFig = figure(2);
dFig.Position = [45 375 1150 420];
[H,T,outperm] = dendrogram(sstClus,0);
for ii = 1:length(H)
    H(ii).LineWidth = 2; H(ii).Color = 'k';
end
dAx = gca;
dAx.Box = 'off'; dAx.YColor= 'k'; dAx.LineWidth = 1; dAx.TickDir = 'out';

%% PLOT HISTOGRAMS
%plotting histograms of decision tree branches
%colors
grpOne = [.5 .33 .55];
grpTwo = [247 168 184]/255;

%% first cut is max firing rate
%histogram edges.maxFR
edges.maxFR = 0:5:200;

%get data by cluster, (cell#95 is the last cell in group 1)
clusterOneNeurons = outperm(1:find(outperm==95));
clusterTwoNeurons = outperm(find(outperm==95)+1:end);
clusterOneFR = maxFR(clusterOneNeurons);
clusterTwoFR = maxFR(clusterTwoNeurons);
clusterOneCounts = histc(clusterOneFR,edges.maxFR);
clusterTwoCounts = histc(clusterTwoFR,edges.maxFR);
sumCounts = clusterOneCounts + clusterTwoCounts;

%plot
histAllFig = figure(3);
histAllFig.Position = [50 150 300 135];
hold on
bar(edges.maxFR,sumCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.maxFR,clusterTwoCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.maxFR(1) edges.maxFR(end)])
allAx = gca;
setAx(allAx);
xlabel('maximum firing rate (Hz)')

%< 95Hz only
%get data by cluster
clusterOneLessThanFR= clusterOneFR(find(clusterOneFR<95));
clusterTwoLessThanFR= clusterTwoFR(find(clusterTwoFR<95));
clusterOneLessThanFRInd = find(clusterOneFR<95); 
clusterTwoLessThanFRInd = find(clusterTwoFR<95);
clusterOneLessThanCounts = histc(clusterOneLessThanFR,edges.maxFR);
clusterTwoLessThanCounts = histc(clusterTwoLessThanFR,edges.maxFR);
sumLessThanCounts = clusterOneLessThanCounts + clusterTwoLessThanCounts;

%plot
histLessThanFig = figure(4);
histLessThanFig.Position = [360 150 300 135];
hold on
bar(edges.maxFR,sumLessThanCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.maxFR,clusterTwoLessThanCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.maxFR(1) edges.maxFR(end)])
lessAx = gca;
setAx(lessAx);
xlabel('maximum firing rate (hz)')

%>= 95Hz only
%get data by cluster
clusterOneGreaterThanFR = clusterOneFR(find(clusterOneFR>=95));
clusterTwoGreaterThanFR = clusterTwoFR(find(clusterTwoFR>=95));
clusterOneGreaterThanFRInd = find(clusterOneFR>=95); 
clusterTwoGreaterThanFRInd = find(clusterTwoFR>=95);
clusterOneGreaterThanCounts = histc(clusterOneGreaterThanFR,edges.maxFR);
clusterTwoGreaterThanCounts = histc(clusterTwoGreaterThanFR,edges.maxFR);
sumGreaterThanCounts = clusterOneGreaterThanCounts + clusterTwoGreaterThanCounts;

%plot
histGreaterThanFig = figure(5);
histGreaterThanFig.Position = [670 150 300 135];
hold on
bar(edges.maxFR,sumGreaterThanCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.maxFR,clusterTwoGreaterThanCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.maxFR(1) edges.maxFR(end)])
greaterAx = gca;
setAx(greaterAx);
xlabel('maximum firing rate (hz)')

%% second cut is sag
%histogram edges.sag
edges.sag = 0:1:40;

%get data by cluster, only with Hz >= 95
clusterOneFirstCut = clusterOneNeurons(clusterOneGreaterThanFRInd);
clusterTwoFirstCut = clusterTwoNeurons(clusterTwoGreaterThanFRInd);
clusterOneGreaterFRsags = sags(clusterOneFirstCut);
clusterTwoGreaterFRsags = sags(clusterTwoFirstCut);
clusterOneFirstCutGreaterCounts = histc(clusterOneGreaterFRsags,edges.sag);
clusterTwoFirstCutGreaterCounts = histc(clusterTwoGreaterFRsags,edges.sag);
sumCountsFirstCutGreater = clusterOneFirstCutGreaterCounts + clusterTwoFirstCutGreaterCounts;

% %plot
histAllFig = figure(6);
histAllFig.Position = [50 150 300 135];
hold on
bar(edges.sag,sumCountsFirstCutGreater,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
bar(edges.sag,clusterOneFirstCutGreaterCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
xlim([edges.sag(1) edges.sag(end)])
allAx = gca;
setAx(allAx);
xlabel('sag (%)')

%< 20% only
%get data by cluster
clusterOneLessThanSagGreaterFR= clusterOneGreaterFRsags(find(clusterOneGreaterFRsags<20));
clusterTwoLessThanSagGreaterFR= clusterTwoGreaterFRsags(find(clusterTwoGreaterFRsags<20));
clusterOneLessThanSagGreaterFRCounts = histc(clusterOneLessThanSagGreaterFR,edges.sag);
clusterTwoLessThanSagGreaterFRCounts = histc(clusterTwoLessThanSagGreaterFR,edges.sag);
sumLessThanSagGreaterFRCounts = clusterOneLessThanSagGreaterFRCounts + clusterTwoLessThanSagGreaterFRCounts;

%plot
histLessThanFig = figure(7);
histLessThanFig.Position = [360 150 300 135];
hold on
bar(edges.sag,sumLessThanSagGreaterFRCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.sag,clusterTwoLessThanSagGreaterFRCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.sag(1) edges.sag(end)])
lessAx = gca;
setAx(lessAx);
xlabel('sag (%)')

%>= 20% only
%get data by cluster
clusterOneGreaterThanSagGreaterFR = clusterOneGreaterFRsags(find(clusterOneGreaterFRsags>=20));
clusterTwoGreaterThanSagGreaterFR = clusterTwoGreaterFRsags(find(clusterTwoGreaterFRsags>=20));
clusterOneGreaterThanSagGreaterFRCounts = histc(clusterOneGreaterThanSagGreaterFR,edges.sag);
clusterTwoGreaterThanSagGreaterFRCounts = histc(clusterTwoGreaterThanSagGreaterFR,edges.sag);
sumGreaterThanSagGreaterFRCounts = clusterOneGreaterThanSagGreaterFRCounts + clusterTwoGreaterThanSagGreaterFRCounts;

%plot
histGreaterThanFig = figure(8);
histGreaterThanFig.Position = [670 150 300 135];
hold on
bar(edges.sag,sumGreaterThanSagGreaterFRCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.sag,clusterTwoGreaterThanSagGreaterFRCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.sag(1) edges.sag(end)])
greaterAx = gca;
setAx(greaterAx);
xlabel('sag (%)')

%get data by cluster, only with Hz < 95
clusterOneFirstCutLessThan = clusterOneNeurons(clusterOneLessThanFRInd);
clusterTwoFirstCutLessThan = clusterTwoNeurons(clusterTwoLessThanFRInd);
clusterOneLessFRsags = sags(clusterOneFirstCutLessThan);
clusterTwoLessFRsags = sags(clusterTwoFirstCutLessThan);
clusterOneFirstCutLessCounts = histc(clusterOneLessFRsags,edges.sag);
clusterTwoFirstCutLessCounts = histc(clusterTwoLessFRsags,edges.sag);
sumCountsFirstCutLess = clusterOneFirstCutLessCounts + clusterTwoFirstCutLessCounts;

% %plot
histAllFig = figure(13);
histAllFig.Position = [50 150 300 135];
hold on
bar(edges.sag,sumCountsFirstCutLess,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
bar(edges.sag,clusterOneFirstCutLessCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
xlim([edges.sag(1) edges.sag(end)])
allAx = gca;
setAx(allAx);
xlabel('sag (%)')

%< 5.9% only
%get data by cluster
clusterOneLessThanSagLessFR= clusterOneLessFRsags(find(clusterOneLessFRsags<5.9));
clusterTwoLessThanSagLessFR= clusterTwoLessFRsags(find(clusterTwoLessFRsags<5.9));
clusterOneLessThanSagLessFRIndHold = find(clusterOneLessFRsags<5.9);
clusterOneLessThanSagLessFRInd = clusterOneLessThanFRInd(clusterOneLessThanSagLessFRIndHold);
clusterTwoLessThanSagLessFRIndHold = find(clusterTwoLessFRsags<5.9);
clusterTwoLessThanSagLessFRInd = clusterTwoLessThanFRInd(clusterTwoLessThanSagLessFRIndHold);
clusterOneLessThanSagLessFRCounts = histc(clusterOneLessThanSagLessFR,edges.sag);
clusterTwoLessThanSagLessFRCounts = histc(clusterTwoLessThanSagLessFR,edges.sag);
sumLessThanSagLessFRCounts = clusterOneLessThanSagLessFRCounts + clusterTwoLessThanSagLessFRCounts;

%plot
histLessThanFig = figure(14);
histLessThanFig.Position = [360 150 300 135];
hold on
bar(edges.sag,sumLessThanSagLessFRCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.sag,clusterTwoLessThanSagLessFRCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.sag(1) edges.sag(end)])
lessAx = gca;
setAx(lessAx);
xlabel('sag (%)')

%>= 5.9% only
%get data by cluster
clusterOneGreaterThanSagLessFR = clusterOneLessFRsags(find(clusterOneLessFRsags>=5.9));
clusterTwoGreaterThanSagLessFR = clusterTwoLessFRsags(find(clusterTwoLessFRsags>=5.9));
clusterOneGreaterThanSagLessFRCounts = histc(clusterOneGreaterThanSagLessFR,edges.sag);
clusterTwoGreaterThanSagLessFRCounts = histc(clusterTwoGreaterThanSagLessFR,edges.sag);
sumGreaterThanSagLessFRCounts = clusterOneGreaterThanSagLessFRCounts + clusterTwoGreaterThanSagLessFRCounts;

%plot
histGreaterThanFig = figure(15);
histGreaterThanFig.Position = [670 150 300 135];
hold on
bar(edges.sag,sumGreaterThanSagLessFRCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.sag,clusterTwoGreaterThanSagLessFRCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.sag(1) edges.sag(end)])
greaterAx = gca;
setAx(greaterAx);
xlabel('sag (%)')

%% last cut is ahp latency
edges.ahplat = 0:1:55;

%get data by cluster, only with FR < 95Hz and sag < 5.9%
clusterOneLastCut = clusterOneNeurons(clusterOneLessThanSagLessFRInd);
clusterTwoLastCut = clusterTwoNeurons(clusterTwoLessThanSagLessFRInd);
clusterOneahpLat = ahpLat(clusterOneLastCut);
clusterTwoahpLat = ahpLat(clusterTwoLastCut);
clusterOneLastCutCounts = histc(clusterOneahpLat,edges.ahplat);
clusterTwoLastCutCounts = histc(clusterTwoahpLat,edges.ahplat);
sumCountsLastCut = clusterOneLastCutCounts + clusterTwoLastCutCounts;

% %plot
histAllFig = figure(9);
histAllFig.Position = [50 150 300 135];
hold on
bar(edges.ahplat,sumCountsLastCut,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
bar(edges.ahplat,clusterOneLastCutCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
xlim([edges.ahplat(1) edges.ahplat(end)])
allAx = gca;
setAx(allAx);
xlabel('AHP latency (ms)')

%< 8 ms only
%get data by cluster
clusterOneLessThanahp= clusterOneahpLat(find(clusterOneahpLat<8));
clusterTwoLessThanahp= clusterTwoahpLat(find(clusterTwoahpLat<8));
clusterOneLessThanahpCounts = histc(clusterOneLessThanahp,edges.ahplat);
clusterTwoLessThanahpCounts = histc(clusterTwoLessThanahp,edges.ahplat);
if isempty(clusterTwoLessThanahpCounts)
    sumLessThanCountsahp = clusterOneLessThanahpCounts;
else
    sumLessThanCountsahp = clusterOneLessThanahpCounts + clusterTwoLessThanahpCounts;
end

%plot
histLessThanFig = figure(10);
histLessThanFig.Position = [360 150 300 135];
hold on
bar(edges.ahplat,sumLessThanCountsahp,.85,'facecolor',grpOne,'edgecolor',grpOne)
if isempty(clusterTwoLessThanahpCounts) ~=1
    bar(edges.ahplat,clusterTwoLessThanahpCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
end
xlim([edges.ahplat(1) edges.ahplat(end)])
lessAx = gca;
setAx(lessAx);
xlabel('AHP latency (ms)')

%>= 8 ms only
%get data by cluster
clusterOneGreaterThanahp = clusterOneahpLat(find(clusterOneahpLat>=8));
clusterTwoGreaterThanahp = clusterTwoahpLat(find(clusterTwoahpLat>=8));
clusterOneGreaterThanahpCounts = histc(clusterOneGreaterThanahp,edges.ahplat);
clusterTwoGreaterThanahpCounts = histc(clusterTwoGreaterThanahp,edges.ahplat);
sumGreaterThanahpCounts = clusterOneGreaterThanahpCounts + clusterTwoGreaterThanahpCounts;

%plot
histGreaterThanFig = figure(11);
histGreaterThanFig.Position = [670 150 300 135];
hold on
bar(edges.ahplat,sumGreaterThanahpCounts,.85,'facecolor',grpOne,'edgecolor',grpOne)
bar(edges.ahplat,clusterTwoGreaterThanahpCounts,.85,'facecolor',grpTwo,'edgecolor',grpTwo)
xlim([edges.ahplat(1) edges.ahplat(end)])
greaterAx = gca;
setAx(greaterAx);
xlabel('AHP latency (ms)')

%% 3D scatter plot of the 3 Decision Tree Variables
scatterFig = figure (12);
scatterFig.Position = [1000 200 400 450];
hold on
scatter3(maxFR(clusterOneNeurons),sags(clusterOneNeurons),ahpLat(clusterOneNeurons),100,grpOne,'filled')
scatter3(maxFR(clusterTwoNeurons),sags(clusterTwoNeurons),ahpLat(clusterTwoNeurons),100,grpTwo,'filled')
ylabel('sag (%)')
ylim([0 ceil(1.15*max(sags))])
xlabel('max firing rate (Hz)')
xlim([0 ceil(1.15*max(maxFR))])
zlabel('AHP latency (ms)')
zlim([0 ceil(1.15*max(ahpLat))])
scatterAx = gca;
scatterAx.XGrid = 'on';
scatterAx.YGrid = 'on';
scatterAx.ZGrid = 'on';
scatterAx.ZTick = 0:10:ceil(1.15*max(ahpLat));
scatterAx.XTick = 0:40:ceil(1.15*max(maxFR));
scatterAx.YTick = 0:10:ceil(1.15*max(sags));