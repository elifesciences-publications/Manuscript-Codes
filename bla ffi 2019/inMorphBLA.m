%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% in morphology comp %%%%%%%%%%%
%%%%%%%%% Created: 04-14-2019 %%%%%%%%%%
%%%%%%%%%% Edited: 05-13-2019 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INIT VARS
close all
clearvars

fldLabels = {'PV', 'FS', 'nFS'};
shollRadii = 50:50:600;

%% LOAD DATA
%get files
%directory
morphDir = uigetdir;
cd(morphDir);
contents = dir('*.xlsx');
filenames = {contents.name}';
morphFiles = fullfile(cd,filenames);
for ii = 1:length(fldLabels)
    
    %load data
    if ii == 1 %pv
        %morphology
        somaPerim.(fldLabels{ii}) = xlsread(morphFiles{1},1,'C6:C11');
        somaArea.(fldLabels{ii}) = xlsread(morphFiles{1},1,'D6:D11');
        branchPoints.(fldLabels{ii}) = xlsread(morphFiles{1},2,'C5:C10');
        endPoints.(fldLabels{ii}) = xlsread(morphFiles{1},2,'D5:D10');
        dendriteLength.(fldLabels{ii}) = xlsread(morphFiles{1},2,'E5:E10');
        dendriteSA.(fldLabels{ii}) = xlsread(morphFiles{1},2,'G5:G10'); %surface area
        
        %sholl
        shollLength.(fldLabels{ii}) = zeros(6,length(shollRadii));
        shollLengthHold = xlsread(morphFiles{3},1,'C4:N10');
        for jj = 1:length(shollRadii)-1 %data only go to 550 instead of 600 here
            shollLength.(fldLabels{ii})(:,jj) = shollLengthHold(:,jj).*~isnan(shollLengthHold(:,jj));
        end
        shollLength.(fldLabels{ii})(isnan(shollLength.(fldLabels{ii}))) = 0;
        
        shollSA.(fldLabels{ii}) = zeros(6,length(shollRadii)); %surface area
        shollSAHold = xlsread(morphFiles{3},2,'C4:N10');
        for jj = 1:length(shollRadii)-1 %data only go to 550 instead of 600 here
            shollSA.(fldLabels{ii})(:,jj) = shollSAHold(:,jj).*~isnan(shollLengthHold(:,jj));
        end
        shollSA.(fldLabels{ii})(isnan(shollSA.(fldLabels{ii}))) = 0;
        
        clear shollSAHold shollLengthHold
    end
    
    if ii == 2 %fs sst
        %morphology
        somaPerim.(fldLabels{ii}) = xlsread(morphFiles{2},1,'C6:C13');
        somaArea.(fldLabels{ii}) = xlsread(morphFiles{2},1,'D6:D13');
        branchPoints.(fldLabels{ii}) = xlsread(morphFiles{2},2,'C5:C12');
        endPoints.(fldLabels{ii}) = xlsread(morphFiles{2},2,'D5:D12');
        dendriteLength.(fldLabels{ii}) = xlsread(morphFiles{2},2,'E5:E12');
        dendriteSA.(fldLabels{ii}) = xlsread(morphFiles{2},2,'G5:G12'); %surface area
        
        %sholl
        shollLength.(fldLabels{ii}) = zeros(8,length(shollRadii));
        shollLengthHold = xlsread(morphFiles{4},1,'C4:N11');
        for jj = 1:length(shollRadii)
            shollLength.(fldLabels{ii})(:,jj) = shollLengthHold(:,jj).*~isnan(shollLengthHold(:,jj));
        end
        shollLength.(fldLabels{ii})(isnan(shollLength.(fldLabels{ii}))) = 0;
        
        shollSA.(fldLabels{ii}) = zeros(8,length(shollRadii)); %surface area
        shollSAHold = xlsread(morphFiles{4},2,'C4:N11');
        for jj = 1:length(shollRadii)
            shollSA.(fldLabels{ii})(:,jj) = shollSAHold(:,jj).*~isnan(shollLengthHold(:,jj));
        end
        shollSA.(fldLabels{ii})(isnan(shollSA.(fldLabels{ii}))) = 0;
        
        clear shollSAHold shollLengthHold
    end
    
    if ii == 3 %nfs sst
        %morphology
        somaPerim.(fldLabels{ii}) = xlsread(morphFiles{2},1,'C14:C16');
        somaArea.(fldLabels{ii}) = xlsread(morphFiles{2},1,'D14:D16');
        branchPoints.(fldLabels{ii}) = xlsread(morphFiles{2},2,'C13:C15');
        endPoints.(fldLabels{ii}) = xlsread(morphFiles{2},2,'D13:D15');
        dendriteLength.(fldLabels{ii}) = xlsread(morphFiles{2},2,'E13:E15');
        dendriteSA.(fldLabels{ii}) = xlsread(morphFiles{2},2,'G13:G15'); %surface area
        
        %sholl
        shollLength.(fldLabels{ii}) = zeros(3,length(shollRadii));
        shollLengthHold = xlsread(morphFiles{4},1,'C12:N14');
        for jj = 1:length(shollRadii)
            if jj <= size(shollLengthHold,2)
                shollLength.(fldLabels{ii})(:,jj) = shollLengthHold(:,jj).*~isnan(shollLengthHold(:,jj));
            else
                shollLength.(fldLabels{ii})(:,jj) = zeros(3,1);
            end
        end
        shollLength.(fldLabels{ii})(isnan(shollLength.(fldLabels{ii}))) = 0;
        
        shollSA.(fldLabels{ii}) = zeros(3,length(shollRadii)); %surface area
        shollSAHold = xlsread(morphFiles{4},2,'C12:N14');
        for jj = 1:length(shollRadii)
            if jj <= size(shollSAHold,2)
                shollSA.(fldLabels{ii})(:,jj) = shollSAHold(:,jj).*~isnan(shollLengthHold(:,jj));
            else
                shollSA.(fldLabels{ii})(:,jj) = zeros(3,1);
            end
        end
        shollSA.(fldLabels{ii})(isnan(shollSA.(fldLabels{ii}))) = 0;
        
        clear shollSAHold shollLengthHold
    end
    
end

%% TEST DIFFERENCES
%adtest (Anderson-Darling Test) to test for normality
%soma data, normality
for ii = 1:2
    [holdH.(fldLabels{ii}).somaPerim,~] = adtest(somaPerim.(fldLabels{ii}));
    [holdH.(fldLabels{ii}).somaArea,~] = adtest(somaArea.(fldLabels{ii}));
end
hSoma.perimeter = sum([holdH.(fldLabels{1}).somaPerim holdH.(fldLabels{2}).somaPerim]);
hSoma.area = sum([holdH.(fldLabels{1}).somaArea holdH.(fldLabels{2}).somaArea]);
clear holdH

%soma data, comparison tests
%perimeter
if hSoma.perimeter == 0 %if groups normal, run ttest
    [~,p.somaPerimTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(somaPerim.(fldLabels{1}),somaPerim.(fldLabels{2}));
elseif hSoma.perimeter ~= 0 %if data not normal, run MWU test
    [p.somaPerimMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(somaPerim.(fldLabels{1}),somaPerim.(fldLabels{2}));
end
%area
if hSoma.area == 0 %if groups normal, run ttest
    [~,p.somaAreaTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(somaArea.(fldLabels{1}),somaArea.(fldLabels{2}));
elseif hSoma.area ~= 0 %if data not normal, run MWU test
    [p.somaAreaMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(somaArea.(fldLabels{1}),somaArea.(fldLabels{2}));
end

%dendrite data, normality
for ii = 1:2
    [holdH.(fldLabels{ii}).branchpoints,~] = adtest(branchPoints.(fldLabels{ii}));
    [holdH.(fldLabels{ii}).endpoints,~] = adtest(endPoints.(fldLabels{ii}));
    [holdH.(fldLabels{ii}).dndLength,~] = adtest(dendriteLength.(fldLabels{ii}));
    [holdH.(fldLabels{ii}).dndSA,~] = adtest(dendriteSA.(fldLabels{ii}));
end
hDnd.bp = sum([holdH.(fldLabels{1}).branchpoints holdH.(fldLabels{2}).branchpoints]);
hDnd.ep = sum([holdH.(fldLabels{1}).endpoints holdH.(fldLabels{2}).endpoints]);
hDnd.length = sum([holdH.(fldLabels{1}).branchpoints holdH.(fldLabels{2}).dndLength]);
hDnd.sa = sum([holdH.(fldLabels{1}).endpoints holdH.(fldLabels{2}).dndSA]);
clear holdH

%dendrite data, comparison tests
%branchpoints
if hDnd.bp == 0 %if groups normal, run ttest
    [~,p.dndBranchTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(branchPoints.(fldLabels{1}),branchPoints.(fldLabels{2}));
elseif hDnd.bp ~= 0 %if data not normal, run MWU test
    [p.dndBranchMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(branchPoints.(fldLabels{1}),branchPoints.(fldLabels{2}));
end
%endpoints
if hDnd.ep == 0 %if groups normal, run ttest
    [~,p.dndEndptsTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(endPoints.(fldLabels{1}),endPoints.(fldLabels{2}));
elseif hDnd.ep ~= 0 %if data not normal, run MWU test
    [p.dndEndptshMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(endPoints.(fldLabels{1}),endPoints.(fldLabels{2}));
end
%length
if hDnd.length == 0 %if groups normal, run ttest
    [~,p.dndLengthTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(dendriteLength.(fldLabels{1}),dendriteLength.(fldLabels{2}));
elseif hDnd.length ~= 0 %if data not normal, run MWU test
    [p.dndLengthMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(dendriteLength.(fldLabels{1}),dendriteLength.(fldLabels{2}));
end
%surface area
if hDnd.sa == 0 %if groups normal, run ttest
    [~,p.dndSATtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(dendriteSA.(fldLabels{1}),dendriteSA.(fldLabels{2}));
elseif hDnd.sa ~= 0 %if data not normal, run MWU test
    [p.dndSAMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(dendriteSA.(fldLabels{1}),dendriteSA.(fldLabels{2}));
end

%% DISPLAY DATA
%color
lightgray=[.75 .75 .75]; %light gray for individual data points
pColor = [.35 .6 .25; .6 .4 .67; .9686 .6588 .7216];
nfsColor = [247 168 184]/255;

%soma figure, perimeter
somaFig.perimeter=figure(1);
somaFig.perimeter.Position = [415 545 165 250];
hold on
for ii = 1:2
    scatter(ii.*ones(1,length(somaPerim.(fldLabels{ii}))),somaPerim.(fldLabels{ii}),40,lightgray,'filled')
    if hSoma.perimeter == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
        errorbar(ii,mean(somaPerim.(fldLabels{ii})),sem(somaPerim.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(somaPerim.(fldLabels{ii})),125,pColor(ii,:),'filled')
    elseif hSoma.perimeter ~= 0 %if data non-normal, plot median w IQR
        [lowerError.somaPerim.(fldLabels{ii}),upperError.somaPerim.(fldLabels{ii})] = iqrError(somaPerim.(fldLabels{ii}),1);
        errorbar(ii,median(somaPerim.(fldLabels{ii})),lowerError.somaPerim.(fldLabels{ii}),upperError.somaPerim.(fldLabels{ii}),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(somaPerim.(fldLabels{ii})),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allSoma.perimeter = somaPerim.(fldLabels{ii});
    else
        allSoma.perimeter = [allSoma.perimeter; somaPerim.(fldLabels{ii})];
    end
end
ylabel('soma perimeter (\mum)')
somaAx.perimeter=gca;
setAx(somaAx.perimeter);
somaAx.perimeter.XTick = [1 2];
somaAx.perimeter.XTickLabel = {fldLabels{1} fldLabels{2}};
xlim([0.5 2.5])
ylim([0 ceil(1.25*max(allSoma.perimeter))])

%nfs
nfsSomaFig.perimeter=figure(9);
nfsSomaFig.perimeter.Position = [415 545 100 250];
hold on
for ii = 1:3
    scatter(((rand(1)+9.5)/10),somaPerim.(fldLabels{3})(ii),125,nfsColor,'filled')
end
ylabel('soma perimeter (\mum)')
nfsSomaAx.perimeter=gca;
setAx(nfsSomaAx.perimeter);
nfsSomaAx.perimeter.XTick = 1;
nfsSomaAx.perimeter.XTickLabel = fldLabels{3};
xlim([.5 1.5])
ylim(somaAx.perimeter.YLim)

%soma figure, area
somaFig.area=figure(2);
somaFig.area.Position = [590 545 165 250];
hold on
for ii = 1:2
    scatter(ii.*ones(1,length(somaArea.(fldLabels{ii}))),somaArea.(fldLabels{ii}),40,lightgray,'filled')
    if hSoma.perimeter == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
        errorbar(ii,mean(somaArea.(fldLabels{ii})),sem(somaArea.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(somaArea.(fldLabels{ii})),125,pColor(ii,:),'filled')
    elseif hSoma.perimeter ~= 0 %if data non-normal, plot median w IQR
        [lowerError.somaArea.(fldLabels{ii}),upperError.somaArea.(fldLabels{ii})] = iqrError(somaArea.(fldLabels{ii}),1);
        errorbar(ii,median(somaArea.(fldLabels{ii})),lowerError.somaArea.(fldLabels{ii}),upperError.somaArea.(fldLabels{ii}),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(somaArea.(fldLabels{ii})),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allSoma.area = somaArea.(fldLabels{ii});
    else
        allSoma.area = [allSoma.area; somaArea.(fldLabels{ii})];
    end
end
ylabel('soma area (\mum^{2})')
somaAx.area=gca;
setAx(somaAx.area);
somaAx.area.XTick = [1 2];
somaAx.area.XTickLabel = {fldLabels{1} fldLabels{2}};
xlim([0.5 2.5])
ylim([0 ceil(1.25*max(allSoma.area))])

%nfs
nfsSomaFig.area=figure(10);
nfsSomaFig.area.Position = [590 545 100 250];
hold on
for ii = 1:3
    scatter(((rand(1)+9.5)/10),somaArea.(fldLabels{3})(ii),125,nfsColor,'filled')
end
ylabel('soma area (\mum^{2})')
nfsSomaAx.area=gca;
setAx(nfsSomaAx.area);
nfsSomaAx.area.XTick = 1;
nfsSomaAx.area.XTickLabel = fldLabels{3};
xlim([.5 1.5])
ylim(somaAx.area.YLim)

%dendrite figure, branchpoints
dndFig.branchpoints=figure(3);
dndFig.branchpoints.Position = [765 545 165 250];
hold on
for ii = 1:2
    scatter(ii.*ones(1,length(branchPoints.(fldLabels{ii}))),branchPoints.(fldLabels{ii}),40,lightgray,'filled')
    if hDnd.bp == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
        errorbar(ii,mean(branchPoints.(fldLabels{ii})),sem(branchPoints.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(branchPoints.(fldLabels{ii})),125,pColor(ii,:),'filled')
    elseif hDnd.bp ~= 0 %if data non-normal, plot median w IQR
        [lowerError.branchPoints.(fldLabels{ii}),upperError.branchPoints.(fldLabels{ii})] = iqrError(branchPoints.(fldLabels{ii}),1);
        errorbar(ii,median(branchPoints.(fldLabels{ii})),lowerError.branchPoints.(fldLabels{ii}),upperError.branchPoints.(fldLabels{ii}),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(branchPoints.(fldLabels{ii})),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allDnd.bp = branchPoints.(fldLabels{ii});
    else
        allDnd.bp = [allDnd.bp; branchPoints.(fldLabels{ii})];
    end
end
ylabel('branchpoints')
dndAx.branchpoints=gca;
setAx(dndAx.branchpoints);
dndAx.branchpoints.XTick = [1 2];
dndAx.branchpoints.XTickLabel = {fldLabels{1} fldLabels{2}};
xlim([0.5 2.5])
ylim([0 ceil(1.25*max(allDnd.bp))])

%dnfs
nfsDndFig.branchpoints=figure(11);
nfsDndFig.branchpoints.Position = [765 545 100 250];
hold on
for ii = 1:3
    scatter((rand(1)+9.5)/10,branchPoints.(fldLabels{3})(ii),125,nfsColor,'filled')
end
ylabel('branchpoints')
nfsDndAx.branchpoints=gca;
setAx(nfsDndAx.branchpoints);
nfsDndAx.branchpoints.XTick = 1;
dndAx.branchpoints.XTickLabel = fldLabels{3};
xlim([0.5 1.5])
ylim(dndAx.branchpoints.YLim)

%dendrite figure, endpoints
dndFig.endpoints=figure(4);
dndFig.endpoints.Position = [940 545 165 250];
hold on
for ii = 1:2
    scatter(ii.*ones(1,length(endPoints.(fldLabels{ii}))),endPoints.(fldLabels{ii}),40,lightgray,'filled')
    if hDnd.ep == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
        errorbar(ii,mean(endPoints.(fldLabels{ii})),sem(endPoints.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(endPoints.(fldLabels{ii})),125,pColor(ii,:),'filled')
    elseif hDnd.ep ~= 0 %if data non-normal, plot median w IQR
        [lowerError.endPoints.(fldLabels{ii}),upperError.endPoints.(fldLabels{ii})] = iqrError(endPoints.(fldLabels{ii}),1);
        errorbar(ii,median(endPoints.(fldLabels{ii})),lowerError.endPoints.(fldLabels{ii}),upperError.endPoints.(fldLabels{ii}),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(endPoints.(fldLabels{ii})),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allDnd.ep = endPoints.(fldLabels{ii});
    else
        allDnd.ep = [allDnd.ep; endPoints.(fldLabels{ii})];
    end
end
ylabel('endpoints')
dndAx.endpoints=gca;
setAx(dndAx.endpoints);
dndAx.endpoints.XTick = [1 2];
dndAx.endpoints.XTickLabel = {fldLabels{1} fldLabels{2}};
xlim([0.5 2.5])
ylim([0 ceil(1.25*max(allDnd.ep))])

%nfs
nfsDndFig.endpoints=figure(12);
nfsDndFig.endpoints.Position = [940 545 100 250];
hold on
for ii = 1:3
    scatter((rand(1)+9.5)/10,endPoints.(fldLabels{3})(ii),125,nfsColor,'filled')
end
ylabel('endpoints')
nfsDndAx.endpoints=gca;
setAx(nfsDndAx.endpoints);
nfsDndAx.endpoints.XTick = 1;
nfsDndAx.endpoints.XTickLabel = fldLabels{3};
xlim([0.5 1.5])
ylim(dndAx.endpoints.YLim)

%dendrite figure, length
dndFig.length=figure(5);
dndFig.length.Position = [1115 545 165 250];
hold on
for ii = 1:2
    scatter(ii.*ones(1,length(dendriteLength.(fldLabels{ii}))),dendriteLength.(fldLabels{ii}),40,lightgray,'filled')
    if hDnd.length == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
        errorbar(ii,mean(dendriteLength.(fldLabels{ii})),sem(dendriteLength.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(dendriteLength.(fldLabels{ii})),125,pColor(ii,:),'filled')
    elseif hDnd.length ~= 0 %if data non-normal, plot median w IQR
        [lowerError.dendriteLength.(fldLabels{ii}),upperError.dendriteLength.(fldLabels{ii})] = iqrError(dendriteLength.(fldLabels{ii}),1);
        errorbar(ii,median(dendriteLength.(fldLabels{ii})),lowerError.dendriteLength.(fldLabels{ii}),upperError.dendriteLength.(fldLabels{ii}),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(dendriteLength.(fldLabels{ii})),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allDnd.length = dendriteLength.(fldLabels{ii});
    else
        allDnd.length = [allDnd.length; dendriteLength.(fldLabels{ii})];
    end
end
ylabel('dendrite length (\mum)')
dndAx.length=gca;
setAx(dndAx.length);
dndAx.length.XTick = [1 2];
dndAx.length.XTickLabel = {fldLabels{1} fldLabels{2}};
xlim([0.5 2.5])
ylim([0 ceil(1.25*max(allDnd.length))])

%nfs
nfsDndFig.length=figure(13);
nfsDndFig.length.Position = [1115 545 100 250];
hold on
for ii = 1:3
    scatter((rand(1)+9.5)/10,dendriteLength.(fldLabels{3})(ii),125,nfsColor,'filled')
end
ylabel('dendrite length (\mum)')
nfsDndAx.length=gca;
setAx(nfsDndAx.length);
nfsDndAx.length.XTick = 1;
nfsDndAx.length.XTickLabel = fldLabels{3};
xlim([0.5 1.5])
ylim(dndAx.length.YLim)

%dendrite figure, surface area
dndFig.sa=figure(6);
dndFig.sa.Position = [1290 545 165 250];
hold on
for ii = 1:2
    scatter(ii.*ones(1,length(dendriteSA.(fldLabels{ii}))),dendriteSA.(fldLabels{ii}),40,lightgray,'filled')
    if hDnd.sa == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
        errorbar(ii,mean(dendriteSA.(fldLabels{ii})),sem(dendriteSA.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(dendriteSA.(fldLabels{ii})),125,pColor(ii,:),'filled')
    elseif hDnd.sa ~= 0 %if data non-normal, plot median w IQR
        [lowerError.dendriteSA.(fldLabels{ii}),upperError.dendriteSA.(fldLabels{ii})] = iqrError(dendriteSA.(fldLabels{ii}),1);
        errorbar(ii,median(dendriteSA.(fldLabels{ii})),lowerError.dendriteSA.(fldLabels{ii}),upperError.dendriteSA.(fldLabels{ii}),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(dendriteSA.(fldLabels{ii})),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allDnd.sa = dendriteSA.(fldLabels{ii});
    else
        allDnd.sa = [allDnd.sa; dendriteSA.(fldLabels{ii})];
    end
end
ylabel('dendrite surface area (\mum^{2})')
dndAx.sa=gca;
setAx(dndAx.sa);
dndAx.sa.XTick = [1 2];
dndAx.sa.XTickLabel = {fldLabels{1} fldLabels{2}};
xlim([0.5 2.5])
ylim([0 ceil(1.25*max(allDnd.sa))])

%nfs
nfsDndFig.sa=figure(14);
nfsDndFig.sa.Position = [1290 545 100 250];
hold on
for ii = 1:3
    scatter((rand(1)+9.5)/10,dendriteSA.(fldLabels{3})(ii),125,nfsColor,'filled')
end
ylabel('dendrite surface area (\mum^{2})')
nfsDndAx.sa=gca;
setAx(nfsDndAx.sa);
nfsDndAx.sa.XTick = 1;
nfsDndAx.sa.XTickLabel = fldLabels{3};
xlim([0.5 1.5])
ylim(dndAx.sa.YLim)

%sholl analyses
%zeros for NaN
%dendrite length
shollFig.length=figure(7);
shollFig.length.Position = [415 100 415 165];
hold on
for ii = 1:2
    errorbar(shollRadii,mean(shollLength.(fldLabels{ii})),sem(shollLength.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
end
ylabel('dendrite length (\mum)')
shollAx.length=gca;
setAx(shollAx.length);

%dendrite surface area
shollFig.sa=figure(8);
shollFig.sa.Position = [850 100 415 165];
hold on
for ii = 1:2
    errorbar(shollRadii,mean(shollSA.(fldLabels{ii})),sem(shollSA.(fldLabels{ii}),1),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
end
ylabel('dendrite surface area (\mum^{2})')
shollAx.sa=gca;
setAx(shollAx.sa);

%nfs
%dendrite length
nfsShollFig.length=figure(15);
nfsShollFig.length.Position = [415 100 415 165];
hold on
for ii = 1:3
    plot(shollRadii,shollLength.(fldLabels{3})(ii,:),'color',nfsColor,'linewidth',2)
end
ylabel('dendrite length (\mum)')
nfsShollAx.length=gca;
setAx(nfsShollAx.length);

%dendrite surface area
nfsShollFig.sa=figure(16);
nfsShollFig.sa.Position = [850 100 415 165];
hold on
for ii = 1:3
    plot(shollRadii,shollSA.(fldLabels{3})(ii,:),'color',nfsColor,'linewidth',2)
end
ylabel('dendrite surface area (\mum^{2})')
nfsShollAx.sa=gca;
setAx(nfsShollAx.sa);