%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% drug change in Vm display/stats %%%
%%%%%%%%% Created: 09-21-2016 %%%%%%%%%%
%%%%%%%%%% Edited: 01-12-2016 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INIT VARS
close all
clearvars
t=(-3):(1/2):9.5;

%% LOAD DATA
%sst
vm_dir=uigetdir;
cd(vm_dir);
contents = dir('*.mat');
filenames = {contents.name}';
vm_files.sst = fullfile(cd,filenames);

% %pv
vm_dir=uigetdir;
cd(vm_dir);
contents = dir('*.mat');
filenames = {contents.name}';
vm_files.pv = fullfile(cd,filenames);

for ii = 1:length(vm_files.sst)
    %load sst data
    load(vm_files.sst{ii}, 'dVm', 'baseVm','mean_rawdVm')
    deltaV.sst(ii,:)=dVm; dVm=[];
    baseV.sst(ii,:)=baseVm; baseVm=[];
    mean_dV.sst(ii)=mean_rawdVm; mean_rawdVm=[];
end

for ii = 1:length(vm_files.pv)
    %load pv data
    load(vm_files.pv{ii}, 'dVm', 'baseVm','mean_rawdVm')
    deltaV.pv(ii,:)=dVm; dVm=[];
    baseV.pv(ii,:)=baseVm; baseVm=[];
    mean_dV.pv(ii)=mean_rawdVm; mean_rawdVm=[];
end

%% DATA ANALYSIS -- edit once I have n > 4
%check normality
%sst
[h_dV.sst,~] = adtest(mean_dV.sst);

%check for difference in base vs drug
if h_dV.sst == 0
    [~,p_sst]=ttest(mean_dV.sst)
else
    p_sst_wilc=signrank(mean_dV.sst)
end

%pv
[h_dV.pv,~] = adtest(mean_dV.pv);

%check for difference in base vs drug
if h_dV.pv == 0
    [~,p_pv]=ttest(mean_dV.pv) %#ok<*NOPTS>
else
    p_pv_wilc=signrank(mean_dV.pv)
end


%% DISPLAY DATA
%colors
%pick color
gray=[.33 .33 .33];
lightgray=[.75 .75 .75];
orange=[.5 .33 .55]; %orange
dorange=[.4 .2 .05]; %darker orange
green=[.35 .6 .25]; %green
dgreen=[.1 .3 .15]; %darker green

%sst figure
sst_dVFig=figure(1);
sst_dVFig.Position=[800 550 625 250];
subplot(1,2,1)
hold on
line([t(7) t(end)],[6 6],'linewidth',2,'color','k')
line([-3 t(end)],[0 0],'linewidth',1,'color','k','linestyle','--')
plot(t,[mean(baseV.sst) mean(deltaV.sst)],'color',orange,'markerfacecolor',orange,'markeredgecolor',orange,'linewidth',1.5);
for ii = 1:6
    baseSem.sst(ii)=sem(baseV.sst(:,ii),1);
end
for ii = 1:20
    deltaSem.sst(ii)=sem(deltaV.sst(:,ii),1);
end
errorbar(t,[mean(baseV.sst) mean(deltaV.sst)],[baseSem.sst deltaSem.sst],'color',orange,'linewidth',1.5)
xlim([-3 10.5])
ylim([-22.5 7.5])
xlabel('Time (min)')
ylabel('\Delta Vm (mV)')
subplot(1,2,2)
hold on
bar(mean(mean_dV.sst),'facecolor',orange,'edgecolor',orange)
errorbar(mean(mean_dV.sst),sem(mean_dV.sst,2),'color',orange,'linewidth',2)
scatter(ones(length(vm_files.sst),1),mean_dV.sst,'markerfacecolor','k','markeredgecolor','k')
xlim([0 2])
ylim([-22.5 7.5])
sstAx=gca;
sstAx.XAxisLocation='origin';
sstAx.XTickLabel=[];
sstAx.XTick=[];
ylabel('\Delta Vm (mV)')

%pv figure
pv_dVFig=figure(2);
pv_dVFig.Position=[800 200 625 250];
subplot(1,2,1)
hold on
line([t(7) t(end)],[6 6],'linewidth',2,'color','k')
line([-3 t(end)],[0 0],'linewidth',1,'color','k','linestyle','--')
plot(t,[mean(baseV.pv) mean(deltaV.pv)],'color',green,'markerfacecolor',green,'markeredgecolor',green,'linewidth',1.5);
for ii = 1:6
    baseSem.pv(ii)=sem(baseV.pv(:,ii),1);
end
for ii = 1:20
    deltaSem.pv(ii)=sem(deltaV.pv(:,ii),1);
end
errorbar(t,[mean(baseV.pv) mean(deltaV.pv)],[baseSem.pv deltaSem.pv],'color',green,'linewidth',1.5)
xlim([-3 10.5])
ylim([-22.5 7.5])
xlabel('Time (min)')
ylabel('\Delta Vm (mV)')
subplot(1,2,2)
hold on
bar(mean(mean_dV.pv),'facecolor',green,'edgecolor',green)
errorbar(mean(mean_dV.pv),sem(mean_dV.pv,2),'color',green,'linewidth',2)
scatter(ones(length(vm_files.pv),1),mean_dV.pv,'markerfacecolor','k','markeredgecolor','k')
xlim([0 2])
ylim([-22.5 7.5])
pvAx=gca;
pvAx.XAxisLocation='origin';
pvAx.XTickLabel=[];
pvAx.XTick=[];
ylabel('\Delta Vm (mV)')

%comparison figure
compFig.tPlot=figure(3);
compFig.tPlot.Position=[50 650 275 125];

%time course graph
hold on
line([t(7) t(end)],[6 6],'linewidth',1.5,'color','k')
line([-3 t(end)],[0 0],'linewidth',1.5,'color','k','linestyle','--')
%scatter and errorbars
scatter(t,[mean(baseV.pv) mean(deltaV.pv)],80,green,'filled')
errorbar(t,[mean(baseV.pv) mean(deltaV.pv)],[baseSem.pv deltaSem.pv],'color',green,'linewidth',2,'linestyle','none','capsize',0)
scatter(t,[mean(baseV.sst) mean(deltaV.sst)],80,orange,'filled')
errorbar(t,[mean(baseV.sst) mean(deltaV.sst)],[baseSem.sst deltaSem.sst],'color',orange,'linewidth',2,'linestyle','none','capsize',0)

xlim([-3 10.5])
ylim([-25 7.5])
tAx=gca;
setAx(tAx);

% summary graph
compFig.dPlot=figure(4);
compFig.dPlot.Position=[350 650 100 125];
hold on
scatter(ones(length(mean_dV.sst),1),mean_dV.sst,40,lightgray,'filled')
scatter(2.*ones(length(mean_dV.pv),1),mean_dV.pv,40,lightgray,'filled')
if sum([h_dV.sst h_dV.pv]) == 0 %if data normal, pliot as mean +/- Sem
    errorbar(1,mean(mean_dV.sst),sem(mean_dV.sst,2),'color',orange,'linewidth',4)
    errorbar(2,mean(mean_dV.pv),sem(mean_dV.pv,2),'color',green,'linewidth',4)
    scatter(1,mean(mean_dV.sst),125,orange,'filled')
    scatter(2,mean(mean_dV.pv),125,green,'filled')
else %plot as median w IQR
    errorbar(1,median(mean_dV.sst),median(mean_dV.sst)-quantile(mean_dV.sst,.25),quantile(mean_dV.sst,.75)-median(mean_dV.sst),'color',gray,'linewidth',1.5)
    errorbar(2,median(mean_dV.pv),median(mean_dV.pv)-quantile(mean_dV.pv,.25),quantile(mean_dV.pv,.75)-median(mean_dV.pv),'color',color,'linewidth',1.5)
    scatter(1,median(mean_dV.sst),75,orange,'filled')
    scatter(2,median(mean_dV.pv),75,green,'filled')
end
ylim([-25 7.5])
xlim([.5 2.5])
sAx=gca;
sAx.LineWidth=1; sAx.XAxisLocation='origin'; sAx.YAxisLocation='origin'; sAx.TickDir='out'; sAx.Box ='off'; sAx.XColor='k'; sAx.YColor='k';