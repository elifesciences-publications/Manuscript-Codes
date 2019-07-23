%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% stim. induced spiking comp for Suite %%%
%%%%%%%%%%% Created: 05-12-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 05-12-2018 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = stimSpkCompFxn(inputs)
%modeled after pn_excite_stats.m // sst_pv_excite.m

%%%% make so functional for interneuron analysis too
%inits
samplerate=inputs.SamplerateHzEditField.Value;
noGroups = sum([inputs.ExpgrponeCheckBox.Value inputs.ExpgrptwoCheckBox.Value]);
for ii = 1:noGroups
    if ii == 1
        grpLabels{ii} = inputs.Grp1EditField.Value;
    elseif ii == 2
        grpLabels{ii} = inputs.Grp2EditField.Value;
    end
end
for ii = 1:noGroups
    fldLabels{ii} = grpLabels{ii};
    dashes = find(fldLabels{ii}=='-');
    fldLabels{ii}(dashes) = [];
    clear dashes
end

%% LOAD DATA
% Directory
cd(inputs.stimSpkDir);
for ii = 1:noGroups
    %get files,
    if inputs.PairedAnalysisCheckBox.Value == 1
        cd([inputs.stimSpkDir,'/',grpLabels{ii}]);
    elseif inputs.PairedAnalysisCheckBox.Value ~= 1
        thisDir = uigetdir;
        cd(thisDir);
    end
    contents.(fldLabels{ii}) = dir('*.mat');
    filenames.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
    stimSpkFiles.(fldLabels{ii}) = fullfile(cd,filenames.(fldLabels{ii}));
    
    %load
    for jj = 1:length(stimSpkFiles.(fldLabels{ii}))
        load(stimSpkFiles.(fldLabels{ii}){jj},'app')
        Vm.(fldLabels{ii})(jj) = app.stimSpkData.blVm;
        spkPerStim.byStim.(fldLabels{ii})(jj,:) = app.stimSpkData.mSpkPerStim.byStim;
        spkPerStim.bySweep.(fldLabels{ii})(jj) = app.stimSpkData.mSpkPerStim.bySweep;
        vTraces.(fldLabels{ii})(jj,:,:) = app.stimSpkData.vData - app.stimSpkData.baseVm(jj);
        tstim.(fldLabels{ii})(jj,:) = app.stimSpkData.stimstart;
        if jj == 1
            allstEPSP.(fldLabels{ii})(:) = app.stimSpkData.stEPSP;
        else
            allstEPSP.(fldLabels{ii}) = [allstEPSP.(fldLabels{ii}) app.stimSpkData.stEPSP];
        end
        cellMeanstEPSP.(fldLabels{ii})(jj) = mean(app.stimSpkData.stEPSP);
        if ii == 1
            if inputs.ExCellGrp1EditField.Value == jj
                mstTrace.(fldLabels{ii}) = app.stimSpkData.mEPSPtrace';
            end
        elseif ii == 2
            if inputs.ExCellGrp2EditField.Value == jj
                mstTrace.(fldLabels{ii}) = app.stimSpkData.mEPSPtrace';
            end
        end
    end
    cellMeanstEPSP.(fldLabels{ii})(isnan(cellMeanstEPSP.(fldLabels{ii}))) = [];
end
if inputs.PairedAnalysisCheckBox.Value == 1
    if inputs.ExCellGrp1EditField.Value ~= inputs.ExCellGrp2EditField.Value
        error('for paired analysis, pick the same cell')
    end
end

%% ANALYZE DATA
%check for differences in membrane voltage
%normality
for ii = 1:noGroups
    [hHold.(fldLabels{ii}),~]=adtest(Vm.(fldLabels{ii}));
end
if noGroups == 1
    hVm = hHold.(fldLabels{1});
elseif noGroups == 2
    hVm = sum([hHold.(fldLabels{1}) hHold.(fldLabels{2})]);
end
clear hHold

%comparisons
if noGroups == 2
    if inputs.PairedAnalysisCheckBox.Value == 1
        if hVm == 0 %if data are normal, run paired ttest
            [~,output.p.Vmpairedttest]=ttest(Vm.(fldLabels{1}),Vm.(fldLabels{2}));
        else %if data are not normal, run wilcoxon signed rank
            [output.p.VmWSR,~]=signrank(Vm.(fldLabels{1}),Vm.(fldLabels{2}));
        end
    else
        if hVm == 0 %if data are normal, run unpaired ttest
            [~,output.p.Vmunpairedttest]=ttest2(Vm.(fldLabels{1}),Vm.(fldLabels{2}));
        else %if data are not normal, run mann whitney u test
            [output.p.VmMW,~]=ranksum(Vm.(fldLabels{1}),Vm.(fldLabels{2}));
        end
    end
end

%spikes per stimulus
%normality
for ii = 1:noGroups
    [hHold.(fldLabels{ii}),~]=adtest(spkPerStim.bySweep.(fldLabels{ii}));
end
if noGroups == 1
    hSPS = hHold.(fldLabels{1});
elseif noGroups == 2
    hSPS = sum([hHold.(fldLabels{1}) hHold.(fldLabels{2})]);
end
clear hHold

%comparisons
if noGroups == 2
    if inputs.PairedAnalysisCheckBox.Value == 1
        if hSPS == 0 %if data are normal, run paired ttest
            [~,output.p.SPSpairedttest]=ttest(spkPerStim.bySweep.(fldLabels{1}),spkPerStim.bySweep.(fldLabels{2}));
        else %if data are not normal, run wilcoxon signed rank
            [output.p.SPSWSR,~]=signrank(spkPerStim.bySweep.(fldLabels{1}),spkPerStim.bySweep.(fldLabels{2}));
        end
    else
        if hSPS == 0 %if data are normal, run unpaired ttest
            [~,output.p.SPSunpairedttest]=ttest2(spkPerStim.bySweep.(fldLabels{1}),spkPerStim.bySweep.(fldLabels{2}));
        else %if data are not normal, run mann whitney u test
            [output.p.SPSMW,~]=ranksum(spkPerStim.bySweep.(fldLabels{1}),spkPerStim.bySweep.(fldLabels{2}));
        end
    end
end

%subthreshold EPSP amplitude
%normality, cellular means
for ii = 1:noGroups
    [hHold.(fldLabels{ii}),~]=adtest(cellMeanstEPSP.(fldLabels{ii}));
end
if noGroups == 1
    hstEPSP = hHold.(fldLabels{1});
elseif noGroups == 2
    hstEPSP = sum([hHold.(fldLabels{1}) hHold.(fldLabels{2})]);
end
clear hHold

%comparisons, cellular means
if noGroups == 2
    if inputs.PairedAnalysisCheckBox.Value == 1
        if hstEPSP == 0 %if data are normal, run paired ttest
            [~,output.p.stEPSPpairedttest]=ttest(cellMeanstEPSP.(fldLabels{1}),cellMeanstEPSP.(fldLabels{2}));
        else %if data are not normal, run wilcoxon signed rank
            [output.p.stEPSPWSR,~]=signrank(cellMeanstEPSP.(fldLabels{1}),cellMeanstEPSP.(fldLabels{2}));
        end
    else
        if hstEPSP == 0 %if data are normal, run unpaired ttest
            [~,output.p.stEPSPunpairedttest]=ttest2(cellMeanstEPSP.(fldLabels{1}),cellMeanstEPSP.(fldLabels{2}));
        else %if data are not normal, run mann whitney u test
            [output.p.stEPSPMW,~]=ranksum(cellMeanstEPSP.(fldLabels{1}),cellMeanstEPSP.(fldLabels{2}));
        end
    end
end

%distribution of all stEPSP and k-s test
for ii = 1:noGroups
    [fEPSP.(fldLabels{ii}),xEPSP.(fldLabels{ii})] = ecdf(allstEPSP.(fldLabels{ii}));
end
if noGroups == 2
    [~,output.p.stEPSPks] = kstest2(allstEPSP.(fldLabels{1}),allstEPSP.(fldLabels{2}));
end

%% PLOT DATA
%additional colors
lightgray=[.75 .75 .75]; %light gray for individual data points
darkgray=[.33 .33 .33]; %dark gray for connecting mean/medians in paired analysis

%example traces
%inits
for ii = 1:noGroups
    if ii == 1
        exNeuron(ii) = inputs.ExCellGrp1EditField.Value;
    elseif ii == 2
        exNeuron(ii) = inputs.ExCellGrp2EditField.Value;
    end
end
t.vData = 1000.*(1/samplerate:1/samplerate:size(vTraces.(fldLabels{1})(exNeuron(ii),:,:),2)/samplerate); %time in ms
t.stEPSP = 1000.*(1/samplerate:1/samplerate:length(mstTrace.(fldLabels{1}))/samplerate);

%plot traces
for ii = 1:noGroups
    tracesFig.(fldLabels{ii})=figure(ii);
    tracesFig.(fldLabels{ii}).Position=[100+(ii-1)*410 125 400 225];
    hold on
    for jj = 1:size(vTraces.(fldLabels{1})(exNeuron(ii),:,:),3)
        if ii == 1
            plot(t.vData,vTraces.(fldLabels{ii})(exNeuron(ii),:,jj)+(jj-1)*round(100/size(vTraces.(fldLabels{1})(exNeuron(ii),:,:),3)),'color',inputs.plotColorGrpOne,'linewidth',1.5)
        elseif ii == 2
            plot(t.vData,vTraces.(fldLabels{ii})(exNeuron(ii),:,jj)+(jj-1)*round(100/size(vTraces.(fldLabels{1})(exNeuron(ii),:,:),3)),'color',inputs.plotColorGrpTwo,'linewidth',1.5)
        end
    end
    xlim((1000/samplerate).*[tstim.(fldLabels{ii})(exNeuron(ii),1)-.05*samplerate tstim.(fldLabels{ii})(exNeuron(ii),end)+.1*samplerate])
    ylabel('Vm (mV // trials separated)')
    xlabel('time (ms)')
    title(fldLabels{ii})
    tracesAx.(fldLabels{ii}) = gca;
    setAx(tracesAx.(fldLabels{ii}));
    ylims(ii,:) = tracesAx.(fldLabels{ii}).YLim;
end
for ii = 1:noGroups
    figure(ii)
    ylim([.5*min(ylims(:)) max(ylims(:))])
    if tracesAx.(fldLabels{ii}).YLim(2) < 110
        tracesAx.(fldLabels{ii}).YLim(2) = 110;
    end
end

%spikes per stimulus, by stim
spkPerStimFig.byStim = figure(noGroups+1);
spkPerStimFig.byStim.Position = [450 450 175 275];
hold on
for ii = 1:noGroups
    if hSPS == 0 %if data normal, plot mean + sem
        if ii == 1
            errorbar(1:1:size(spkPerStim.byStim.(fldLabels{ii}),2),mean(spkPerStim.byStim.(fldLabels{ii}),1),sem(spkPerStim.byStim.(fldLabels{ii}),1),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            plot(mean(spkPerStim.byStim.(fldLabels{ii}),1),'linewidth',3,'color',inputs.plotColorGrpOne)
        elseif ii == 2
            errorbar(1:1:size(spkPerStim.byStim.(fldLabels{ii}),2),mean(spkPerStim.byStim.(fldLabels{ii}),1),sem(spkPerStim.byStim.(fldLabels{ii}),1),'linewidth',3,'color',inputs.plotColorGrpTwo,'capsize',0)
            plot(mean(spkPerStim.byStim.(fldLabels{ii}),1),'linewidth',3,'color',inputs.plotColorGrpTwo)
        end
    elseif hSPS ~=0 %if data not normal, plot median + iqr
        [lowerError.spkPerStim.byStim.(fldLabels{ii}),upperError.spkPerStim.byStim.(fldLabels{ii})] = iqrError(spkPerStim.byStim.(fldLabels{ii}),1);
        if ii == 1
            errorbar(1:1:size(spkPerStim.byStim.(fldLabels{ii}),2),median(spkPerStim.byStim.(fldLabels{ii}),1),lowerError.spkPerStim.byStim.(fldLabels{ii}),upperError.spkPerStim.byStim.(fldLabels{ii}),'color',inputs.plotColorGrpOne,'linewidth',3,'CapSize',0)
            plot(median(spkPerStim.byStim.(fldLabels{ii}),1),'linewidth',3,'color',inputs.plotColorGrpOne)
        elseif ii == 2
            errorbar(1:1:size(spkPerStim.byStim.(fldLabels{ii}),2),median(spkPerStim.byStim.(fldLabels{ii}),1),lowerError.spkPerStim.byStim.(fldLabels{ii}),upperError.spkPerStim.byStim.(fldLabels{ii}),'color',inputs.plotColorGrpTwo,'linewidth',3,'CapSize',0)
            plot(median(spkPerStim.byStim.(fldLabels{ii}),1),'linewidth',3,'color',inputs.plotColorGrpTwo)
        end
    end
end
spsAx.byStim = gca;
setAx(spsAx.byStim);
xlim([.5 size(spkPerStim.byStim.(fldLabels{1}),2)+.5])
xlabel('Stimulus Number')
spsAx.byStim.XTick = 1:1:size(spkPerStim.byStim.(fldLabels{1}),2);
if spsAx.byStim.YLim(2) < 1
    spsAx.byStim.YLim(2) = 1;
end
spsAx.byStim.YLim(1) = 0;
ylabel('Spikes per stimulus')

%spikes per stimulus, by sweep
spkPerStimFig.bySweep = figure(noGroups+2);
spkPerStimFig.bySweep.Position = [635 450 100 275];
hold on
if noGroups == 2
    if inputs.PairedAnalysisCheckBox.Value == 1
        for ii = 1:length(spkPerStim.bySweep.(fldLabels{1}))
            plot([1 2],[spkPerStim.bySweep.(fldLabels{1})(ii),spkPerStim.bySweep.(fldLabels{2})(ii)],'color',lightgray,'linewidth',1)
        end
        if hSPS == 0 %if data normal, connect means
            plot([1 2],[mean(spkPerStim.bySweep.(fldLabels{1})) mean(spkPerStim.bySweep.(fldLabels{2}))],'linewidth',2,'color',darkgray)
        elseif hSPS ~= 0 %if data not normal, connect medians
            plot([1 2],[median(spkPerStim.bySweep.(fldLabels{1})) median(spkPerStim.bySweep.(fldLabels{2}))],'linewidth',2,'color',darkgray)
        end
    end
end
for ii = 1:noGroups
    scatter(ii.*ones(1,length(spkPerStim.bySweep.(fldLabels{ii}))),spkPerStim.bySweep.(fldLabels{ii}),40,lightgray,'filled')
    if hSPS == 0 %if data normal, plot mean + sem
        if ii == 1
            errorbar(ii,mean(spkPerStim.bySweep.(fldLabels{ii})),sem(spkPerStim.bySweep.(fldLabels{ii}),2),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,mean(spkPerStim.bySweep.(fldLabels{ii})),125,inputs.plotColorGrpOne,'filled')
        elseif ii == 2
            errorbar(ii,mean(spkPerStim.bySweep.(fldLabels{ii})),sem(spkPerStim.bySweep.(fldLabels{ii}),2),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,mean(spkPerStim.bySweep.(fldLabels{ii})),125,inputs.plotColorGrpTwo,'filled')
        end
    elseif hSPS ~=0 %if data not normal, plot median + iqr
        [lowerError.spkPerStim.bySweep.(fldLabels{ii}),upperError.spkPerStim.bySweep.(fldLabels{ii})] = iqrError(spkPerStim.bySweep.(fldLabels{ii}),2);
        if ii == 1
            errorbar(ii,median(spkPerStim.bySweep.(fldLabels{ii})),lowerError.spkPerStim.bySweep.(fldLabels{ii}),upperError.spkPerStim.bySweep.(fldLabels{ii}),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,median(spkPerStim.bySweep.(fldLabels{ii})),125,inputs.plotColorGrpOne,'filled')
        elseif ii == 2
            errorbar(ii,median(spkPerStim.bySweep.(fldLabels{ii})),lowerError.spkPerStim.bySweep.(fldLabels{ii}),upperError.spkPerStim.bySweep.(fldLabels{ii}),'linewidth',3,'color',inputs.plotColorGrpTwo,'capsize',0)
            scatter(ii,median(spkPerStim.bySweep.(fldLabels{ii})),125,inputs.plotColorGrpTwo,'filled')
        end
    end
end
spsAx.bySweep = gca;
setAx(spsAx.bySweep);
xlim([.5 noGroups+.5])
for ii = 1:noGroups
    spsAx.bySweep.XTickLabel{ii} = fldLabels{ii};
end
if spsAx.bySweep.YLim(2) ~= spsAx.byStim.YLim(2)
    if spsAx.bySweep.YLim(2) > spsAx.byStim.YLim(2)
        spsAx.byStim.YLim(2) =  spsAx.bySweep.YLim(2);
    elseif spsAx.bySweep.YLim(2) < spsAx.byStim.YLim(2)
        spsAx.bySweep.YLim(2) = spsAx.byStim.YLim(2);
    end
end
spsAx.bySweep.YLim(1) = 0;
theseTicks = 0:.25:spsAx.bySweep.YLim(2);
spsAx.bySweep.YTick = theseTicks;
spsAx.byStim.YTick = theseTicks;
ylabel('Spikes per stimulus')

%membrane voltage
vmFig = figure(noGroups+3);
vmFig.Position = [745 450 100 275];
hold on
if noGroups == 2
    if inputs.PairedAnalysisCheckBox.Value == 1
        for ii = 1:length(Vm.(fldLabels{1}))
            plot([1 2],[Vm.(fldLabels{1})(ii),Vm.(fldLabels{2})(ii)],'color',lightgray,'linewidth',1)
        end
        if hVm == 0 %if data normal, connect means
            plot([1 2],[mean(Vm.(fldLabels{1})) mean(Vm.(fldLabels{2}))],'linewidth',2,'color',darkgray)
        elseif hVm ~= 0 %if data not normal, connect medians
            plot([1 2],[median(Vm.(fldLabels{1})) median(Vm.(fldLabels{2}))],'linewidth',2,'color',darkgray)
        end
    end
end
for ii = 1:noGroups
    scatter(ii.*ones(1,length(Vm.(fldLabels{ii}))),Vm.(fldLabels{ii}),40,lightgray,'filled')
    if hVm == 0 %if data normal, plot mean + sem
        if ii == 1
            errorbar(ii,mean(Vm.(fldLabels{ii})),sem(Vm.(fldLabels{ii}),2),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,mean(Vm.(fldLabels{ii})),125,inputs.plotColorGrpOne,'filled')
        elseif ii == 2
            errorbar(ii,mean(Vm.(fldLabels{ii})),sem(Vm.(fldLabels{ii}),2),'linewidth',3,'color',inputs.plotColorGrpTwo,'capsize',0)
            scatter(ii,mean(Vm.(fldLabels{ii})),125,inputs.plotColorGrpTwo,'filled')
        end
    elseif hVm ~=0 %if data not normal, plot median + iqr
        [lowerError.Vm.(fldLabels{ii}),upperError.Vm.(fldLabels{ii})] = iqrError(Vm.(fldLabels{ii}),2);
        if ii == 1
            errorbar(ii,median(Vm.(fldLabels{ii})),lowerError.Vm.(fldLabels{ii}),upperError.Vm.(fldLabels{ii}),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,median(Vm.(fldLabels{ii})),125,inputs.plotColorGrpOne,'filled')
        elseif ii == 2
            errorbar(ii,median(Vm.(fldLabels{ii})),lowerError.Vm.(fldLabels{ii}),upperError.Vm.(fldLabels{ii}),'linewidth',3,'color',inputs.plotColorGrpTwo,'capsize',0)
            scatter(ii,median(Vm.(fldLabels{ii})),125,inputs.plotColorGrpTwo,'filled')
        end
    end
end
vmAx = gca;
setAx(vmAx);
xlim([.5 noGroups+.5])
for ii = 1:noGroups
    vmAx.XTickLabel{ii} = fldLabels{ii};
end
lowY = floor(vmAx.YLim(1)/10);
ylim([10*lowY -30])
theseTicks = [];
theseTicks = 10*lowY:10:-30;
vmAx.YTick = theseTicks;
ylabel('Vm (mV)')

%stEPSP
%plot cumulative distribution and cell subthreshold EPSPs
stEPSPFig=figure(noGroups+4);
stEPSPFig.Position=[855 450 475 275];
hold on
for ii = 1:noGroups
    if ii == 1
        plot(xEPSP.(fldLabels{ii}),fEPSP.(fldLabels{ii}),'linewidth',2,'color',inputs.plotColorGrpOne)
    elseif ii == 2
        plot(xEPSP.(fldLabels{ii}),fEPSP.(fldLabels{ii}),'linewidth',2,'color',inputs.plotColorGrpTwo)
    end
end
cdEPSPAx = gca;
setAx(cdEPSPAx);
xlabel('stEPSP amplitude (mV)')
theseTicks = [];
theseTicks = 0:cdEPSPAx.XTick(end)/4:cdEPSPAx.XTick(end);
cdEPSPAx.XTick = theseTicks;
ylabel('P(stEPSP)')
title('stEPSP data')

%plot cellular mean stEPSP
axes('Position',[.7 .2 .125 .4])
cmEPSPAx = gca;
setAx(cmEPSPAx);
hold on
if noGroups == 2
    if inputs.PairedAnalysisCheckBox.Value == 1
        for ii = 1:length(cellMeanstEPSP.(fldLabels{1}))
            plot([1 2],[cellMeanstEPSP.(fldLabels{1})(ii),cellMeanstEPSP.(fldLabels{2})(ii)],'color',lightgray,'linewidth',1)
        end
        if hstEPSP == 0 %if data normal, connect means
            plot([1 2],[mean(cellMeanstEPSP.(fldLabels{1})) mean(cellMeanstEPSP.(fldLabels{2}))],'linewidth',2,'color',darkgray)
        elseif hstEPSP ~= 0 %if data not normal, connect medians
            plot([1 2],[median(cellMeanstEPSP.(fldLabels{1})) median(cellMeanstEPSP.(fldLabels{2}))],'linewidth',2,'color',darkgray)
        end
    end
end
for ii = 1:noGroups
    scatter(ii.*ones(1,length(cellMeanstEPSP.(fldLabels{ii}))),cellMeanstEPSP.(fldLabels{ii}),40,lightgray,'filled')
    if hstEPSP == 0 %if data normal, plot mean + sem
        if ii == 1
            errorbar(ii,mean(cellMeanstEPSP.(fldLabels{ii})),sem(cellMeanstEPSP.(fldLabels{ii}),2),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,mean(cellMeanstEPSP.(fldLabels{ii})),125,inputs.plotColorGrpOne,'filled')
        elseif ii == 2
            errorbar(ii,mean(cellMeanstEPSP.(fldLabels{ii})),sem(cellMeanstEPSP.(fldLabels{ii}),2),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,mean(cellMeanstEPSP.(fldLabels{ii})),125,inputs.plotColorGrpTwo,'filled')
        end
    elseif hstEPSP ~=0 %if data not normal, plot median + iqr
        [lowerError.cellMeanstEPSP.(fldLabels{ii}),upperError.cellMeanstEPSP.(fldLabels{ii})] = iqrError(cellMeanstEPSP.(fldLabels{ii}),2);
        if ii == 1
            errorbar(ii,median(cellMeanstEPSP.(fldLabels{ii})),lowerError.cellMeanstEPSP.(fldLabels{ii}),upperError.cellMeanstEPSP.(fldLabels{ii}),'linewidth',3,'color',inputs.plotColorGrpOne,'capsize',0)
            scatter(ii,median(cellMeanstEPSP.(fldLabels{ii})),125,inputs.plotColorGrpOne,'filled')
        elseif ii == 2
            errorbar(ii,median(cellMeanstEPSP.(fldLabels{ii})),lowerError.cellMeanstEPSP.(fldLabels{ii}),upperError.cellMeanstEPSP.(fldLabels{ii}),'linewidth',3,'color',inputs.plotColorGrpTwo,'capsize',0)
            scatter(ii,median(cellMeanstEPSP.(fldLabels{ii})),125,inputs.plotColorGrpTwo,'filled')
        end
    end
end
xlim([.5 noGroups+.5])
for ii = 1:noGroups
    cmEPSPAx.XTickLabel{ii} = fldLabels{ii};
end

%plot example stEPSP traces
for ii = 1:noGroups
    stEPSPTraceFig.(fldLabels{ii}) = figure(noGroups+4+ii);
    stEPSPTraceFig.(fldLabels{ii}).Position = [945 275-185*(ii-1) 400 100];
    if ii == 1
        plot(t.stEPSP,mstTrace.(fldLabels{ii}),'linewidth',2,'color',inputs.plotColorGrpOne)
    elseif ii == 2
        plot(t.stEPSP,mstTrace.(fldLabels{ii}),'linewidth',2,'color',inputs.plotColorGrpTwo)
    end
    stEPSPAx.(fldLabels{ii}) = gca;
    setAx(stEPSPAx.(fldLabels{ii}));
    xlabel('time (ms)')
    xlim(1000.*[1/samplerate length(mstTrace.(fldLabels{ii}))/samplerate])
    ylabel('stEPSP amplitude (mV)')
end
if noGroups == 2
    if stEPSPAx.(fldLabels{1}).YLim(2) > stEPSPAx.(fldLabels{2}).YLim(2)
        stEPSPAx.(fldLabels{2}).YLim(2) = stEPSPAx.(fldLabels{1}).YLim(2);
        stEPSPAx.(fldLabels{1}).YLim(1) = -1;
        stEPSPAx.(fldLabels{2}).YLim(1) = -1;
    elseif stEPSPAx.(fldLabels{1}).YLim(2) < stEPSPAx.(fldLabels{2}).YLim(2)
        stEPSPAx.(fldLabels{1}).YLim(2) = stEPSPAx.(fldLabels{2}).YLim(2);
        stEPSPAx.(fldLabels{1}).YLim(1) = -1;
        stEPSPAx.(fldLabels{2}).YLim(1) = -1;
    end
end

output.p
end
