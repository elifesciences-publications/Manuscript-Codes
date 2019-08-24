%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% spon event Grp Comp %%%%%%%%%%%%
%%%%%%%%%%% Created: 10-20-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 02-15-2019 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outputdata] = seCompFxn(inputs)
close all

%% inits
%% INIT VARS
fldrName = inputs.FolderNameEditField.Value;
noGroups = sum([inputs.ExpgrponeCheckBox.Value inputs.ExpgrptwoCheckBox.Value inputs.ExpgrpthreeCheckBox.Value]);
for ii = 1:noGroups
    if ii == 1
        grpLabels{ii} = inputs.Grp1EditField.Value;
    elseif ii == 2
        grpLabels{ii} = inputs.Grp2EditField.Value;
    elseif ii == 3
        grpLabels{ii} = inputs.Grp3EditField.Value;
    end
end
for ii = 1:noGroups
    fldLabels{ii} = grpLabels{ii};
    dashes = find(fldLabels{ii}=='-');
    fldLabels{ii}(dashes) = [];
    clear dashes
end

samplerate = inputs.SamplerateHzEditField_2.Value;
lightgray=[.75 .75 .75]; %light gray for individual data points

%% LOAD DATA
% Directory
cd(inputs.spscDir);

for ii = 1:noGroups
    %get files
    cd([inputs.spscDir,'/',grpLabels{ii},'/',fldrName]);
    contents.(fldLabels{ii}) = dir('*.mat');
    filenames.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
    spscFiles.(fldLabels{ii}) = fullfile(cd,filenames.(fldLabels{ii}));
    
    for jj = 1:length(spscFiles.(fldLabels{ii}))
        load(spscFiles.(fldLabels{ii}){jj}, 'app')
        
        outputdata.amplitude.(fldLabels{ii})(jj) = abs(app.spscData.mAmp);
        outputdata.frequency.(fldLabels{ii})(jj) = app.spscData.mFreq;
        if strcmp(fldrName,'sEPSC')
            outputdata.rise.(fldLabels{ii})(jj) = app.spscData.sEPSCrt;
            outputdata.decay.(fldLabels{ii})(jj) = app.spscData.sEPSCtau;
        elseif strcmp(fldrName,'sIPSC')
            outputdata.rise.(fldLabels{ii})(jj) = app.spscData.sIPSCrt;
            outputdata.decay.(fldLabels{ii})(jj) = app.spscData.sIPSCtau;
        end
        
        if ii == 1
            exCell.(fldLabels{ii}) = inputs.ExCellEditField.Value;
        elseif ii == 2
            exCell.(fldLabels{ii}) = inputs.ExCellEditField_2.Value;
        elseif ii == 3
            exCell.(fldLabels{ii}) = inputs.ExCellEditField_3.Value;
        end
        if jj == exCell.(fldLabels{ii})
            outputdata.eventTraces.(fldLabels{ii}) = app.spscData.alignedEvents;
            outputdata.ieiDist.(fldLabels{ii}) = app.spscData.eventiei;
            outputdata.frequencyDist.(fldLabels{ii}) = app.spscData.eventFreq;
            outputdata.rawTrace.(fldLabels{ii}) = app.spscData.sponData;
        end
        clear app
    end
end

%time vector
t = samplerate^-1:samplerate^-1:length(outputdata.rawTrace.(fldLabels{ii}))/samplerate;

%% ANALYSIS
%amplitude, normality
if noGroups == 3
    for ii = 1:noGroups
        res.amp.(fldLabels{ii}) = mean(outputdata.amplitude.(fldLabels{ii})) - outputdata.amplitude.(fldLabels{ii});
    end
    allRes.amp = [res.amp.(fldLabels{1})'; res.amp.(fldLabels{2})'; res.amp.(fldLabels{3})'];
    [hAmp,~]=adtest(allRes.amp);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(outputdata.amplitude.(fldLabels{ii}));
    end
    if noGroups == 1
        hAmp = holdH.(fldLabels{1});
    elseif noGroups == 2
        hAmp = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%amplitude, comparison tests
if noGroups == 3
    if hAmp == 0 %if residuals normal, run anova
        [outputdata.p.ampAnova,~,ampstats]=anova1([outputdata.amplitude.(fldLabels{1}) outputdata.amplitude.(fldLabels{2}) outputdata.amplitude.(fldLabels{3})],...
            [ones(1,length(outputdata.amplitude.(fldLabels{1}))) 2.*ones(1,length(outputdata.amplitude.(fldLabels{2}))) 3.*ones(1,length(outputdata.amplitude.(fldLabels{3})))]);
        if outputdata.p.ampAnova <= 0.05
            outputdata.ampPostHocTukey=multcompare(ampstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [outputdata.p.ampKW,~,~]=kruskalwallis([outputdata.amplitude.(fldLabels{1}) outputdata.amplitude.(fldLabels{2}) outputdata.amplitude.(fldLabels{3})],...
            [ones(1,length(outputdata.amplitude.(fldLabels{1}))) 2.*ones(1,length(outputdata.amplitude.(fldLabels{2}))) 3.*ones(1,length(outputdata.amplitude.(fldLabels{3})))]);
        if outputdata.p.ampKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [outputdata.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.amplitude.(fldLabels{1}),outputdata.amplitude.(fldLabels{2}));
            [outputdata.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.amplitude.(fldLabels{1}),outputdata.amplitude.(fldLabels{3}));
            [outputdata.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.amplitude.(fldLabels{2}),outputdata.amplitude.(fldLabels{3}));
            outputdata.pFDR.amp=drsFDRpval([outputdata.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hAmp == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,outputdata.p.ampTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.amplitude.(fldLabels{1}),outputdata.amplitude.(fldLabels{2}));
        end
    elseif hAmp ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [outputdata.p.ampMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.amplitude.(fldLabels{1}),outputdata.amplitude.(fldLabels{2}));
        end
    end
end

%frequency, normality
if noGroups == 3
    for ii = 1:noGroups
        res.freq.(fldLabels{ii}) = mean(outputdata.frequency.(fldLabels{ii})) - outputdata.frequency.(fldLabels{ii});
    end
    allRes.freq = [res.freq.(fldLabels{1})'; res.freq.(fldLabels{2})'; res.freq.(fldLabels{3})'];
    [hFreq,~]=adtest(allRes.freq);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(outputdata.frequency.(fldLabels{ii}));
    end
    if noGroups == 1
        hFreq = holdH.(fldLabels{1});
    elseif noGroups == 2
        hFreq = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%frequency, comparison tests
if noGroups == 3
    if hFreq == 0 %if residuals normal, run anova
        [outputdata.p.freqAnova,~,freqstats]=anova1([outputdata.frequency.(fldLabels{1}) outputdata.frequency.(fldLabels{2}) outputdata.frequency.(fldLabels{3})],...
            [ones(1,length(outputdata.frequency.(fldLabels{1}))) 2.*ones(1,length(outputdata.frequency.(fldLabels{2}))) 3.*ones(1,length(outputdata.frequency.(fldLabels{3})))]);
        if outputdata.p.freqAnova <= 0.05
            outputdata.freqPostHocTukey=multcompare(freqstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [outputdata.p.freqKW,~,~]=kruskalwallis([outputdata.frequency.(fldLabels{1}) outputdata.frequency.(fldLabels{2}) outputdata.frequency.(fldLabels{3})],...
            [ones(1,length(outputdata.frequency.(fldLabels{1}))) 2.*ones(1,length(outputdata.frequency.(fldLabels{2}))) 3.*ones(1,length(outputdata.frequency.(fldLabels{3})))]);
        if outputdata.p.freqKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [outputdata.p.freqPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.frequency.(fldLabels{1}),outputdata.frequency.(fldLabels{2}));
            [outputdata.p.freqPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.frequency.(fldLabels{1}),outputdata.frequency.(fldLabels{3}));
            [outputdata.p.freqPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.frequency.(fldLabels{2}),outputdata.frequency.(fldLabels{3}));
            outputdata.pFDR.freq=drsFDRpval([outputdata.p.freqPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.freqPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.freqPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hFreq == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,outputdata.p.freqTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.frequency.(fldLabels{1}),outputdata.frequency.(fldLabels{2}));
        end
    elseif hFreq ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [outputdata.p.freqMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.frequency.(fldLabels{1}),outputdata.frequency.(fldLabels{2}));
        end
    end
end

%rise time, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.rise.(fldLabels{ii}) = mean(outputdata.rise.(fldLabels{ii})) - outputdata.rise.(fldLabels{ii});
    end
    allRes.rise = [res.rise.(fldLabels{1})'; res.rise.(fldLabels{2})'; res.rise.(fldLabels{3})'];
    [hRT,~]=adtest(allRes.rise);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(rise.(fldLabels{ii}));
    end
    if noGroups == 1
        hRT = holdH.(fldLabels{1});
    elseif noGroups == 2
        hRT = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%rise time, comparison tests
if noGroups == 3
    if hRT == 0 %if residuals normal, run anova
        [outputdata.p.riseAnova,~,risestats]=anova1([outputdata.rise.(fldLabels{1}) outputdata.rise.(fldLabels{2}) outputdata.rise.(fldLabels{3})],...
            [ones(1,length(outputdata.rise.(fldLabels{1}))) 2.*ones(1,length(outputdata.rise.(fldLabels{2}))) 3.*ones(1,length(outputdata.rise.(fldLabels{3})))]);
        if outputdata.p.riseAnova <= 0.05
            outputdata.risePostHocTukey=multcompare(risestats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [outputdata.p.riseKW,~,~]=kruskalwallis([outputdata.rise.(fldLabels{1}) outputdata.rise.(fldLabels{2}) outputdata.rise.(fldLabels{3})],...
            [ones(1,length(outputdata.rise.(fldLabels{1}))) 2.*ones(1,length(outputdata.rise.(fldLabels{2}))) 3.*ones(1,length(outputdata.rise.(fldLabels{3})))]);
        if outputdata.p.riseKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [outputdata.p.risePostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.rise.(fldLabels{1}),outputdata.rise.(fldLabels{2}));
            [outputdata.p.risePostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.rise.(fldLabels{1}),outputdata.rise.(fldLabels{3}));
            [outputdata.p.risePostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.rise.(fldLabels{2}),outputdata.rise.(fldLabels{3}));
            outputdata.pFDR.rise=drsFDRpval([outputdata.p.risePostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.risePostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.risePostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hRT == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,outputdata.p.riseTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.rise.(fldLabels{1}),outputdata.rise.(fldLabels{2}));
        end
    elseif hRT ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [outputdata.p.riseMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.rise.(fldLabels{1}),outputdata.rise.(fldLabels{2}));
        end
    end
end

%decay tau, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.decay.(fldLabels{ii}) = mean(outputdata.decay.(fldLabels{ii})) - outputdata.decay.(fldLabels{ii});
    end
    allRes.decay = [res.decay.(fldLabels{1})'; res.decay.(fldLabels{2})'; res.decay.(fldLabels{3})'];
    [hDecay,~]=adtest(allRes.decay);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(decay.(fldLabels{ii}));
    end
    if noGroups == 1
        hDecay = holdH.(fldLabels{1});
    elseif noGroups == 2
        hDecay = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%decay tau, comparison tests
if noGroups == 3
    if hDecay == 0 %if residuals normal, run anova
        [outputdata.p.decayAnova,~,decaystats]=anova1([outputdata.decay.(fldLabels{1}) outputdata.decay.(fldLabels{2}) outputdata.decay.(fldLabels{3})],...
            [ones(1,length(outputdata.decay.(fldLabels{1}))) 2.*ones(1,length(outputdata.decay.(fldLabels{2}))) 3.*ones(1,length(outputdata.decay.(fldLabels{3})))]);
        if outputdata.p.decayAnova <= 0.05
            outputdata.decayPostHocTukey=multcompare(decaystats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [outputdata.p.decayKW,~,~]=kruskalwallis([outputdata.decay.(fldLabels{1}) outputdata.decay.(fldLabels{2}) outputdata.decay.(fldLabels{3})],...
            [ones(1,length(outputdata.decay.(fldLabels{1}))) 2.*ones(1,length(outputdata.decay.(fldLabels{2}))) 3.*ones(1,length(outputdata.decay.(fldLabels{3})))]);
        if outputdata.p.decayKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [outputdata.p.decayPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.decay.(fldLabels{1}),outputdata.decay.(fldLabels{2}));
            [outputdata.p.decayPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.decay.(fldLabels{1}),outputdata.decay.(fldLabels{3}));
            [outputdata.p.decayPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.decay.(fldLabels{2}),outputdata.decay.(fldLabels{3}));
            outputdata.pFDR.decay=drsFDRpval([outputdata.p.decayPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.decayPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.decayPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hDecay == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,outputdata.p.decayTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.decay.(fldLabels{1}),outputdata.decay.(fldLabels{2}));
        end
    elseif hDecay ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [outputdata.p.decayMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.decay.(fldLabels{1}),outputdata.decay.(fldLabels{2}));
        end
    end
end

%% DISPLAY DATA
close all

%trace figures
for ii = 1:noGroups
    mTraceFig.(fldLabels{ii}) = figure(ii);
    mTraceFig.(fldLabels{ii}).Position = [150+(ii-1)*425 500 425 150];
    plot(t,outputdata.rawTrace.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',1)
    traceAx.(fldLabels{ii}) = gca;
    setAx(traceAx.(fldLabels{ii}));
    ylabel('current (pA)')
    xlabel('t (s)')
    title([fldLabels{ii} ' spontaneous currents'])
%     xTemp = randperm(length(outputdata.rawTrace.(fldLabels{ii}))-.5*samplerate,1);
%     xlim([t(xTemp) t(xTemp+.5*samplerate)])
    if ii == 1
        yRangeHolder = traceAx.(fldLabels{ii}).YLim(2)-traceAx.(fldLabels{ii}).YLim(1);
    else
        yRangeHolder = [yRangeHolder traceAx.(fldLabels{ii}).YLim(2)-traceAx.(fldLabels{ii}).YLim(1)];
    end
end
yRangeMax = max(yRangeHolder);
for ii = 1:noGroups
    traceAx.(fldLabels{ii}).YLim = [traceAx.(fldLabels{ii}).YLim(2)-yRangeMax traceAx.(fldLabels{ii}).YLim(2)];
end

%amplitude figure
ampFig=figure(noGroups+1);
ampFig.Position = [140 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(outputdata.amplitude.(fldLabels{ii}))),outputdata.amplitude.(fldLabels{ii}),40,lightgray,'filled')
    if hAmp == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(outputdata.amplitude.(fldLabels{ii})),sem(outputdata.amplitude.(fldLabels{ii}),2),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(outputdata.amplitude.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hAmp ~= 0 %if data non-normal, plot median w IQR
        [lowerError.outputdata.amplitude.(fldLabels{ii}),upperError.outputdata.amplitude.(fldLabels{ii})] = iqrError(outputdata.amplitude.(fldLabels{ii}),1);
        errorbar(ii,median(outputdata.amplitude.(fldLabels{ii})),lowerError.outputdata.amplitude.(fldLabels{ii}),upperError.outputdata.amplitude.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(outputdata.amplitude.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allAmp = outputdata.amplitude.(fldLabels{ii})';
    else
        allAmp = [allAmp; outputdata.amplitude.(fldLabels{ii})'];
    end
end
ylabel('event amplitude (pA)')
ampAx=gca;
setAx(ampAx); 
ampAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 1.2*max(allAmp)])

%frequency figure
freqFig=figure(noGroups+2);
freqFig.Position = [290 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(outputdata.frequency.(fldLabels{ii}))),outputdata.frequency.(fldLabels{ii}),40,lightgray,'filled')
    if hFreq == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(outputdata.frequency.(fldLabels{ii})),sem(outputdata.frequency.(fldLabels{ii}),2),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(outputdata.frequency.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hFreq ~= 0 %if data non-normal, plot median w IQR
        [lowerError.outputdata.frequency.(fldLabels{ii}),upperError.outputdata.frequency.(fldLabels{ii})] = iqrError(outputdata.frequency.(fldLabels{ii}),1);
        errorbar(ii,median(outputdata.frequency.(fldLabels{ii})),lowerError.outputdata.frequency.(fldLabels{ii}),upperError.outputdata.frequency.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(outputdata.frequency.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allFreq = outputdata.frequency.(fldLabels{ii})';
    else
        allFreq = [allFreq; outputdata.frequency.(fldLabels{ii})'];
    end
end
ylabel('event frequency (Hz)')
ampAx=gca;
setAx(ampAx); 
ampAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 1.2*max(allFreq)])

%freq vs amplitude plot
faFig = figure(noGroups+3);
faFig.Position = [440 100 325 300];
hold on
for ii = 1:noGroups
scatter(outputdata.amplitude.(fldLabels{ii}),...
    outputdata.frequency.(fldLabels{ii}),75,inputs.plotColorGrps(ii,:),'filled')
end
faAx = gca;
setAx(faAx);
xlabel('amplitude (pA)')
ylabel('frequency (Hz)')
xlim([0 1.2*max(allAmp)])
ylim([0 1.2*max(allFreq)])

%rise figure
riseFig=figure(noGroups+4);
riseFig.Position = [780 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(outputdata.rise.(fldLabels{ii}))),outputdata.rise.(fldLabels{ii}),40,lightgray,'filled')
    if hRT == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(outputdata.rise.(fldLabels{ii})),sem(outputdata.rise.(fldLabels{ii}),2),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(outputdata.rise.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hRT ~= 0 %if data non-normal, plot median w IQR
        [lowerError.outputdata.rise.(fldLabels{ii}),upperError.outputdata.rise.(fldLabels{ii})] = iqrError(outputdata.rise.(fldLabels{ii}),2);
        errorbar(ii,median(outputdata.rise.(fldLabels{ii})),lowerError.outputdata.rise.(fldLabels{ii}),upperError.outputdata.rise.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(outputdata.rise.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allRise = outputdata.rise.(fldLabels{ii})';
    else
        allRise = [allRise; outputdata.rise.(fldLabels{ii})'];
    end
end
ylabel('20-80% risetime (ms)')
ampAx=gca;
setAx(ampAx); 
ampAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 1.2*max(allRise)])

%decay figure
decayFig=figure(noGroups+5);
decayFig.Position = [920 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(outputdata.decay.(fldLabels{ii}))),outputdata.decay.(fldLabels{ii}),40,lightgray,'filled')
    if hDecay == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(outputdata.decay.(fldLabels{ii})),sem(outputdata.decay.(fldLabels{ii}),2),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(outputdata.decay.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hDecay ~= 0 %if data non-normal, plot median w IQR
        [lowerError.outputdata.decay.(fldLabels{ii}),upperError.outputdata.decay.(fldLabels{ii})] = iqrError(outputdata.decay.(fldLabels{ii}),2);
        errorbar(ii,median(outputdata.decay.(fldLabels{ii})),lowerError.outputdata.decay.(fldLabels{ii}),upperError.outputdata.decay.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(outputdata.decay.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allTau = outputdata.decay.(fldLabels{ii})';
    else
        allTau = [allTau; outputdata.decay.(fldLabels{ii})'];
    end
end
ylabel('decay \tau (ms)')
ampAx=gca;
setAx(ampAx); 
ampAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 1.2*max(allTau)])

%% REPORT RESULTS
outputdata
outputdata.p

end
