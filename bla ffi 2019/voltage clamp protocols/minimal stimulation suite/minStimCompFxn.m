function [output] = minStimCompFxn(inputs)
%for group comparison of minimal stimulation data
%modeled off ms_dispstats.m code
%created 05-03-2018, edited 05-04-18

%% INIT VARS
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

%% LOAD DATA
% Directory
cd(inputs.minStimDir);

%load responding cell proportions
contents.excel = dir('MinStimRespPercentages.xlsx');
filename.excel = {contents.excel.name}';
excelFile = fullfile(cd,filename.excel);

for ii = 1:noGroups
    %requires file to have the proportion of responding cells in E2 and for sheets to be in order of groups
    proportionResp.(fldLabels{ii}) = xlsread(excelFile{1},ii,'E2'); 
end

for ii = 1:noGroups
    %get files
    %ctrl data
    if strcmp(grpLabels{ii},'PrincipalNeurons')
        cd([inputs.minStimDir,'/Wildtype/minstim/ctrl']);
    else
        cd([inputs.minStimDir,'/',grpLabels{ii},'/minstim/ctrl']);
    end
    contents.(fldLabels{ii}) = dir('*.mat');
    filenames.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
    minstimFiles.(fldLabels{ii}).ctrl = fullfile(cd,filenames.(fldLabels{ii}));
    
    %inits
    amp.(fldLabels{ii})=-1.*ones(length(minstimFiles.(fldLabels{ii}).ctrl),1);
    lat.(fldLabels{ii})=-1.*ones(length(minstimFiles.(fldLabels{ii}).ctrl),1);
    jtr.(fldLabels{ii})=-1.*ones(length(minstimFiles.(fldLabels{ii}).ctrl),1);
    rt.(fldLabels{ii})=-1.*ones(length(minstimFiles.(fldLabels{ii}).ctrl),1);
    dtau.(fldLabels{ii})=-1.*ones(length(minstimFiles.(fldLabels{ii}).ctrl),1);
    stim.(fldLabels{ii})=-1.*ones(length(minstimFiles.(fldLabels{ii}).ctrl),1);
    success.(fldLabels{ii}).ctrl=-1.*ones(length(minstimFiles.(fldLabels{ii}).ctrl),1);
    
    for jj = 1:length(minstimFiles.(fldLabels{ii}).ctrl)
        %load ctrl data
        load(minstimFiles.(fldLabels{ii}).ctrl{jj}, 'app')
        amp.(fldLabels{ii})(jj) = app.minStimData.mAmp;
        lat.(fldLabels{ii})(jj) = app.minStimData.mLat;
        jtr.(fldLabels{ii})(jj) = app.minStimData.jitter;
        rt.(fldLabels{ii})(jj) = app.minStimData.mRT;
        dtau.(fldLabels{ii})(jj) = app.minStimData.decay;
        stim.(fldLabels{ii})(jj) = app.minStimData.stimVal;
        success.(fldLabels{ii}).ctrl(jj) = app.minStimData.sRate;
        
        if ii == 1
            exCell.(fldLabels{ii}) = inputs.ExCellEditField.Value;
        elseif ii == 2
            exCell.(fldLabels{ii}) = inputs.ExCellEditField_2.Value;
        elseif ii == 3
            exCell.(fldLabels{ii}) = inputs.ExCellEditField_3.Value;
        end
        if jj == exCell.(fldLabels{ii})
            meanTrace.(fldLabels{ii}) = app.minStimData.mTrace;
            tracesAligned.(fldLabels{ii}) = app.minStimData.alignedTraces;
        end
        clear app
    end
    
    
    %dnqx/apv data
    if strcmp(grpLabels{ii},'PrincipalNeurons')
        cd([inputs.minStimDir,'/Wildtype/minstim/xv']);
    else
        cd([inputs.minStimDir,'/',grpLabels{ii},'/minstim/xv']);
    end
    contents.(fldLabels{ii}) = dir('*.mat');
    filenames.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
    minstimFiles.(fldLabels{ii}).xv = fullfile(cd,filenames.(fldLabels{ii}));
    
    %init
    success.(fldLabels{ii}).xv=-1.*ones(length(minstimFiles.(fldLabels{ii}).xv),1);
    
    %load xv data
    for jj = 1:length(minstimFiles.(fldLabels{ii}).xv)
        load(minstimFiles.(fldLabels{ii}).xv{jj},'app')
        success.(fldLabels{ii}).xv(jj) = app.minStimData.sRate;
    end
end

%% TEST DIFFERENCES
%adtest (Anderson-Darling Test) to test for normality
%amplitude, normality
if noGroups == 3
    for ii = 1:noGroups
        res.amp.(fldLabels{ii}) = mean(amp.(fldLabels{ii})) - amp.(fldLabels{ii});
    end
    allRes.amp = [res.amp.(fldLabels{1}); res.amp.(fldLabels{2}); res.amp.(fldLabels{3})];
    [hAmp,~]=adtest(allRes.amp);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(amp.(fldLabels{ii}));
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
        [output.p.ampAnova,~,ampstats]=anova1([amp.(fldLabels{1})' amp.(fldLabels{2})' amp.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.ampAnova <= 0.05
            output.ampPostHocTukey=multcompare(ampstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.ampKW,~,~]=kruskalwallis([amp.(fldLabels{1})' amp.(fldLabels{2})' amp.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.ampKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(amp.(fldLabels{1}),amp.(fldLabels{2}));
            [output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(amp.(fldLabels{1}),amp.(fldLabels{3}));
            [output.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(amp.(fldLabels{2}),amp.(fldLabels{3}));
            output.pFDR.amp=drsFDRpval([output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hAmp == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,output.p.ampTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(amp.(fldLabels{1}),amp.(fldLabels{2}));
        end
    elseif hAmp ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [output.p.ampMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(amp.(fldLabels{1}),amp.(fldLabels{2}));
        end
    end
end

%latency, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.lat.(fldLabels{ii}) = mean(lat.(fldLabels{ii})) - lat.(fldLabels{ii});
    end
    allRes.lat = [res.lat.(fldLabels{1}); res.lat.(fldLabels{2}); res.lat.(fldLabels{3})];
    [hLat,~]=adtest(allRes.lat);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(lat.(fldLabels{ii}));
    end
    if noGroups == 1
        hLat = holdH.(fldLabels{1});
    elseif noGroups == 2
        hLat = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%latency, comparison tests
if noGroups == 3
    if hLat == 0 %if residuals normal, run anova
        [output.p.latAnova,~,latstats]=anova1([lat.(fldLabels{1})' lat.(fldLabels{2})' lat.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.latAnova <= 0.05
            output.latPostHocTukey=multcompare(latstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.latKW,~,~]=kruskalwallis([lat.(fldLabels{1})' lat.(fldLabels{2})' lat.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.latKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.latPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(lat.(fldLabels{1}),lat.(fldLabels{2}));
            [output.p.latPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(lat.(fldLabels{1}),lat.(fldLabels{3}));
            [output.p.latPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(lat.(fldLabels{2}),lat.(fldLabels{3}));
            output.pFDR.lat=drsFDRpval([output.p.latPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.latPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.latPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hLat == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,output.p.latTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(lat.(fldLabels{1}),lat.(fldLabels{2}));
        end
    elseif hLat ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [output.p.latMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(lat.(fldLabels{1}),lat.(fldLabels{2}));
        end
    end
end

%jitter, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.jtr.(fldLabels{ii}) = mean(jtr.(fldLabels{ii})) - jtr.(fldLabels{ii});
    end
    allRes.jtr = [res.jtr.(fldLabels{1}); res.jtr.(fldLabels{2}); res.jtr.(fldLabels{3})];
    [hJtr,~]=adtest(allRes.jtr);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(jtr.(fldLabels{ii}));
    end
    if noGroups == 1
        hJtr = holdH.(fldLabels{1});
    elseif noGroups == 2
        hJtr = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%jitter, comparison tests
if noGroups == 3
    if hJtr == 0 %if residuals normal, run anova
        [output.p.jtrAnova,~,jtrstats]=anova1([jtr.(fldLabels{1})' jtr.(fldLabels{2})' jtr.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.jtrAnova <= 0.05
            output.jtrPostHocTukey=multcompare(jtrstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.jtrKW,~,~]=kruskalwallis([jtr.(fldLabels{1})' jtr.(fldLabels{2})' jtr.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.jtrKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.jtrPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(jtr.(fldLabels{1}),jtr.(fldLabels{2}));
            [output.p.jtrPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(jtr.(fldLabels{1}),jtr.(fldLabels{3}));
            [output.p.jtrPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(jtr.(fldLabels{2}),jtr.(fldLabels{3}));
            output.pFDR.jtr=drsFDRpval([output.p.jtrPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.jtrPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.jtrPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hJtr == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,output.p.jtrTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(jtr.(fldLabels{1}),jtr.(fldLabels{2}));
        end
    elseif hJtr ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [output.p.jtrMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(jtr.(fldLabels{1}),jtr.(fldLabels{2}));
        end
    end
end

%rise time, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.rt.(fldLabels{ii}) = mean(rt.(fldLabels{ii})) - rt.(fldLabels{ii});
    end
    allRes.rt = [res.rt.(fldLabels{1}); res.rt.(fldLabels{2}); res.rt.(fldLabels{3})];
    [hRT,~]=adtest(allRes.rt);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(rt.(fldLabels{ii}));
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
        [output.p.rtAnova,~,rtstats]=anova1([rt.(fldLabels{1})' rt.(fldLabels{2})' rt.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.rtAnova <= 0.05
            output.rtPostHocTukey=multcompare(rtstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.rtKW,~,~]=kruskalwallis([rt.(fldLabels{1})' rt.(fldLabels{2})' rt.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.rtKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(rt.(fldLabels{1}),rt.(fldLabels{2}));
            [output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(rt.(fldLabels{1}),rt.(fldLabels{3}));
            [output.p.rtPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(rt.(fldLabels{2}),rt.(fldLabels{3}));
            output.pFDR.rt=drsFDRpval([output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.rtPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hRT == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,output.p.rtTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(rt.(fldLabels{1}),rt.(fldLabels{2}));
        end
    elseif hRT ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [output.p.rtMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(rt.(fldLabels{1}),rt.(fldLabels{2}));
        end
    end
end

%decay tau, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.dtau.(fldLabels{ii}) = mean(dtau.(fldLabels{ii})) - dtau.(fldLabels{ii});
    end
    allRes.dtau = [res.dtau.(fldLabels{1}); res.dtau.(fldLabels{2}); res.dtau.(fldLabels{3})];
    [hDecay,~]=adtest(allRes.dtau);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(dtau.(fldLabels{ii}));
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
        [output.p.dtauAnova,~,dtaustats]=anova1([dtau.(fldLabels{1})' dtau.(fldLabels{2})' dtau.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.dtauAnova <= 0.05
            output.dtauPostHocTukey=multcompare(dtaustats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.dtauKW,~,~]=kruskalwallis([dtau.(fldLabels{1})' dtau.(fldLabels{2})' dtau.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.dtauKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.dtauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(dtau.(fldLabels{1}),dtau.(fldLabels{2}));
            [output.p.dtauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(dtau.(fldLabels{1}),dtau.(fldLabels{3}));
            [output.p.dtauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(dtau.(fldLabels{2}),dtau.(fldLabels{3}));
            output.pFDR.dtau=drsFDRpval([output.p.dtauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.dtauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.dtauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hDecay == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,output.p.dtauTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(dtau.(fldLabels{1}),dtau.(fldLabels{2}));
        end
    elseif hDecay ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [output.p.dtauMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(dtau.(fldLabels{1}),dtau.(fldLabels{2}));
        end
    end
end

%stim intensity, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.stim.(fldLabels{ii}) = mean(stim.(fldLabels{ii})) - stim.(fldLabels{ii});
    end
    allRes.stim = [res.stim.(fldLabels{1}); res.stim.(fldLabels{2}); res.stim.(fldLabels{3})];
    [hStim,~]=adtest(allRes.stim);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(stim.(fldLabels{ii}));
    end
    if noGroups == 1
        hStim = holdH.(fldLabels{1});
    elseif noGroups == 2
        hStim = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%stim intensity, comparison tests
if noGroups == 3
    if hStim == 0 %if residuals normal, run anova
        [output.p.stimAnova,~,stimstats]=anova1([stim.(fldLabels{1})' stim.(fldLabels{2})' stim.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.stimAnova <= 0.05
            output.stimPostHocTukey=multcompare(stimstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.stimKW,~,~]=kruskalwallis([stim.(fldLabels{1})' stim.(fldLabels{2})' stim.(fldLabels{3})'],...
            [ones(1,length(minstimFiles.(fldLabels{1}).ctrl)) 2.*ones(1,length(minstimFiles.(fldLabels{2}).ctrl)) 3.*ones(1,length(minstimFiles.(fldLabels{3}).ctrl))]);
        if output.p.stimKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.stimPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(stim.(fldLabels{1}),stim.(fldLabels{2}));
            [output.p.stimPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(stim.(fldLabels{1}),stim.(fldLabels{3}));
            [output.p.stimPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(stim.(fldLabels{2}),stim.(fldLabels{3}));
            output.pFDR.stim=drsFDRpval([output.p.stimPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.stimPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.stimPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hStim == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,output.p.stimTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(stim.(fldLabels{1}),stim.(fldLabels{2}));
        end
    elseif hStim ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [output.p.stimMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(stim.(fldLabels{1}),stim.(fldLabels{2}));
        end
    end
end

%success rate, effect of glutamate receptor blockade
for ii = 1:noGroups
    hHold.(fldLabels{ii}).ctrl = adtest(success.(fldLabels{ii}).ctrl);
    hHold.(fldLabels{ii}).xv = adtest(success.(fldLabels{ii}).xv);
    hSuccess.(fldLabels{ii}) = sum([hHold.(fldLabels{ii}).ctrl hHold.(fldLabels{ii}).xv]);
    if hSuccess.(fldLabels{ii}) == 0 %if data normal, run ttest
        [~,output.p.sRateTtest.(fldLabels{ii})] = ttest2(success.(fldLabels{ii}).ctrl,success(fldLabels{ii}).xv);
    elseif hSuccess.(fldLabels{ii}) ~= 0 %if data not normal, run MW U test
        [output.p.sRateMW.(fldLabels{ii}),~] = ranksum(success.(fldLabels{ii}).ctrl,success.(fldLabels{ii}).xv);
    end
end
clear holdH

close all %closes anova and kw figures

%% DISPLAY DATA
%color
lightgray=[.75 .75 .75]; %light gray for individual data points

%trace figures
t=(1000/inputs.SamplerateHzEditField.Value).*(1:1:length(meanTrace.(fldLabels{ii})));
for ii = 1:noGroups
    mTraceFig.(fldLabels{ii}) = figure(ii);
    mTraceFig.(fldLabels{ii}).Position = [650+(ii-1)*250 670 225 125];
    hold on
    plot(t,meanTrace.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',2)
    traceAx.(fldLabels{ii}) = gca;
    setAx(traceAx.(fldLabels{ii}));
    ylabel('current (pA)')
    xlabel('t (ms)')
    title([fldLabels{ii} ' example mean response'])
    xlim([12.5 42.5])
    if ii == 1
        yLimHolder = traceAx.(fldLabels{ii}).YLim;
    else
        yLimHolder = [yLimHolder traceAx.(fldLabels{ii}).YLim];
    end
end
yLimMin = min(yLimHolder);
yLimMax = max(yLimHolder);
for ii = 1:noGroups
    traceAx.(fldLabels{ii}).YLim = [yLimMin yLimMax];
end

%amplitude figure
ampFig=figure(noGroups+1);
ampFig.Position = [140 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(amp.(fldLabels{ii}))),amp.(fldLabels{ii}),40,lightgray,'filled')
    if hAmp == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(amp.(fldLabels{ii})),sem(amp.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(amp.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hAmp ~= 0 %if data non-normal, plot median w IQR
        [lowerError.amp.(fldLabels{ii}),upperError.amp.(fldLabels{ii})] = iqrError(amp.(fldLabels{ii}),1);
        errorbar(ii,median(amp.(fldLabels{ii})),lowerError.amp.(fldLabels{ii}),upperError.amp.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(amp.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allAmp = amp.(fldLabels{ii});
    else
        allAmp = [allAmp; amp.(fldLabels{ii})];
    end
end
ylabel('uEPSC amplitude (pA)')
ampAx=gca;
setAx(ampAx); 
ampAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 1.2*max(allAmp)])


%latency figure
latFig=figure(noGroups+2);
latFig.Position = [290 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(lat.(fldLabels{ii}))),lat.(fldLabels{ii}),40,lightgray,'filled')
    if hLat == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(lat.(fldLabels{ii})),sem(lat.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(lat.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hLat ~= 0 %if data non-normal, plot median w IQR
        [lowerError.lat.(fldLabels{ii}),upperError.lat.(fldLabels{ii})] = iqrError(lat.(fldLabels{ii}),1);
        errorbar(ii,median(lat.(fldLabels{ii})),lowerError.lat.(fldLabels{ii}),upperError.lat.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(lat.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allLat = lat.(fldLabels{ii});
    else
        allLat = [allLat; lat.(fldLabels{ii})];
    end
end
ylabel('uEPSC latency (ms)')
latAx=gca;
setAx(latAx); 
latAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 1.2*max(allLat)])

%jitter figure
jtrFig=figure(noGroups+3);
jtrFig.Position = [440 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(jtr.(fldLabels{ii}))),jtr.(fldLabels{ii}),40,lightgray,'filled')
    if hJtr == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(jtr.(fldLabels{ii})),sem(jtr.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(jtr.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hJtr ~= 0 %if data non-normal, plot median w IQR
        [lowerError.jtr.(fldLabels{ii}),upperError.jtr.(fldLabels{ii})] = iqrError(jtr.(fldLabels{ii}),1);
        errorbar(ii,median(jtr.(fldLabels{ii})),lowerError.jtr.(fldLabels{ii}),upperError.jtr.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(jtr.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allJtr = jtr.(fldLabels{ii});
    else
        allJtr = [allJtr; jtr.(fldLabels{ii})];
    end
end
ylabel('uEPSC jitter (ms)')
jtrAx=gca;
setAx(jtrAx); 
jtrAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 2*ceil(max(allJtr))])

%rise time figure
rtFig=figure(noGroups+4);
rtFig.Position = [590 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(rt.(fldLabels{ii}))),rt.(fldLabels{ii}),40,lightgray,'filled')
    if hRT == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(rt.(fldLabels{ii})),sem(rt.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(rt.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hRT ~= 0 %if data non-normal, plot median w IQR
        [lowerError.rt.(fldLabels{ii}),upperError.rt.(fldLabels{ii})] = iqrError(rt.(fldLabels{ii}),1);
        errorbar(ii,median(rt.(fldLabels{ii})),lowerError.rt.(fldLabels{ii}),upperError.rt.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(rt.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allRT = rt.(fldLabels{ii});
    else
        allRT = [allRT; rt.(fldLabels{ii})];
    end
end
ylabel('uEPSC 20%-80% risetime (ms)')
rtAx=gca;
setAx(rtAx); 
rtAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 2*ceil(max(allRT))])

%decay tau figure
dtauFig=figure(noGroups+5);
dtauFig.Position = [740 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(dtau.(fldLabels{ii}))),dtau.(fldLabels{ii}),40,lightgray,'filled')
    if hDecay == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(dtau.(fldLabels{ii})),sem(dtau.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(dtau.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hDecay ~= 0 %if data non-normal, plot median w IQR
        [lowerError.dtau.(fldLabels{ii}),upperError.dtau.(fldLabels{ii})] = iqrError(dtau.(fldLabels{ii}),1);
        errorbar(ii,median(dtau.(fldLabels{ii})),lowerError.dtau.(fldLabels{ii}),upperError.dtau.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(dtau.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allDecay = dtau.(fldLabels{ii});
    else
        allDecay = [allDecay; dtau.(fldLabels{ii})];
    end
end
ylabel('uEPSC decay \tau (ms)')
dtauAx=gca;
setAx(dtauAx); 
dtauAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 1.2*max(allDecay)])

%stimulus intensity figure
stimFig=figure(noGroups+6);
stimFig.Position = [890 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(stim.(fldLabels{ii}))),stim.(fldLabels{ii}),40,lightgray,'filled')
    if hStim == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(stim.(fldLabels{ii})),sem(stim.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(stim.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hStim ~= 0 %if data non-normal, plot median w IQR
        [lowerError.stim.(fldLabels{ii}),upperError.stim.(fldLabels{ii})] = iqrError(stim.(fldLabels{ii}),1);
        errorbar(ii,median(stim.(fldLabels{ii})),lowerError.stim.(fldLabels{ii}),upperError.stim.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(stim.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allStim = stim.(fldLabels{ii});
    else
        allStim = [allStim; stim.(fldLabels{ii})];
    end
end
ylabel('uEPSC threshold stimulation (\muA*ms)')
stimAx=gca;
setAx(stimAx); 
stimAx.XTick=[];
stimAx.YScale = 'log';
xlim([0.5 noGroups+.5])
ylim([10^0 10^(ceil(log10(max(allStim))))])

%confirmation of glutamatergic origin of events
for ii = 1:noGroups
    %lighter color
    lightColor = (1+(.85-max(inputs.plotColorGrps(ii,:)))./max(inputs.plotColorGrps(ii,:))).*inputs.plotColorGrps(ii,:); %linearly scales RGB values to a point where max value is .85 for one of R,G,B
    %plot
    glutFig.(fldLabels{ii}) = figure(noGroups+6+ii);
    glutFig.(fldLabels{ii}).Position = [1030+(ii-1)*125 385 120 170];
    hold on
    scatter(ones(1,length(success.(fldLabels{ii}).ctrl)),success.(fldLabels{ii}).ctrl,40,lightgray,'filled')
    scatter(2.*ones(1,length(success.(fldLabels{ii}).xv)),success.(fldLabels{ii}).xv,40,lightgray,'filled')
    if hSuccess.(fldLabels{ii}) == 0 %if data normal, plot mean with sem
        errorbar(1,mean(success.(fldLabels{ii}).ctrl),sem(success.(fldLabels{ii}).ctrl,1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(1,mean(success.(fldLabels{ii}).ctrl),125,inputs.plotColorGrps(ii,:),'filled')
        errorbar(2,mean(success.(fldLabels{ii}).xv),sem(success.(fldLabels{ii}).xv,1),'color',lightColor,'linewidth',3,'CapSize',0)
        scatter(2,mean(success.(fldLabels{ii}).xv),125,lightColor,'filled')
    elseif hSuccess.(fldLabels{ii}) ~= 0 %if data not normal, plot median with iqr
        [lowerError.success.(fldLabels{ii}).ctrl,upperError.success.(fldLabels{ii}).ctrl] = iqrError(success.(fldLabels{ii}).ctrl,1);
        [lowerError.success.(fldLabels{ii}).xv,upperError.success.(fldLabels{ii}).xv] = iqrError(success.(fldLabels{ii}).xv,1);
        errorbar(1,median(success.(fldLabels{ii}).ctrl),lowerError.success.(fldLabels{ii}).ctrl,upperError.success.(fldLabels{ii}).ctrl,'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(1,median(success.(fldLabels{ii}).ctrl),125,inputs.plotColorGrps(ii,:),'filled')
        errorbar(2,median(success.(fldLabels{ii}).xv),lowerError.success.(fldLabels{ii}).xv,upperError.success.(fldLabels{ii}).xv,'color',lightColor,'linewidth',3,'CapSize',0)
        scatter(2,median(success.(fldLabels{ii}).xv),125,lightColor,'filled')
    end
    xlim([.5 2.5])
    ylim([0 100])
    glutAx.(fldLabels{ii}) = gca;
    setAx(glutAx.(fldLabels{ii}));
    glutAx.(fldLabels{ii}).YTick = [0 25 50 75 100];
    glutAx.(fldLabels{ii}).XTickLabel = {'ctrl', '+dqnx/apv'};
    ylabel('Success Rate')
    title(fldLabels{ii})
end

%proportion responding for each cell type
prFig = figure(noGroups+10);
prFig.Position = [1030 55 133 225];
hold on
for ii = 1:noGroups
    bar(ii,proportionResp.(fldLabels{ii}),'facecolor',inputs.plotColorGrps(ii,:),'edgecolor',inputs.plotColorGrps(ii,:),'linewidth',1)
    prAx = gca;
    setAx(prAx);
    prAx.XTick = [];
    prAx.YTick = [0 .25 .5 .75 1];
    ylim([0 1])
    ylabel('Proportion of cells responding to stimulus')
end
end