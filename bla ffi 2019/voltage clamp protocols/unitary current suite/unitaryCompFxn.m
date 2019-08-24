function [output] = unitaryCompFxn(inputs)
%for group comparison of minimal stimulation data
%modeled off ms_dispstats.m code
%created 02-06-2019, edited 04-06-19

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

%% LOAD DATA
% Directory
cd(inputs.upscDir);

for ii = 1:noGroups
    %get files
    cd([inputs.upscDir,'/',grpLabels{ii},'/',fldrName]);
    contents.(fldLabels{ii}) = dir('*.mat');
    filenames.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
    upscFiles.(fldLabels{ii}) = fullfile(cd,filenames.(fldLabels{ii}));
    
    %inits
    amp.(fldLabels{ii})=-1.*ones(length(upscFiles.(fldLabels{ii})),1);
    lat.(fldLabels{ii})=-1.*ones(length(upscFiles.(fldLabels{ii})),1);
    jtr.(fldLabels{ii})=-1.*ones(length(upscFiles.(fldLabels{ii})),1);
    rt.(fldLabels{ii})=-1.*ones(length(upscFiles.(fldLabels{ii})),1);
    dtau.(fldLabels{ii})=-1.*ones(length(upscFiles.(fldLabels{ii})),1);
    success.(fldLabels{ii})=-1.*ones(length(upscFiles.(fldLabels{ii})),1);
    
    count.(fldLabels{ii}) = 0;
    for jj = 1:length(upscFiles.(fldLabels{ii}))
        %load ctrl data
        load(upscFiles.(fldLabels{ii}){jj}, 'app')
        success.(fldLabels{ii})(jj) = app.upscData.sRate;
        if success.(fldLabels{ii})(jj) == 0
            count.(fldLabels{ii}) = count.(fldLabels{ii}) +1;
            failures.(fldLabels{ii})(count.(fldLabels{ii})) = jj;
        elseif success.(fldLabels{ii})(jj) ~= 0
            amp.(fldLabels{ii})(jj) = app.upscData.mAmp;
            lat.(fldLabels{ii})(jj) = app.upscData.mLat;
            jtr.(fldLabels{ii})(jj) = app.upscData.jitter;
            rt.(fldLabels{ii})(jj) = app.upscData.mRT;
            dtau.(fldLabels{ii})(jj) = app.upscData.decay;
        end
        
        if ii == 1
            exCell.(fldLabels{ii}) = inputs.ExCellEditField.Value;
        elseif ii == 2
            exCell.(fldLabels{ii}) = inputs.ExCellEditField_2.Value;
        elseif ii == 3
            exCell.(fldLabels{ii}) = inputs.ExCellEditField_3.Value;
        end
        if jj == exCell.(fldLabels{ii})
            traces.(fldLabels{ii}) = app.upscData.wholeTrace;
            exCellAmps.(fldLabels{ii}) = app.upscData.amplitude;
            meanTrace.(fldLabels{ii}) = app.upscData.mTrace;
            spkTrace.(fldLabels{ii}) = app.upscData.meanSpkTrace;
        end
        clear app
    end
    %get rid of "data" for cells without a response
    if sum(success.(fldLabels{ii})==0) > 0
        amp.(fldLabels{ii})(failures.(fldLabels{ii})) = [];
        lat.(fldLabels{ii})(failures.(fldLabels{ii})) = [];
        jtr.(fldLabels{ii})(failures.(fldLabels{ii})) = [];
        rt.(fldLabels{ii})(failures.(fldLabels{ii})) = [];
        dtau.(fldLabels{ii})(failures.(fldLabels{ii})) = [];
        success.(fldLabels{ii})(failures.(fldLabels{ii})) = [];
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
            [ones(1,length(amp.(fldLabels{1}))) 2.*ones(1,length(amp.(fldLabels{2}))) 3.*ones(1,length(amp.(fldLabels{3})))]);
        if output.p.ampAnova <= 0.05
            output.ampPostHocTukey=multcompare(ampstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.ampKW,~,~]=kruskalwallis([amp.(fldLabels{1})' amp.(fldLabels{2})' amp.(fldLabels{3})'],...
            [ones(1,length(amp.(fldLabels{1}))) 2.*ones(1,length(amp.(fldLabels{2}))) 3.*ones(1,length(amp.(fldLabels{3})))]);
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
            [ones(1,length(lat.(fldLabels{1}))) 2.*ones(1,length(lat.(fldLabels{2}))) 3.*ones(1,length(lat.(fldLabels{3})))]);
        if output.p.latAnova <= 0.05
            output.latPostHocTukey=multcompare(latstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.latKW,~,~]=kruskalwallis([lat.(fldLabels{1})' lat.(fldLabels{2})' lat.(fldLabels{3})'],...
            [ones(1,length(lat.(fldLabels{1}))) 2.*ones(1,length(lat.(fldLabels{2}))) 3.*ones(1,length(lat.(fldLabels{3})))]);
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
            [ones(1,length(jtr.(fldLabels{1}))) 2.*ones(1,length(jtr.(fldLabels{2}))) 3.*ones(1,length(jtr.(fldLabels{3})))]);
        if output.p.jtrAnova <= 0.05
            output.jtrPostHocTukey=multcompare(jtrstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.jtrKW,~,~]=kruskalwallis([jtr.(fldLabels{1})' jtr.(fldLabels{2})' jtr.(fldLabels{3})'],...
            [ones(1,length(jtr.(fldLabels{1}))) 2.*ones(1,length(jtr.(fldLabels{2}))) 3.*ones(1,length(jtr.(fldLabels{3})))]);
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
            [ones(1,length(rt.(fldLabels{1}))) 2.*ones(1,length(rt.(fldLabels{2}))) 3.*ones(1,length(rt.(fldLabels{3})))]);
        if output.p.rtAnova <= 0.05
            output.rtPostHocTukey=multcompare(rtstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.rtKW,~,~]=kruskalwallis([rt.(fldLabels{1})' rt.(fldLabels{2})' rt.(fldLabels{3})'],...
            [ones(1,length(rt.(fldLabels{1}))) 2.*ones(1,length(rt.(fldLabels{2}))) 3.*ones(1,length(rt.(fldLabels{3})))]);
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
            [ones(1,length(dtau.(fldLabels{1}))) 2.*ones(1,length(dtau.(fldLabels{2}))) 3.*ones(1,length(dtau.(fldLabels{3})))]);
        if output.p.dtauAnova <= 0.05
            output.dtauPostHocTukey=multcompare(dtaustats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.dtauKW,~,~]=kruskalwallis([dtau.(fldLabels{1})' dtau.(fldLabels{2})' dtau.(fldLabels{3})'],...
            [ones(1,length(dtau.(fldLabels{1}))) 2.*ones(1,length(dtau.(fldLabels{2}))) 3.*ones(1,length(dtau.(fldLabels{3})))]);
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

%success rate, normality tests
if noGroups == 3
    for ii = 1:noGroups
        res.success.(fldLabels{ii}) = mean(success.(fldLabels{ii})) - success.(fldLabels{ii});
    end
    allRes.success = [res.success.(fldLabels{1}); res.success.(fldLabels{2}); res.success.(fldLabels{3})];
    [hSuccess,~]=adtest(allRes.success);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(success.(fldLabels{ii}));
    end
    if noGroups == 1
        hSuccess = holdH.(fldLabels{1});
    elseif noGroups == 2
        hSuccess = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%success rate
if noGroups == 3
    if hSuccess == 0 %if residuals normal, run anova
        [output.p.successAnova,~,successstats]=anova1([success.(fldLabels{1})' success.(fldLabels{2})' success.(fldLabels{3})'],...
            [ones(1,length(success.(fldLabels{1}))) 2.*ones(1,length(success.(fldLabels{2}))) 3.*ones(1,length(success.(fldLabels{3})))]);
        if output.p.successAnova <= 0.05
            output.successPostHocTukey=multcompare(successstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.successKW,~,~]=kruskalwallis([success.(fldLabels{1})' success.(fldLabels{2})' success.(fldLabels{3})'],...
            [ones(1,length(success.(fldLabels{1}))) 2.*ones(1,length(success.(fldLabels{2}))) 3.*ones(1,length(success.(fldLabels{3})))]);
        if output.p.successKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.successPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(success.(fldLabels{1}),success.(fldLabels{2}));
            [output.p.successPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(success.(fldLabels{1}),success.(fldLabels{3}));
            [output.p.successPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(success.(fldLabels{2}),success.(fldLabels{3}));
            output.pFDR.success=drsFDRpval([output.p.successPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.successPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.successPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hSuccess == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,output.p.successTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(success.(fldLabels{1}),success.(fldLabels{2}));
        end
    elseif hSuccess ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [output.p.successMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(success.(fldLabels{1}),success.(fldLabels{2}));
        end
    end
end

close all %closes anova and kw figures

%% DISPLAY DATA
%color
lightgray=[.75 .75 .75]; %light gray for individual data points

%trace figures
t=(1000/inputs.SamplerateHzEditField.Value).*(1:1:size(traces.(fldLabels{ii}),1));
for ii = 1:noGroups
    mTraceFig.(fldLabels{ii}) = figure(ii);
    mTraceFig.(fldLabels{ii}).Position = [150+(ii-1)*450 670 425 425];
    subplot(2,1,1)
    hold on
    plot(t,traces.(fldLabels{ii})(:,isnan(exCellAmps.(fldLabels{ii}))),'color','k','linewidth',2)
    plot(t,traces.(fldLabels{ii})(:,~isnan(exCellAmps.(fldLabels{ii}))),'color',inputs.plotColorGrps(ii,:),'linewidth',2)
    traceAx.(fldLabels{ii}) = gca;
    setAx(traceAx.(fldLabels{ii}));
    ylabel('current (pA)')
    xlabel('t (ms)')
    title([fldLabels{ii} ' example unitary response'])
    xlim([2060 2090])
    if ii == 1
        yLimHolder = traceAx.(fldLabels{ii}).YLim;
    else
        yLimHolder = [yLimHolder traceAx.(fldLabels{ii}).YLim];
    end
    
    subplot(2,1,2)
    plot(t,spkTrace.(fldLabels{ii}),'color',[.2 .2 .2],'linewidth',2)
    ylabel('voltage (mV)')
    xlabel('t (ms)')
    xlim([2056 2106])
    ylim([-75 45])
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
ylim([1.2*min(allAmp) 0])


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

%success rate figure
successFig=figure(noGroups+6);
successFig.Position = [890 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(success.(fldLabels{ii}))),success.(fldLabels{ii}),40,lightgray,'filled')
    if hSuccess == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(success.(fldLabels{ii})),sem(success.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(success.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hSuccess ~= 0 %if data non-normal, plot median w IQR
        [lowerError.success.(fldLabels{ii}),upperError.success.(fldLabels{ii})] = iqrError(success.(fldLabels{ii}),1);
        errorbar(ii,median(success.(fldLabels{ii})),lowerError.success.(fldLabels{ii}),upperError.success.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(success.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allSuccess = success.(fldLabels{ii});
    else
        allSuccess = [allSuccess; success.(fldLabels{ii})];
    end
end
ylabel('Success Rate (%)')
successAx=gca;
setAx(successAx); 
successAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([0 100])

output
output.p %#ok<*NOPRT>
end