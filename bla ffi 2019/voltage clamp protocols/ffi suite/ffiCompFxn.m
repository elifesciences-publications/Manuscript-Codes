%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FFI comp For Suite %%%%%%%%%%%%%
%%%%%%%%%%% Created: 05-06-2018 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 05-14-2019 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = ffiCompFxn(inputs)
%modeled after ffi_dispstats_dreadds.m

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
cd(inputs.ffiDir);
for ii = 1:noGroups
    %get files,
    cd([inputs.ffiDir,'/',grpLabels{ii}]);
    contents.(fldLabels{ii}) = dir('*.mat');
    filenames.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
    ffiFiles.(fldLabels{ii}) = fullfile(cd,filenames.(fldLabels{ii}));
    
    for jj = 1:length(ffiFiles.(fldLabels{ii}))
        %load data
        load(ffiFiles.(fldLabels{ii}){jj}, 'app')
        if isnan(app.ffiData.decay.epsc)
            if isnan(app.ffiData.rise.epsc) %if these two things occur in same file, it is due to lack of events
                amp.(fldLabels{ii}).epsc(jj) = app.ffiData.mAmp.epsc;
            else %if events occurred, the EPSC was recorded with a positive amplitude due to nature of code, switch to a negative value
                amp.(fldLabels{ii}).epsc(jj) = -1.*app.ffiData.mAmp.epsc;
            end
        else
            amp.(fldLabels{ii}).epsc(jj) = -1.*app.ffiData.mAmp.epsc;
        end
        amp.(fldLabels{ii}).ipsc(jj) = app.ffiData.mAmp.ipsc; %always returned as positive value
        ieRat.(fldLabels{ii})(jj) = amp.(fldLabels{ii}).ipsc(jj)./abs(amp.(fldLabels{ii}).epsc(jj));
        if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
            if strcmp(fldLabels{ii},'ctrl')
                lat.(fldLabels{ii}).epsc(jj) = app.ffiData.mLat.epsc;
                lat.(fldLabels{ii}).ipsc(jj) = app.ffiData.mLat.ipsc;
                jitter.(fldLabels{ii}).epsc(jj) = app.ffiData.jtr.epsc;
                jitter.(fldLabels{ii}).ipsc(jj) = app.ffiData.jtr.ipsc;
            end
        else
            lat.(fldLabels{ii}).epsc(jj) = app.ffiData.mLat.epsc;
            lat.(fldLabels{ii}).ipsc(jj) = app.ffiData.mLat.ipsc;
            jitter.(fldLabels{ii}).epsc(jj) = app.ffiData.jtr.epsc;
            jitter.(fldLabels{ii}).ipsc(jj) = app.ffiData.jtr.ipsc;
        end
        if inputs.KineticsAnalysisCheckBox.Value == 1
            rt.(fldLabels{ii}).epsc(jj) = mean(app.ffiData.rise.epsc(~isnan(app.ffiData.rise.epsc)));
            rt.(fldLabels{ii}).ipsc(jj) = mean(app.ffiData.rise.ipsc(~isnan(app.ffiData.rise.ipsc)));
            tau.(fldLabels{ii}).epsc(jj) = app.ffiData.decay.epsc;
            tau.(fldLabels{ii}).ipsc(jj) = app.ffiData.decay.ipsc;
        end
        clear app
    end
end

%% SELECT EXAMPLE TRACE
exNeuron = inputs.ExCellEditField.Value;
for ii = 1:noGroups
    for jj = 1:length(ffiFiles.(fldLabels{ii}))
        thisFile = filenames.(fldLabels{ii})(jj);
        shortNames.(fldLabels{ii}){jj} = thisFile{:}(1:end-4-length(fldLabels{ii}));
        clear thisFile
    end
end
%check if example neuron is present in each condition
for ii = 2:noGroups
    for jj = 1:length(shortNames.(fldLabels{ii}))
        sameName.(fldLabels{ii})(jj) = strcmp(shortNames.(fldLabels{1}){exNeuron},shortNames.(fldLabels{ii}){jj});
    end
    if sum(sameName.(fldLabels{ii})) == 0
        error(['example neuron not part of ',fldLabels{ii},' condition, choose another'])
    end
end


%pick example neuron and grab its mean traces
for ii = 1:noGroups
    %get files
    cd([inputs.ffiDir,'/',grpLabels{ii}]);
    
    %load mean traces
    if ii == 1
        load(ffiFiles.(fldLabels{ii}){exNeuron}, 'app')
    else
        thisCell = find(sameName.(fldLabels{ii}) == 1);
        load(ffiFiles.(fldLabels{ii}){thisCell}, 'app')
    end
    exTrace.(fldLabels{ii}).epsc = app.ffiData.mTrace.epsc';
    exTrace.(fldLabels{ii}).ipsc = app.ffiData.mTrace.ipsc';
    clear app
end


%% TEST DIFFERENCES
%adtest (Anderson-Darling Test) to test for normality
%amplitude, normality
if noGroups == 3
    for ii = 1:noGroups
        res.amp.(fldLabels{ii}).epsc = mean(amp.(fldLabels{ii}).epsc) - amp.(fldLabels{ii}).epsc;
        res.amp.(fldLabels{ii}).ipsc = mean(amp.(fldLabels{ii}).ipsc) - amp.(fldLabels{ii}).ipsc;
    end
    allRes.amp.epsc = [res.amp.(fldLabels{1}).epsc res.amp.(fldLabels{2}).epsc res.amp.(fldLabels{3}).epsc];
    allRes.amp.ipsc = [res.amp.(fldLabels{1}).ipsc res.amp.(fldLabels{2}).ipsc res.amp.(fldLabels{3}).ipsc];
    [hAmp.epsc,~]=adtest(allRes.amp.epsc);
    [hAmp.ipsc,~]=adtest(allRes.amp.ipsc);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}).epsc,~] = adtest(amp.(fldLabels{ii}).epsc);
        [holdH.(fldLabels{ii}).ipsc,~] = adtest(amp.(fldLabels{ii}).ipsc);
    end
    if noGroups == 1
        hAmp.epsc = holdH.(fldLabels{1}).epsc;
        hAmp.ipsc = holdH.(fldLabels{1}).ipsc;
    elseif noGroups == 2
        hAmp.epsc = sum([holdH.(fldLabels{1}).epsc holdH.(fldLabels{2}).epsc]);
        hAmp.ipsc = sum([holdH.(fldLabels{1}).ipsc holdH.(fldLabels{2}).ipsc]);
    end
    clear holdH
end

%amplitude, comparison tests
if noGroups == 3
    if hAmp.epsc == 0 %if residuals normal, run anova
        [output.p.ampAnova.epsc,~,ampstats.epsc]=anova1([amp.(fldLabels{1}).epsc amp.(fldLabels{2}).epsc amp.(fldLabels{3}).epsc],...
            [ones(1,length(amp.(fldLabels{1}).epsc)) 2.*ones(1,length(amp.(fldLabels{2}).epsc)) 3.*ones(1,length(amp.(fldLabels{3}).epsc))]);
        if output.p.ampAnova.epsc <= 0.05
            output.ampPostHocTukey.epsc=multcompare(ampstats.epsc); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.ampKW.epsc,~,~]=kruskalwallis([amp.(fldLabels{1}).epsc amp.(fldLabels{2}).epsc amp.(fldLabels{3}).epsc],...
            [ones(1,length(amp.(fldLabels{1}).epsc)) 2.*ones(1,length(amp.(fldLabels{2}).epsc)) 3.*ones(1,length(amp.(fldLabels{3}).epsc))]);
        if output.p.ampKW.epsc <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=ranksum(amp.(fldLabels{1}).epsc,amp.(fldLabels{2}).epsc);
            [output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).epsc,~]=ranksum(amp.(fldLabels{1}).epsc,amp.(fldLabels{3}).epsc);
            [output.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).epsc,~]=ranksum(amp.(fldLabels{2}).epsc,amp.(fldLabels{3}).epsc);
            output.pFDR.amp.epsc=drsFDRpval([output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).epsc output.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).epsc]);
        end
    end
    if hAmp.ipsc == 0 %if residuals normal, run anova
        [output.p.ampAnova.ipsc,~,ampstats.ipsc]=anova1([amp.(fldLabels{1}).ipsc amp.(fldLabels{2}).ipsc amp.(fldLabels{3}).ipsc],...
            [ones(1,length(amp.(fldLabels{1}).ipsc)) 2.*ones(1,length(amp.(fldLabels{2}).ipsc)) 3.*ones(1,length(amp.(fldLabels{3}).ipsc))]);
        if output.p.ampAnova.ipsc <= 0.05
            output.ampPostHocTukey.ipsc=multcompare(ampstats.ipsc); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.ampKW.ipsc,~,~]=kruskalwallis([amp.(fldLabels{1}).ipsc amp.(fldLabels{2}).ipsc amp.(fldLabels{3}).ipsc],...
            [ones(1,length(amp.(fldLabels{1}).ipsc)) 2.*ones(1,length(amp.(fldLabels{2}).ipsc)) 3.*ones(1,length(amp.(fldLabels{3}).ipsc))]);
        if output.p.ampKW.ipsc <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=ranksum(amp.(fldLabels{1}).ipsc,amp.(fldLabels{2}).ipsc);
            [output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).ipsc,~]=ranksum(amp.(fldLabels{1}).ipsc,amp.(fldLabels{3}).ipsc);
            [output.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).ipsc,~]=ranksum(amp.(fldLabels{2}).ipsc,amp.(fldLabels{3}).ipsc);
            output.pFDR.amp.ipsc=drsFDRpval([output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc output.p.ampPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).ipsc output.p.ampPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).ipsc]);
        end
    end
elseif noGroups < 3
    if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 0
        if hAmp.epsc == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,output.p.ampTtest.([fldLabels{1} 'vs' fldLabels{2}]).epsc]=ttest2(amp.(fldLabels{1}).epsc,amp.(fldLabels{2}).epsc);
            end
        elseif hAmp.epsc ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [output.p.ampMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=ranksum(amp.(fldLabels{1}).epsc,amp.(fldLabels{2}).epsc);
            end
        end
        if hAmp.ipsc == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,output.p.ampTtest.([fldLabels{1} 'vs' fldLabels{2}]).ipsc]=ttest2(amp.(fldLabels{1}).ipsc,amp.(fldLabels{2}).ipsc);
            end
        elseif hAmp.ipsc ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [output.p.ampMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=ranksum(amp.(fldLabels{1}).ipsc,amp.(fldLabels{2}).ipsc);
            end
        end
    elseif inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
        if hAmp.epsc == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,output.p.ampPairedttest.([fldLabels{1} 'vs' fldLabels{2}]).epsc]=ttest(amp.(fldLabels{1}).epsc,amp.(fldLabels{2}).epsc);
            end
        elseif hAmp.epsc ~= 0 %if data not normal, run wilcoxon signed rank test
            if noGroups == 2
                [output.p.ampWSR.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=signrank(amp.(fldLabels{1}).epsc,amp.(fldLabels{2}).epsc);
            end
        end
        if hAmp.ipsc == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,output.p.ampPairedttest.([fldLabels{1} 'vs' fldLabels{2}]).ipsc]=ttest(amp.(fldLabels{1}).ipsc,amp.(fldLabels{2}).ipsc);
            end
        elseif hAmp.ipsc ~= 0 %if data not normal, run wilcoxon signed rank test
            if noGroups == 2
                [output.p.ampWSR.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=signrank(amp.(fldLabels{1}).ipsc,amp.(fldLabels{2}).ipsc);
            end
        end
    end
end

%ie ratio, normality
if noGroups == 3
    for ii = 1:noGroups
        res.ieRat.(fldLabels{ii}) = mean(ieRat.(fldLabels{ii})) - ieRat.(fldLabels{ii});
    end
    allRes.ieRat = [res.ieRat.(fldLabels{1}) res.ieRat.(fldLabels{2}) res.ieRat.(fldLabels{3})];
    [hieRat,~]=adtest(allRes.ieRat);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}).ieRat,~] = adtest(ieRat.(fldLabels{ii}));
    end
    if noGroups == 1
        hieRat = holdH.(fldLabels{1}).ieRat;
    elseif noGroups == 2
        hieRat = sum([holdH.(fldLabels{1}).ieRat holdH.(fldLabels{2}).ieRat]);
    end
    clear holdH
end

%ie ratio, comparison tests
if noGroups == 3
    if hieRat == 0 %if residuals normal, run anova
        [output.p.ieRatAnova,~,ieRatstats]=anova1([ieRat.(fldLabels{1}) ieRat.(fldLabels{2}) ieRat.(fldLabels{3})],...
            [ones(1,length(ieRat.(fldLabels{1}))) 2.*ones(1,length(ieRat.(fldLabels{2}))) 3.*ones(1,length(ieRat.(fldLabels{3})))]);
        if output.p.ieRatAnova <= 0.05
            output.ieRatPostHocTukey=multcompare(ieRatstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [output.p.ieRatKW,~,~]=kruskalwallis([ieRat.(fldLabels{1}) ieRat.(fldLabels{2}) ieRat.(fldLabels{3})],...
            [ones(1,length(ieRat.(fldLabels{1}))) 2.*ones(1,length(ieRat.(fldLabels{2}))) 3.*ones(1,length(ieRat.(fldLabels{3})))]);
        if output.p.ieRatKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [output.p.ieRatPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(ieRat.(fldLabels{1}),ieRat.(fldLabels{2}));
            [output.p.ieRatPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(ieRat.(fldLabels{1}),ieRat.(fldLabels{3}));
            [output.p.ieRatPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(ieRat.(fldLabels{2}),ieRat.(fldLabels{3}));
            output.pFDR.ieRat=drsFDRpval([output.p.ieRatPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) output.p.ieRatPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) output.p.ieRatPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 0
        if hieRat == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,output.p.ieRatTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(ieRat.(fldLabels{1}),ieRat.(fldLabels{2}));
            end
        elseif hieRat ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [output.p.ieRatMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(ieRat.(fldLabels{1}),ieRat.(fldLabels{2}));
            end
        end
    elseif inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
        if hieRat == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,output.p.ieRatPairedttest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest(ieRat.(fldLabels{1}),ieRat.(fldLabels{2}));
            end
        elseif hieRat ~= 0 %if data not normal, run wilcoxon signed rank test
            if noGroups == 2
                [output.p.ieRatWSR.([fldLabels{1} 'vs' fldLabels{2}]),~]=signrank(ieRat.(fldLabels{1}),ieRat.(fldLabels{2}));
            end
        end
    end
end

if inputs.KineticsAnalysisCheckBox.Value == 1
    %rise time, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.rt.(fldLabels{ii}).epsc = mean(rt.(fldLabels{ii}).epsc) - rt.(fldLabels{ii}).epsc;
            res.rt.(fldLabels{ii}).ipsc = mean(rt.(fldLabels{ii}).ipsc) - rt.(fldLabels{ii}).ipsc;
        end
        allRes.rt.epsc = [res.rt.(fldLabels{1}).epsc res.rt.(fldLabels{2}).epsc res.rt.(fldLabels{3}).epsc];
        allRes.rt.ipsc = [res.rt.(fldLabels{1}).ipsc res.rt.(fldLabels{2}).ipsc res.rt.(fldLabels{3}).ipsc];
        [hRT.epsc,~]=adtest(allRes.rt.epsc);
        [hRT.ipsc,~]=adtest(allRes.rt.ipsc);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}).epsc,~] = adtest(rt.(fldLabels{ii}).epsc);
            [holdH.(fldLabels{ii}).ipsc,~] = adtest(rt.(fldLabels{ii}).ipsc);
        end
        if noGroups == 1
            hRT.epsc = holdH.(fldLabels{1}).epsc;
            hRT.ipsc = holdH.(fldLabels{1}).ipsc;
        elseif noGroups == 2
            hRT.epsc = sum([holdH.(fldLabels{1}).epsc holdH.(fldLabels{2}).epsc]);
            hRT.ipsc = sum([holdH.(fldLabels{1}).ipsc holdH.(fldLabels{2}).ipsc]);
        end
        clear holdH
    end
    
    %rise time, comparison tests
    if noGroups == 3
        if hRT.epsc == 0 %if residuals normal, run anova
            [output.p.rtAnova.epsc,~,rtstats.epsc]=anova1([rt.(fldLabels{1}).epsc rt.(fldLabels{2}).epsc rt.(fldLabels{3}).epsc],...
                [ones(1,length(rt.(fldLabels{1}).epsc)) 2.*ones(1,length(rt.(fldLabels{2}).epsc)) 3.*ones(1,length(rt.(fldLabels{3}).epsc))]);
            if output.p.rtAnova.epsc <= 0.05
                output.rtPostHocTukey.epsc=multcompare(rtstats.epsc); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [output.p.rtKW.epsc,~,~]=kruskalwallis([rt.(fldLabels{1}).epsc rt.(fldLabels{2}).epsc rt.(fldLabels{3}).epsc],...
                [ones(1,length(rt.(fldLabels{1}).epsc)) 2.*ones(1,length(rt.(fldLabels{2}).epsc)) 3.*ones(1,length(rt.(fldLabels{3}).epsc))]);
            if output.p.rtKW.epsc <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=ranksum(rt.(fldLabels{1}).epsc,rt.(fldLabels{2}).epsc);
                [output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).epsc,~]=ranksum(rt.(fldLabels{1}).epsc,rt.(fldLabels{3}).epsc);
                [output.p.rtPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).epsc,~]=ranksum(rt.(fldLabels{2}).epsc,rt.(fldLabels{3}).epsc);
                output.pFDR.rt.epsc=drsFDRpval([output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).epsc output.p.rtPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).epsc]);
            end
        end
        if hRT.ipsc == 0 %if residuals normal, run anova
            [output.p.rtAnova.ipsc,~,rtstats.ipsc]=anova1([rt.(fldLabels{1}).ipsc rt.(fldLabels{2}).ipsc rt.(fldLabels{3}).ipsc],...
                [ones(1,length(rt.(fldLabels{1}).ipsc)) 2.*ones(1,length(rt.(fldLabels{2}).ipsc)) 3.*ones(1,length(rt.(fldLabels{3}).ipsc))]);
            if output.p.rtAnova.ipsc <= 0.05
                output.rtPostHocTukey.ipsc=multcompare(rtstats.ipsc); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [output.p.rtKW.ipsc,~,~]=kruskalwallis([rt.(fldLabels{1}).ipsc rt.(fldLabels{2}).ipsc rt.(fldLabels{3}).ipsc],...
                [ones(1,length(rt.(fldLabels{1}).ipsc)) 2.*ones(1,length(rt.(fldLabels{2}).ipsc)) 3.*ones(1,length(rt.(fldLabels{3}).ipsc))]);
            if output.p.rtKW.ipsc <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=ranksum(rt.(fldLabels{1}).ipsc,rt.(fldLabels{2}).ipsc);
                [output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).ipsc,~]=ranksum(rt.(fldLabels{1}).ipsc,rt.(fldLabels{3}).ipsc);
                [output.p.rtPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).ipsc,~]=ranksum(rt.(fldLabels{2}).ipsc,rt.(fldLabels{3}).ipsc);
                output.pFDR.rt.ipsc=drsFDRpval([output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc output.p.rtPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).ipsc output.p.rtPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).ipsc]);
            end
        end
    elseif noGroups < 3
        if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 0
            if hRT.epsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.rtTtest.([fldLabels{1} 'vs' fldLabels{2}]).epsc]=ttest2(rt.(fldLabels{1}).epsc,rt.(fldLabels{2}).epsc);
                end
            elseif hRT.epsc ~= 0 %if data not normal, run MWU test
                if noGroups == 2
                    [output.p.rtMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=ranksum(rt.(fldLabels{1}).epsc,rt.(fldLabels{2}).epsc);
                end
            end
            if hRT.ipsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.rtTtest.([fldLabels{1} 'vs' fldLabels{2}]).ipsc]=ttest2(rt.(fldLabels{1}).ipsc,rt.(fldLabels{2}).ipsc);
                end
            elseif hRT.ipsc ~= 0 %if data not normal, run MWU test
                if noGroups == 2
                    [output.p.rtMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=ranksum(rt.(fldLabels{1}).ipsc,rt.(fldLabels{2}).ipsc);
                end
            end
        elseif inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
            if hRT.epsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.rtPairedttest.([fldLabels{1} 'vs' fldLabels{2}]).epsc]=ttest(rt.(fldLabels{1}).epsc,rt.(fldLabels{2}).epsc);
                end
            elseif hRT.epsc ~= 0 %if data not normal, run wilcoxon signed rank test
                if noGroups == 2
                    [output.p.rtWSR.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=signrank(rt.(fldLabels{1}).epsc,rt.(fldLabels{2}).epsc);
                end
            end
            if hRT.ipsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.rtPairedttest.([fldLabels{1} 'vs' fldLabels{2}]).ipsc]=ttest(rt.(fldLabels{1}).ipsc,rt.(fldLabels{2}).ipsc);
                end
            elseif hRT.ipsc ~= 0 %if data not normal, run wilcoxon signed rank test
                if noGroups == 2
                    [output.p.rtWSR.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=signrank(rt.(fldLabels{1}).ipsc,rt.(fldLabels{2}).ipsc);
                end
            end
        end
    end
    
    %decay, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.tau.(fldLabels{ii}).epsc = mean(tau.(fldLabels{ii}).epsc) - tau.(fldLabels{ii}).epsc;
            res.tau.(fldLabels{ii}).ipsc = mean(tau.(fldLabels{ii}).ipsc) - tau.(fldLabels{ii}).ipsc;
        end
        allRes.tau.epsc = [res.tau.(fldLabels{1}).epsc res.tau.(fldLabels{2}).epsc res.tau.(fldLabels{3}).epsc];
        allRes.tau.ipsc = [res.tau.(fldLabels{1}).ipsc res.tau.(fldLabels{2}).ipsc res.tau.(fldLabels{3}).ipsc];
        [hTau.epsc,~]=adtest(allRes.tau.epsc);
        [hTau.ipsc,~]=adtest(allRes.tau.ipsc);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}).epsc,~] = adtest(tau.(fldLabels{ii}).epsc);
            [holdH.(fldLabels{ii}).ipsc,~] = adtest(tau.(fldLabels{ii}).ipsc);
        end
        if noGroups == 1
            hTau.epsc = holdH.(fldLabels{1}).epsc;
            hTau.ipsc = holdH.(fldLabels{1}).ipsc;
        elseif noGroups == 2
            hTau.epsc = sum([holdH.(fldLabels{1}).epsc holdH.(fldLabels{2}).epsc]);
            hTau.ipsc = sum([holdH.(fldLabels{1}).ipsc holdH.(fldLabels{2}).ipsc]);
        end
        clear holdH
    end
    
    %decay tau, comparison tests
    if noGroups == 3
        if hTau.epsc == 0 %if residuals normal, run anova
            [output.p.tauAnova.epsc,~,taustats.epsc]=anova1([tau.(fldLabels{1}).epsc tau.(fldLabels{2}).epsc tau.(fldLabels{3}).epsc],...
                [ones(1,length(tau.(fldLabels{1}).epsc)) 2.*ones(1,length(tau.(fldLabels{2}).epsc)) 3.*ones(1,length(tau.(fldLabels{3}).epsc))]);
            if output.p.tauAnova.epsc <= 0.05
                output.tauPostHocTukey.epsc=multcompare(taustats.epsc); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [output.p.tauKW.epsc,~,~]=kruskalwallis([tau.(fldLabels{1}).epsc tau.(fldLabels{2}).epsc tau.(fldLabels{3}).epsc],...
                [ones(1,length(tau.(fldLabels{1}).epsc)) 2.*ones(1,length(tau.(fldLabels{2}).epsc)) 3.*ones(1,length(tau.(fldLabels{3}).epsc))]);
            if output.p.tauKW.epsc <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=ranksum(tau.(fldLabels{1}).epsc,tau.(fldLabels{2}).epsc);
                [output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).epsc,~]=ranksum(tau.(fldLabels{1}).epsc,tau.(fldLabels{3}).epsc);
                [output.p.tauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).epsc,~]=ranksum(tau.(fldLabels{2}).epsc,tau.(fldLabels{3}).epsc);
                output.pFDR.tau.epsc=drsFDRpval([output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).epsc output.p.tauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).epsc]);
            end
        end
        if hTau.ipsc == 0 %if residuals normal, run anova
            [output.p.tauAnova.ipsc,~,taustats.ipsc]=anova1([tau.(fldLabels{1}).ipsc' tau.(fldLabels{2}).ipsc' tau.(fldLabels{3}).ipsc'],...
                [ones(1,length(tau.(fldLabels{1}).ipsc)) 2.*ones(1,length(tau.(fldLabels{2}).ipsc)) 3.*ones(1,length(tau.(fldLabels{3}).ipsc))]);
            if output.p.tauAnova.ipsc <= 0.05
                output.tauPostHocTukey.ipsc=multcompare(taustats.ipsc); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [output.p.tauKW.ipsc,~,~]=kruskalwallis([tau.(fldLabels{1}).ipsc tau.(fldLabels{2}).ipsc tau.(fldLabels{3}).ipsc],...
                [ones(1,length(tau.(fldLabels{1}).ipsc)) 2.*ones(1,length(tau.(fldLabels{2}).ipsc)) 3.*ones(1,length(tau.(fldLabels{3}).ipsc))]);
            if output.p.tauKW.ipsc <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=ranksum(tau.(fldLabels{1}).ipsc,tau.(fldLabels{2}).ipsc);
                [output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).ipsc,~]=ranksum(tau.(fldLabels{1}).ipsc,tau.(fldLabels{3}).ipsc);
                [output.p.tauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).ipsc,~]=ranksum(tau.(fldLabels{2}).ipsc,tau.(fldLabels{3}).ipsc);
                output.pFDR.tau.ipsc=drsFDRpval([output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc output.p.tauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]).ipsc output.p.tauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]).ipsc]);
            end
        end
    elseif noGroups < 3
        if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 0
            if hTau.epsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.tauTtest.([fldLabels{1} 'vs' fldLabels{2}]).epsc]=ttest2(tau.(fldLabels{1}).epsc,tau.(fldLabels{2}).epsc);
                end
            elseif hTau.epsc ~= 0 %if data not normal, run MWU test
                if noGroups == 2
                    [output.p.tauMW.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=ranksum(tau.(fldLabels{1}).epsc,tau.(fldLabels{2}).epsc);
                end
            end
            if hTau.ipsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.tauTtest.([fldLabels{1} 'vs' fldLabels{2}]).ipsc]=ttest2(tau.(fldLabels{1}).ipsc,tau.(fldLabels{2}).ipsc);
                end
            elseif hTau.ipsc ~= 0 %if data not normal, run MWU test
                if noGroups == 2
                    [output.p.tauMW.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=ranksum(tau.(fldLabels{1}).ipsc,tau.(fldLabels{2}).ipsc);
                end
            end
        elseif inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
            if hTau.epsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.tauPairedttest.([fldLabels{1} 'vs' fldLabels{2}]).epsc]=ttest(tau.(fldLabels{1}).epsc,tau.(fldLabels{2}).epsc);
                end
            elseif hTau.epsc ~= 0 %if data not normal, run wilcoxon signed rank test
                if noGroups == 2
                    [output.p.tauWSR.([fldLabels{1} 'vs' fldLabels{2}]).epsc,~]=signrank(tau.(fldLabels{1}).epsc,tau.(fldLabels{2}).epsc);
                end
            end
            if hTau.ipsc == 0 %if groups normal, run ttest
                if noGroups == 2
                    [~,output.p.tauPairedttest.([fldLabels{1} 'vs' fldLabels{2}]).ipsc]=ttest(tau.(fldLabels{1}).ipsc,tau.(fldLabels{2}).ipsc);
                end
            elseif hTau.ipsc ~= 0 %if data not normal, run wilcoxon signed rank test
                if noGroups == 2
                    [output.p.tauWSR.([fldLabels{1} 'vs' fldLabels{2}]).ipsc,~]=signrank(tau.(fldLabels{1}).ipsc,tau.(fldLabels{2}).ipsc);
                end
            end
        end
    end
end

%latency, normality tests
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            [holdH.(fldLabels{ii}).epsc,~] = adtest(lat.(fldLabels{ii}).epsc);
            [holdH.(fldLabels{ii}).ipsc,~] = adtest(lat.(fldLabels{ii}).ipsc);
            hLat.(fldLabels{ii}) = sum([holdH.(fldLabels{ii}).epsc holdH.(fldLabels{ii}).ipsc]);
        end
    else
        [holdH.(fldLabels{ii}).epsc,~] = adtest(lat.(fldLabels{ii}).epsc);
        [holdH.(fldLabels{ii}).ipsc,~] = adtest(lat.(fldLabels{ii}).ipsc);
        hLat.(fldLabels{ii}) = sum([holdH.(fldLabels{ii}).epsc holdH.(fldLabels{ii}).ipsc]);
    end
end
clear holdH

%latency, comparison tests
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            if hLat.(fldLabels{ii}) == 0 %if groups normal, run ttest
                [~,output.p.latPairedttest.(fldLabels{ii})]=ttest(lat.(fldLabels{ii}).epsc,lat.(fldLabels{ii}).ipsc);
            elseif hLat.(fldLabels{ii}) ~=0 %if data not normal, run wilcoxon signed rank test
                [output.p.latWSR.(fldLabels{ii}),~]=signrank(lat.(fldLabels{ii}).epsc,lat.(fldLabels{ii}).ipsc);
            end
        end
    else
        if hLat.(fldLabels{ii}) == 0 %if groups normal, run ttest
            [~,output.p.latPairedttest.(fldLabels{ii})]=ttest(lat.(fldLabels{ii}).epsc,lat.(fldLabels{ii}).ipsc);
        elseif hLat.(fldLabels{ii}) ~=0 %if data not normal, run wilcoxon signed rank test
            [output.p.latWSR.(fldLabels{ii}),~]=signrank(lat.(fldLabels{ii}).epsc,lat.(fldLabels{ii}).ipsc);
        end
    end
end

%jitter, normality tests
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            [holdH.(fldLabels{ii}).epsc,~] = adtest(jitter.(fldLabels{ii}).epsc);
            [holdH.(fldLabels{ii}).ipsc,~] = adtest(jitter.(fldLabels{ii}).ipsc);
            hJtr.(fldLabels{ii}) = sum([holdH.(fldLabels{ii}).epsc holdH.(fldLabels{ii}).ipsc]);
        end
    else
        [holdH.(fldLabels{ii}).epsc,~] = adtest(jitter.(fldLabels{ii}).epsc);
        [holdH.(fldLabels{ii}).ipsc,~] = adtest(jitter.(fldLabels{ii}).ipsc);
        hJtr.(fldLabels{ii}) = sum([holdH.(fldLabels{ii}).epsc holdH.(fldLabels{ii}).ipsc]);
    end
end
clear holdH

%jitter, comparison tests
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            if hJtr.(fldLabels{ii}) == 0 %if groups normal, run ttest
                [~,output.p.jtrPairedttest.(fldLabels{ii})]=ttest(jitter.(fldLabels{ii}).epsc,jitter.(fldLabels{ii}).ipsc);
            elseif hJtr.(fldLabels{ii}) ~=0 %if data not normal, run wilcoxon signed rank test
                [output.p.jtrWSR.(fldLabels{ii}),~]=signrank(jitter.(fldLabels{ii}).epsc,jitter.(fldLabels{ii}).ipsc);
            end
        end
    else
        if hJtr.(fldLabels{ii}) == 0 %if groups normal, run ttest
            [~,output.p.jtrPairedttest.(fldLabels{ii})]=ttest(jitter.(fldLabels{ii}).epsc,jitter.(fldLabels{ii}).ipsc);
        elseif hJtr.(fldLabels{ii}) ~=0 %if data not normal, run wilcoxon signed rank test
            [output.p.jtrWSR.(fldLabels{ii}),~]=signrank(jitter.(fldLabels{ii}).epsc,jitter.(fldLabels{ii}).ipsc);
        end
    end
end

close all

%% DISPLAY DATA
%color
lightgray=[.75 .75 .75]; %light gray for individual data points
darkgray=[.33 .33 .33]; %dark gray for connecting mean/medians in paired analysis
lightColorMax = (1+(.85-max(inputs.plotColor))./max(inputs.plotColor)).*inputs.plotColor; %linearly scales RGB values to a point where max value is .85 for one of R,G,B
pColor = -1.*ones(noGroups,3);
for ii = 1:noGroups
    if ii == 1
        pColor(ii,:) = inputs.plotColor;
    elseif ii == 2
        if noGroups == 2
            pColor(ii,:) = lightColorMax;
        else
            pColor(ii,:) = .5.*(pColor(1,:) + lightColorMax);
        end
    elseif ii == 3
        pColor(ii,:) = lightColorMax;
    end
end

%amplitude figure, epsc
ampFig.epsc=figure(noGroups*2+1);
ampFig.epsc.Position = [415 545 165 250];
hold on
line([0 noGroups+1],[0 0],'linewidth',1,'color','k','linestyle','--')
%plot connecting lines for paired analysis display
if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
    for ii = 1:length(amp.(fldLabels{1}).epsc)
        plot([1 2],[amp.(fldLabels{1}).epsc(ii),amp.(fldLabels{2}).epsc(ii)],'color',lightgray,'linewidth',1)
    end
    if hAmp.epsc == 0 %if data are normal, plot mean
        plot([1 2],[mean(amp.(fldLabels{1}).epsc),mean(amp.(fldLabels{2}).epsc)],'color',darkgray,'linewidth',2)
    elseif hAmp.epsc ~=0 %if data non-normal, plot median
        plot([1 2],[median(amp.(fldLabels{1}).epsc),median(amp.(fldLabels{2}).epsc)],'color',darkgray,'linewidth',2)
    end
end
for ii = 1:noGroups
    scatter(ii.*ones(1,length(amp.(fldLabels{ii}).epsc)),amp.(fldLabels{ii}).epsc,40,lightgray,'filled')
    if hAmp.epsc == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
        errorbar(ii,mean(amp.(fldLabels{ii}).epsc),sem(amp.(fldLabels{ii}).epsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(amp.(fldLabels{ii}).epsc),125,pColor(ii,:),'filled')
    elseif hAmp.epsc ~= 0 %if data non-normal, plot median w IQR
        [lowerError.amp.(fldLabels{ii}).epsc,upperError.amp.(fldLabels{ii}).epsc] = iqrError(amp.(fldLabels{ii}).epsc,2);
        errorbar(ii,median(amp.(fldLabels{ii}).epsc),lowerError.amp.(fldLabels{ii}).epsc,upperError.amp.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(amp.(fldLabels{ii}).epsc),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allAmp.epsc = amp.(fldLabels{ii}).epsc;
    else
        allAmp.epsc = [allAmp.epsc amp.(fldLabels{ii}).epsc];
    end
end
ylabel('EPSC amplitude (pA)')
ampAx.epsc=gca;
setAx(ampAx.epsc);
if noGroups == 1
    ampAx.epsc.XTick = 1;
    ampAx.epsc.XTickLabel = fldLabels{1};
elseif noGroups == 2
    ampAx.epsc.XTick = [1 2];
    ampAx.epsc.XTickLabel = {fldLabels{1} fldLabels{2}};
elseif noGroups == 3
    ampAx.epsc.XTick = [1 2 3];
    ampAx.epsc.XTickLabel = {fldLabels{1} fldLabels{2} fldLabels{3}};
end
xlim([0.5 noGroups+.5])
ylim([floor(1.5*min(allAmp.epsc))-.25*min(allAmp.epsc) abs(.1*(floor(1.5*min(allAmp.epsc))-.25*min(allAmp.epsc)))])


%amplitude figure, ipsc
ampFig.ipsc=figure(noGroups*2+2);
ampFig.ipsc.Position = [595 545 165 250];
hold on
line([0 noGroups+1],[0 0],'linewidth',1,'color','k','linestyle','--')
%plot connecting lines for paired analysis display
if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
    for ii = 1:length(amp.(fldLabels{1}).ipsc)
        plot([1 2],[amp.(fldLabels{1}).ipsc(ii),amp.(fldLabels{2}).ipsc(ii)],'color',lightgray,'linewidth',1)
    end
    if hAmp.ipsc == 0 %if data are normal, plot mean
        plot([1 2],[mean(amp.(fldLabels{1}).ipsc),mean(amp.(fldLabels{2}).ipsc)],'color',darkgray,'linewidth',2)
    elseif hAmp.ipsc ~=0 %if data non-normal, plot median
        plot([1 2],[median(amp.(fldLabels{1}).ipsc),median(amp.(fldLabels{2}).ipsc)],'color',darkgray,'linewidth',2)
    end
end
for ii = 1:noGroups
    scatter(ii.*ones(1,length(amp.(fldLabels{ii}).ipsc)),amp.(fldLabels{ii}).ipsc,40,lightgray,'filled')
    if hAmp.ipsc == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(amp.(fldLabels{ii}).ipsc),sem(amp.(fldLabels{ii}).ipsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(amp.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
    elseif hAmp.ipsc ~= 0 %if data non-normal, plot median w IQR
        [lowerError.amp.(fldLabels{ii}).ipsc,upperError.amp.(fldLabels{ii}).ipsc] = iqrError(amp.(fldLabels{ii}).ipsc,2);
        errorbar(ii,median(amp.(fldLabels{ii}).ipsc),lowerError.amp.(fldLabels{ii}).ipsc,upperError.amp.(fldLabels{ii}).ipsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(amp.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allAmp.ipsc = amp.(fldLabels{ii}).ipsc;
    else
        allAmp.ipsc = [allAmp.ipsc amp.(fldLabels{ii}).ipsc];
    end
end
ylabel('IPSC amplitude (pA)')
ampAx.ipsc=gca;
setAx(ampAx.ipsc);
if noGroups == 1
    ampAx.ipsc.XTick = 1;
    ampAx.ipsc.XTickLabel = fldLabels{1};
elseif noGroups == 2
    ampAx.ipsc.XTick = [1 2];
    ampAx.ipsc.XTickLabel = {fldLabels{1} fldLabels{2}};
elseif noGroups == 3
    ampAx.ipsc.XTick = [1 2 3];
    ampAx.ipsc.XTickLabel = {fldLabels{1} fldLabels{2} fldLabels{3}};
end
xlim([0.5 noGroups+.5])
ylim([round(-.1*ceil(1.1*max(allAmp.ipsc))) ceil(1.1*max(allAmp.ipsc))])

%ie ratio figure
ieRatFig=figure(noGroups*2+3);
ieRatFig.Position = [775 545 165 250];
hold on
%plot connecting lines for paired analysis display
if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
    for ii = 1:length(ieRat.(fldLabels{1}))
        plot([1 2],[ieRat.(fldLabels{1})(ii),ieRat.(fldLabels{2})(ii)],'color',lightgray,'linewidth',1)
    end
    if hieRat == 0 %if data are normal, plot mean
        plot([1 2],[mean(ieRat.(fldLabels{1})),mean(ieRat.(fldLabels{2}))],'color',darkgray,'linewidth',2)
    elseif hieRat ~=0 %if data non-normal, plot median
        plot([1 2],[median(ieRat.(fldLabels{1})),median(ieRat.(fldLabels{2}))],'color',darkgray,'linewidth',2)
    end
end
for ii = 1:noGroups
    scatter(ii.*ones(1,length(ieRat.(fldLabels{ii}))),ieRat.(fldLabels{ii}),40,lightgray,'filled')
    if hieRat == 0 %if data are normal, plot mean w sem
        errorbar(ii,mean(ieRat.(fldLabels{ii})),sem(ieRat.(fldLabels{ii}),2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(ieRat.(fldLabels{ii})),125,pColor(ii,:),'filled')
    elseif hieRat ~= 0 %if data non-normal, plot median w IQR
        [lowerError.ieRat.(fldLabels{ii}),upperError.ieRat.(fldLabels{ii})] = iqrError(ieRat.(fldLabels{ii}),2);
        errorbar(ii,median(ieRat.(fldLabels{ii})),lowerError.ieRat.(fldLabels{ii}),upperError.ieRat.(fldLabels{ii}),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(ieRat.(fldLabels{ii})),125,pColor(ii,:),'filled')
    end
    if ii == 1
        allie= ieRat.(fldLabels{ii});
    else
        allie= [allie ieRat.(fldLabels{ii})];
    end
end
ylabel('Feedforward I/E Ratio')
ieAx=gca;
setAx(ieAx);
if noGroups == 1
    ieAx.XTick = 1;
    ieAx.XTickLabel = fldLabels{1};
elseif noGroups == 2
    ieAx.XTick = [1 2];
    ieAx.XTickLabel = {fldLabels{1} fldLabels{2}};
elseif noGroups == 3
    ieAx.XTick = [1 2 3];
    ieAx.XTickLabel = {fldLabels{1} fldLabels{2} fldLabels{3}};
end
xlim([0.5 noGroups+.5])
ylim([floor(1.25*min(allie))-.25*min(allie) ceil(1.1*max(allie))])

if inputs.KineticsAnalysisCheckBox.Value == 1
    %risetime figure, epsc
    rtFig.epsc=figure(noGroups*2+4);
    rtFig.epsc.Position = [115 65 165 250];
    hold on
    %plot connecting lines for paired analysis display
    if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
        for ii = 1:length(rt.(fldLabels{1}).epsc)
            plot([1 2],[rt.(fldLabels{1}).epsc(ii),rt.(fldLabels{2}).epsc(ii)],'color',lightgray,'linewidth',1)
        end
        if hRT.epsc == 0 %if data are normal, plot mean
            plot([1 2],[mean(rt.(fldLabels{1}).epsc),mean(rt.(fldLabels{2}).epsc)],'color',darkgray,'linewidth',2)
        elseif hRT.epsc ~=0 %if data non-normal, plot median
            plot([1 2],[median(rt.(fldLabels{1}).epsc),median(rt.(fldLabels{2}).epsc)],'color',darkgray,'linewidth',2)
        end
    end
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(rt.(fldLabels{ii}).epsc)),rt.(fldLabels{ii}).epsc,40,lightgray,'filled')
        if hRT.epsc == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
            errorbar(ii,mean(rt.(fldLabels{ii}).epsc),sem(rt.(fldLabels{ii}).epsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(rt.(fldLabels{ii}).epsc),125,pColor(ii,:),'filled')
        elseif hRT.epsc ~= 0 %if data non-normal, plot median w IQR
            [lowerError.rt.(fldLabels{ii}).epsc,upperError.rt.(fldLabels{ii}).epsc] = iqrError(rt.(fldLabels{ii}).epsc,2);
            errorbar(ii,median(rt.(fldLabels{ii}).epsc),lowerError.rt.(fldLabels{ii}).epsc,upperError.rt.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(rt.(fldLabels{ii}).epsc),125,pColor(ii,:),'filled')
        end
        if ii == 1
            allRT.epsc = rt.(fldLabels{ii}).epsc;
        else
            allRT.epsc = [allRT.epsc rt.(fldLabels{ii}).epsc];
        end
    end
    ylabel('EPSC 20%-80% risetime (ms)')
    rtAx.epsc=gca;
    setAx(rtAx.epsc);
    if noGroups == 1
        rtAx.epsc.XTick = 1;
        rtAx.epsc.XTickLabel = fldLabels{1};
    elseif noGroups == 2
        rtAx.epsc.XTick = [1 2];
        rtAx.epsc.XTickLabel = {fldLabels{1} fldLabels{2}};
    elseif noGroups == 3
        rtAx.epsc.XTick = [1 2 3];
        rtAx.epsc.XTickLabel = {fldLabels{1} fldLabels{2} fldLabels{3}};
    end
    xlim([0.5 noGroups+.5])
    ylim([0 ceil(1.1*max(allRT.epsc))])
    
    %risetime figure, ipsc
    rtFig.ipsc=figure(noGroups*2+5);
    rtFig.ipsc.Position = [295 65 165 250];
    hold on
    %plot connecting lines for paired analysis display
    if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
        for ii = 1:length(rt.(fldLabels{1}).ipsc)
            plot([1 2],[rt.(fldLabels{1}).ipsc(ii),rt.(fldLabels{2}).ipsc(ii)],'color',lightgray,'linewidth',1)
        end
        if hRT.ipsc == 0 %if data are normal, plot mean
            plot([1 2],[mean(rt.(fldLabels{1}).ipsc),mean(rt.(fldLabels{2}).ipsc)],'color',darkgray,'linewidth',2)
        elseif hRT.ipsc ~=0 %if data non-normal, plot median
            plot([1 2],[median(rt.(fldLabels{1}).ipsc),median(rt.(fldLabels{2}).ipsc)],'color',darkgray,'linewidth',2)
        end
    end
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(rt.(fldLabels{ii}).ipsc)),rt.(fldLabels{ii}).ipsc,40,lightgray,'filled')
        if hRT.ipsc == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(rt.(fldLabels{ii}).ipsc),sem(rt.(fldLabels{ii}).ipsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(rt.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        elseif hRT.ipsc ~= 0 %if data non-normal, plot median w IQR
            [lowerError.rt.(fldLabels{ii}).ipsc,upperError.rt.(fldLabels{ii}).ipsc] = iqrError(rt.(fldLabels{ii}).ipsc,2);
            errorbar(ii,median(rt.(fldLabels{ii}).ipsc),lowerError.rt.(fldLabels{ii}).ipsc,upperError.rt.(fldLabels{ii}).ipsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(rt.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        end
        if ii == 1
            allRT.ipsc = rt.(fldLabels{ii}).ipsc;
        else
            allRT.ipsc = [allRT.ipsc rt.(fldLabels{ii}).ipsc];
        end
    end
    ylabel('IPSC 20%-80% risetime (ms)')
    rtAx.ipsc=gca;
    setAx(rtAx.ipsc);
    if noGroups == 1
        rtAx.ipsc.XTick = 1;
        rtAx.ipsc.XTickLabel = fldLabels{1};
    elseif noGroups == 2
        rtAx.ipsc.XTick = [1 2];
        rtAx.ipsc.XTickLabel = {fldLabels{1} fldLabels{2}};
    elseif noGroups == 3
        rtAx.ipsc.XTick = [1 2 3];
        rtAx.ipsc.XTickLabel = {fldLabels{1} fldLabels{2} fldLabels{3}};
    end
    xlim([0.5 noGroups+.5])
    ylim([0 ceil(1.1*max(allRT.ipsc))])
    
    %decay tau figure, epsc
    tauFig.epsc=figure(noGroups*2+6);
    tauFig.epsc.Position = [475 65 165 250];
    hold on
    %plot connecting lines for paired analysis display
    if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
        for ii = 1:length(tau.(fldLabels{1}).epsc)
            plot([1 2],[tau.(fldLabels{1}).epsc(ii),tau.(fldLabels{2}).epsc(ii)],'color',lightgray,'linewidth',1)
        end
        if hTau.epsc == 0 %if data are normal, plot mean
            plot([1 2],[mean(tau.(fldLabels{1}).epsc),mean(tau.(fldLabels{2}).epsc)],'color',darkgray,'linewidth',2)
        elseif hTau.epsc ~=0 %if data non-normal, plot median
            plot([1 2],[median(tau.(fldLabels{1}).epsc),median(tau.(fldLabels{2}).epsc)],'color',darkgray,'linewidth',2)
        end
    end
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(tau.(fldLabels{ii}).epsc)),tau.(fldLabels{ii}).epsc,40,lightgray,'filled')
        if hTau.epsc == 0 %if data are normal, plot mean w sem (dir = 1 for sem if colm vector, 2 if row vector)
            errorbar(ii,mean(tau.(fldLabels{ii}).epsc),sem(tau.(fldLabels{ii}).epsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(tau.(fldLabels{ii}).epsc),125,pColor(ii,:),'filled')
        elseif hTau.epsc ~= 0 %if data non-normal, plot median w IQR
            [lowerError.tau.(fldLabels{ii}).epsc,upperError.tau.(fldLabels{ii}).epsc] = iqrError(tau.(fldLabels{ii}).epsc,2);
            errorbar(ii,median(tau.(fldLabels{ii}).epsc),lowerError.tau.(fldLabels{ii}).epsc,upperError.tau.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(tau.(fldLabels{ii}).epsc),125,pColor(ii,:),'filled')
        end
        if ii == 1
            allTau.epsc = tau.(fldLabels{ii}).epsc;
        else
            allTau.epsc = [allTau.epsc tau.(fldLabels{ii}).epsc];
        end
    end
    ylabel('EPSC decay \tau (ms)')
    tauAx.epsc=gca;
    setAx(tauAx.epsc);
    if noGroups == 1
        tauAx.epsc.XTick = 1;
        tauAx.epsc.XTickLabel = fldLabels{1};
    elseif noGroups == 2
        tauAx.epsc.XTick = [1 2];
        tauAx.epsc.XTickLabel = {fldLabels{1} fldLabels{2}};
    elseif noGroups == 3
        tauAx.epsc.XTick = [1 2 3];
        tauAx.epsc.XTickLabel = {fldLabels{1} fldLabels{2} fldLabels{3}};
    end
    xlim([0.5 noGroups+.5])
    ylim([0 ceil(1.1*max(allTau.epsc))])
    
    %decay tau figure, ipsc
    tauFig.ipsc=figure(noGroups*2+7);
    tauFig.ipsc.Position = [655 65 165 250];
    hold on
    %plot connecting lines for paired analysis display
    if inputs.PairedAnalysisTwogroupsonlyCheckBox.Value == 1
        for ii = 1:length(tau.(fldLabels{1}).ipsc)
            plot([1 2],[tau.(fldLabels{1}).ipsc(ii),tau.(fldLabels{2}).ipsc(ii)],'color',lightgray,'linewidth',1)
        end
        if hTau.ipsc == 0 %if data are normal, plot mean
            plot([1 2],[mean(tau.(fldLabels{1}).ipsc),mean(tau.(fldLabels{2}).ipsc)],'color',darkgray,'linewidth',2)
        elseif hTau.ipsc ~=0 %if data non-normal, plot median
            plot([1 2],[median(tau.(fldLabels{1}).ipsc),median(tau.(fldLabels{2}).ipsc)],'color',darkgray,'linewidth',2)
        end
    end
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(tau.(fldLabels{ii}).ipsc)),tau.(fldLabels{ii}).ipsc,40,lightgray,'filled')
        if hTau.ipsc == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(tau.(fldLabels{ii}).ipsc),sem(tau.(fldLabels{ii}).ipsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(tau.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        elseif hTau.ipsc ~= 0 %if data non-normal, plot median w IQR
            [lowerError.tau.(fldLabels{ii}).ipsc,upperError.tau.(fldLabels{ii}).ipsc] = iqrError(tau.(fldLabels{ii}).ipsc,2);
            errorbar(ii,median(tau.(fldLabels{ii}).ipsc),lowerError.tau.(fldLabels{ii}).ipsc,upperError.tau.(fldLabels{ii}).ipsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(tau.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        end
        if ii == 1
            allTau.ipsc = tau.(fldLabels{ii}).ipsc;
        else
            allTau.ipsc = [allTau.ipsc tau.(fldLabels{ii}).ipsc];
        end
    end
    ylabel('IPSC decay \tau (ms)')
    tauAx.ipsc=gca;
    setAx(tauAx.ipsc);
    if noGroups == 1
        tauAx.ipsc.XTick = 1;
        tauAx.ipsc.XTickLabel = fldLabels{1};
    elseif noGroups == 2
        tauAx.ipsc.XTick = [1 2];
        tauAx.ipsc.XTickLabel = {fldLabels{1} fldLabels{2}};
    elseif noGroups == 3
        tauAx.ipsc.XTick = [1 2 3];
        tauAx.ipsc.XTickLabel = {fldLabels{1} fldLabels{2} fldLabels{3}};
    end
    xlim([0.5 noGroups+.5])
    ylim([0 ceil(1.1*max(allTau.ipsc))])
end

%latency figures
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            latFig.(fldLabels{ii})=figure(noGroups*2+7+ii);
            latFig.(fldLabels{ii}).Position = [775+ii*185 545 165 250];
            hold on
            %plot connecting lines for paired analysis display
            for jj = 1:length(lat.(fldLabels{ii}).epsc)
                plot([1 2],[lat.(fldLabels{ii}).epsc(jj),lat.(fldLabels{ii}).ipsc(jj)],'color',lightgray,'linewidth',1)
                if hLat.(fldLabels{ii}) == 0 %if data are normal, plot mean
                    plot([1 2],[mean(lat.(fldLabels{ii}).epsc),mean(lat.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
                elseif hLat.(fldLabels{ii}) ~=0 %if data non-normal, plot median
                    plot([1 2],[median(lat.(fldLabels{ii}).epsc),median(lat.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
                end
            end
            scatter(ones(1,length(lat.(fldLabels{ii}).epsc)),lat.(fldLabels{ii}).epsc,40,lightgray,'filled')
            scatter(2.*ones(1,length(lat.(fldLabels{ii}).ipsc)),lat.(fldLabels{ii}).ipsc,40,lightgray,'filled')
            if hLat.(fldLabels{ii}) == 0 %if data are normal, plot mean w sem
                errorbar(1,mean(lat.(fldLabels{ii}).epsc),sem(lat.(fldLabels{ii}).epsc,2),'color',darkgray,'linewidth',3,'CapSize',0)
                scatter(1,mean(lat.(fldLabels{ii}).epsc),125,darkgray,'filled')
                errorbar(2,mean(lat.(fldLabels{ii}).ipsc),sem(lat.(fldLabels{ii}).ipsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
                scatter(2,mean(lat.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
            elseif hLat.(fldLabels{ii}) ~= 0 %if data non-normal, plot median w IQR
                [lowerError.lat.(fldLabels{ii}).epsc,upperError.lat.(fldLabels{ii}).epsc] = iqrError(lat.(fldLabels{ii}).epsc,2);
                [lowerError.lat.(fldLabels{ii}).ipsc,upperError.lat.(fldLabels{ii}).ipsc] = iqrError(lat.(fldLabels{ii}).ipsc,2);
                errorbar(1,median(lat.(fldLabels{ii}).epsc),lowerError.lat.(fldLabels{ii}).epsc,upperError.lat.(fldLabels{ii}).epsc,'color',darkgray,'linewidth',3,'CapSize',0)
                scatter(1,median(lat.(fldLabels{ii}).epsc),125,darkgray,'filled')
                errorbar(2,median(lat.(fldLabels{ii}).ipsc),lowerError.lat.(fldLabels{ii}).ipsc,upperError.lat.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
                scatter(2,median(lat.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
            end
            if ii == 1
                allLat = [lat.(fldLabels{ii}).epsc lat.(fldLabels{ii}).ipsc];
            else
                allLat = [allLat lat.(fldLabels{ii}).epsc lat.(fldLabels{ii}).ipsc];
            end
            ylabel(['Latency, ',fldLabels{ii},' group'])
            latAx.(fldLabels{ii})=gca;
            setAx(latAx.(fldLabels{ii}));
            latAx.(fldLabels{ii}).XTick = [1 2];
            latAx.(fldLabels{ii}).XTickLabel = {'epsc' 'ipsc'};
            xlim([0.5 2.5])
        end
    else
        latFig.(fldLabels{ii})=figure(noGroups*2+7+ii);
        latFig.(fldLabels{ii}).Position = [775+ii*185 545 165 250];
        hold on
        %plot connecting lines for paired analysis display
        for jj = 1:length(lat.(fldLabels{ii}).epsc)
            plot([1 2],[lat.(fldLabels{ii}).epsc(jj),lat.(fldLabels{ii}).ipsc(jj)],'color',lightgray,'linewidth',1)
            if hLat.(fldLabels{ii}) == 0 %if data are normal, plot mean
                plot([1 2],[mean(lat.(fldLabels{ii}).epsc),mean(lat.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
            elseif hLat.(fldLabels{ii}) ~=0 %if data non-normal, plot median
                plot([1 2],[median(lat.(fldLabels{ii}).epsc),median(lat.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
            end
        end
        scatter(ones(1,length(lat.(fldLabels{ii}).epsc)),lat.(fldLabels{ii}).epsc,40,lightgray,'filled')
        scatter(2.*ones(1,length(lat.(fldLabels{ii}).ipsc)),lat.(fldLabels{ii}).ipsc,40,lightgray,'filled')
        if hLat.(fldLabels{ii}) == 0 %if data are normal, plot mean w sem
            errorbar(1,mean(lat.(fldLabels{ii}).epsc),sem(lat.(fldLabels{ii}).epsc,2),'color',darkgray,'linewidth',3,'CapSize',0)
            scatter(1,mean(lat.(fldLabels{ii}).epsc),125,darkgray,'filled')
            errorbar(2,mean(lat.(fldLabels{ii}).ipsc),sem(lat.(fldLabels{ii}).ipsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(2,mean(lat.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        elseif hLat.(fldLabels{ii}) ~= 0 %if data non-normal, plot median w IQR
            [lowerError.lat.(fldLabels{ii}).epsc,upperError.lat.(fldLabels{ii}).epsc] = iqrError(lat.(fldLabels{ii}).epsc,2);
            [lowerError.lat.(fldLabels{ii}).ipsc,upperError.lat.(fldLabels{ii}).ipsc] = iqrError(lat.(fldLabels{ii}).ipsc,2);
            errorbar(1,median(lat.(fldLabels{ii}).epsc),lowerError.lat.(fldLabels{ii}).epsc,upperError.lat.(fldLabels{ii}).epsc,'color',darkgray,'linewidth',3,'CapSize',0)
            scatter(1,median(lat.(fldLabels{ii}).epsc),125,darkgray,'filled')
            errorbar(2,median(lat.(fldLabels{ii}).ipsc),lowerError.lat.(fldLabels{ii}).ipsc,upperError.lat.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(2,median(lat.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        end
        if ii == 1
            allLat = [lat.(fldLabels{ii}).epsc lat.(fldLabels{ii}).ipsc];
        else
            allLat = [allLat lat.(fldLabels{ii}).epsc lat.(fldLabels{ii}).ipsc];
        end
        ylabel(['Latency, ',fldLabels{ii},' group'])
        latAx.(fldLabels{ii})=gca;
        setAx(latAx.(fldLabels{ii}));
        latAx.(fldLabels{ii}).XTick = [1 2];
        latAx.(fldLabels{ii}).XTickLabel = {'epsc' 'ipsc'};
        xlim([0.5 2.5])
    end
end
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            thisFig = latFig.(fldLabels{ii}).Number;
            figure(thisFig);
            ylim([0 ceil(1.1*max(allLat))])
        end
    else
        thisFig = latFig.(fldLabels{ii}).Number;
        figure(thisFig);
        ylim([0 ceil(1.1*max(allLat))])
    end
end

%jitter figures
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            jtrFig.(fldLabels{ii})=figure(noGroups*3+7+ii);
            jtrFig.(fldLabels{ii}).Position = [1450 545 165 250];
            hold on
            %plot connecting lines for paired analysis display
            for jj = 1:length(jitter.(fldLabels{ii}).epsc)
                plot([1 2],[jitter.(fldLabels{ii}).epsc(jj),jitter.(fldLabels{ii}).ipsc(jj)],'color',lightgray,'linewidth',1)
                if hJtr.(fldLabels{ii}) == 0 %if data are normal, plot mean
                    plot([1 2],[mean(jitter.(fldLabels{ii}).epsc),mean(jitter.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
                elseif hJtr.(fldLabels{ii}) ~=0 %if data non-normal, plot median
                    plot([1 2],[median(jitter.(fldLabels{ii}).epsc),median(jitter.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
                end
            end
            scatter(ones(1,length(jitter.(fldLabels{ii}).epsc)),jitter.(fldLabels{ii}).epsc,40,lightgray,'filled')
            scatter(2.*ones(1,length(jitter.(fldLabels{ii}).ipsc)),jitter.(fldLabels{ii}).ipsc,40,lightgray,'filled')
            if hJtr.(fldLabels{ii}) == 0 %if data are normal, plot mean w sem
                errorbar(1,mean(jitter.(fldLabels{ii}).epsc),sem(jitter.(fldLabels{ii}).epsc,2),'color',darkgray,'linewidth',3,'CapSize',0)
                scatter(1,mean(jitter.(fldLabels{ii}).epsc),125,darkgray,'filled')
                errorbar(2,mean(jitter.(fldLabels{ii}).ipsc),sem(jitter.(fldLabels{ii}).ipsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
                scatter(2,mean(jitter.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
            elseif hJtr.(fldLabels{ii}) ~= 0 %if data non-normal, plot median w IQR
                [lowerError.jitter.(fldLabels{ii}).epsc,upperError.jitter.(fldLabels{ii}).epsc] = iqrError(jitter.(fldLabels{ii}).epsc,2);
                [lowerError.jitter.(fldLabels{ii}).ipsc,upperError.jitter.(fldLabels{ii}).ipsc] = iqrError(jitter.(fldLabels{ii}).ipsc,2);
                errorbar(1,median(jitter.(fldLabels{ii}).epsc),lowerError.jitter.(fldLabels{ii}).epsc,upperError.jitter.(fldLabels{ii}).epsc,'color',darkgray,'linewidth',3,'CapSize',0)
                scatter(1,median(jitter.(fldLabels{ii}).epsc),125,darkgray,'filled')
                errorbar(2,median(jitter.(fldLabels{ii}).ipsc),lowerError.jitter.(fldLabels{ii}).ipsc,upperError.jitter.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
                scatter(2,median(jitter.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
            end
            if ii == 1
                allJtr = [jitter.(fldLabels{ii}).epsc jitter.(fldLabels{ii}).ipsc];
            else
                allJtr = [allJtr jitter.(fldLabels{ii}).epsc jitter.(fldLabels{ii}).ipsc];
            end
            ylabel(['Jitter, ',fldLabels{ii},' group'])
            jitterAx.(fldLabels{ii})=gca;
            setAx(jitterAx.(fldLabels{ii}));
            jitterAx.(fldLabels{ii}).XTick = [1 2];
            jitterAx.(fldLabels{ii}).XTickLabel = {'epsc' 'ipsc'};
            xlim([0.5 2.5])
        end
    else
        jtrFig.(fldLabels{ii})=figure(noGroups*3+7+ii);
        jtrFig.(fldLabels{ii}).Position = [1450+(ii-1)*185 545 165 250];
        hold on
        %plot connecting lines for paired analysis display
        for jj = 1:length(jitter.(fldLabels{ii}).epsc)
            plot([1 2],[jitter.(fldLabels{ii}).epsc(jj),jitter.(fldLabels{ii}).ipsc(jj)],'color',lightgray,'linewidth',1)
            if hJtr.(fldLabels{ii}) == 0 %if data are normal, plot mean
                plot([1 2],[mean(jitter.(fldLabels{ii}).epsc),mean(jitter.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
            elseif hJtr.(fldLabels{ii}) ~=0 %if data non-normal, plot median
                plot([1 2],[median(jitter.(fldLabels{ii}).epsc),median(jitter.(fldLabels{ii}).ipsc)],'color',darkgray,'linewidth',2)
            end
        end
        scatter(ones(1,length(jitter.(fldLabels{ii}).epsc)),jitter.(fldLabels{ii}).epsc,40,lightgray,'filled')
        scatter(2.*ones(1,length(jitter.(fldLabels{ii}).ipsc)),jitter.(fldLabels{ii}).ipsc,40,lightgray,'filled')
        if hJtr.(fldLabels{ii}) == 0 %if data are normal, plot mean w sem
            errorbar(1,mean(jitter.(fldLabels{ii}).epsc),sem(jitter.(fldLabels{ii}).epsc,2),'color',darkgray,'linewidth',3,'CapSize',0)
            scatter(1,mean(jitter.(fldLabels{ii}).epsc),125,darkgray,'filled')
            errorbar(2,mean(jitter.(fldLabels{ii}).ipsc),sem(jitter.(fldLabels{ii}).ipsc,2),'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(2,mean(jitter.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        elseif hJtr.(fldLabels{ii}) ~= 0 %if data non-normal, plot median w IQR
            [lowerError.jitter.(fldLabels{ii}).epsc,upperError.jitter.(fldLabels{ii}).epsc] = iqrError(jitter.(fldLabels{ii}).epsc,2);
            [lowerError.jitter.(fldLabels{ii}).ipsc,upperError.jitter.(fldLabels{ii}).ipsc] = iqrError(jitter.(fldLabels{ii}).ipsc,2);
            errorbar(1,median(jitter.(fldLabels{ii}).epsc),lowerError.jitter.(fldLabels{ii}).epsc,upperError.jitter.(fldLabels{ii}).epsc,'color',darkgray,'linewidth',3,'CapSize',0)
            scatter(1,median(jitter.(fldLabels{ii}).epsc),125,darkgray,'filled')
            errorbar(2,median(jitter.(fldLabels{ii}).ipsc),lowerError.jitter.(fldLabels{ii}).ipsc,upperError.jitter.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',3,'CapSize',0)
            scatter(2,median(jitter.(fldLabels{ii}).ipsc),125,pColor(ii,:),'filled')
        end
        if ii == 1
            allJtr = [jitter.(fldLabels{ii}).epsc jitter.(fldLabels{ii}).ipsc];
        else
            allJtr = [allJtr jitter.(fldLabels{ii}).epsc jitter.(fldLabels{ii}).ipsc];
        end
        ylabel(['Jitter, ',fldLabels{ii},' group'])
        jitterAx.(fldLabels{ii})=gca;
        setAx(jitterAx.(fldLabels{ii}));
        jitterAx.(fldLabels{ii}).XTick = [1 2];
        jitterAx.(fldLabels{ii}).XTickLabel = {'epsc' 'ipsc'};
        xlim([0.5 2.5])
    end
end
for ii = 1:noGroups
    if strcmp(inputs.LatencyTestDropDown.Value,'Ctrl only')
        if strcmp(fldLabels{ii},'ctrl')
            thisFig = jtrFig.(fldLabels{ii}).Number;
            figure(thisFig);
            ylim([0 ceil(1.1*max(allJtr))])
        end
    else
        thisFig = jtrFig.(fldLabels{ii}).Number;
        figure(thisFig);
        ylim([0 ceil(1.1*max(allJtr))])
    end
end

%example traces, ipsc
for ii = 1:noGroups
    t.ipsc=(1000/inputs.SamplerateHzEditField.Value).*(1:1:length(exTrace.(fldLabels{ii}).ipsc));
    ipscTraceFig.(fldLabels{ii}) = figure(ii);
    ipscTraceFig.(fldLabels{ii}).Position = [480+315*(ii-1) 300 300 125];
    hold on
    plot(t.ipsc,exTrace.(fldLabels{ii}).ipsc,'color',pColor(ii,:),'linewidth',2)
    traceAx.(fldLabels{ii}).ipsc = gca;
    setAx(traceAx.(fldLabels{ii}).ipsc);
    ylabel('current (pA)')
    xlabel('t (ms)')
    title([fldLabels{ii} ' example feedforward ipsc response'])
    xlim([0 300])
    if ii == 1
        yLimHolder = traceAx.(fldLabels{ii}).ipsc.YLim;
    else
        yLimHolder = [yLimHolder traceAx.(fldLabels{ii}).ipsc.YLim];
    end
    clear t.ipsc
end
yLimMin = min(yLimHolder);
yLimMax = max(yLimHolder);
for ii = 1:noGroups
    traceAx.(fldLabels{ii}).ipsc.YLim = [yLimMin yLimMax];
end

%example traces, epsc
for ii = 1:noGroups
    t.epsc=(1000/inputs.SamplerateHzEditField.Value).*(1:1:length(exTrace.(fldLabels{ii}).epsc));
    epscTraceFig.(fldLabels{ii}) = figure(ii+noGroups);
    epscTraceFig.(fldLabels{ii}).Position = [480+315*(ii-1) 125 300 125];
    hold on
    plot(t.epsc,exTrace.(fldLabels{ii}).epsc,'color',pColor(ii,:),'linewidth',2)
    traceAx.(fldLabels{ii}).epsc = gca;
    setAx(traceAx.(fldLabels{ii}).epsc);
    ylabel('current (pA)')
    xlabel('t (ms)')
    title([fldLabels{ii} ' example epsc response'])
    xlim([0 300])
    if ii == 1
        yLimHolder = traceAx.(fldLabels{ii}).epsc.YLim;
    else
        yLimHolder = [yLimHolder traceAx.(fldLabels{ii}).epsc.YLim];
    end
    clear t.epsc
end
yLimMin = min(yLimHolder);
yLimMax = max(yLimHolder);
for ii = 1:noGroups
    traceAx.(fldLabels{ii}).epsc.YLim = [yLimMin yLimMax];
end

%% DISPLAY p VALS
output.p
end