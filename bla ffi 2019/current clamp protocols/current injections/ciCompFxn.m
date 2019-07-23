%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% current injections Grp Comp %%%%%%%%%
%%%%%%%%%%% Created: 02-16-2019 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 03-08-2019 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outputdata] = ciCompFxn(inputs)
close all

%% INIT VARS
if strcmp(inputs.FolderNameEditField.Value,'tbd') ~= 1
    fldrName = inputs.FolderNameEditField.Value;
end
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

samplerate = inputs.SamplerateHzEditField.Value;
lightgray=[.75 .75 .75]; %light gray for individual data points

%% LOAD DATA
% Directory
cd(inputs.ciDir);

for ii = 1:noGroups
    %get files
    if inputs.SquareinjectionprotocolCheckBox.Value == 1
        if strcmp(inputs.FolderNameEditField.Value,'tbd') ~= 1
            cd([inputs.ciDir,'/',grpLabels{ii},'/',fldrName,'/current injections (final)']);
        else
            cd([inputs.ciDir,'/current injections (final)/',grpLabels{ii}]);
        end
        contents.(fldLabels{ii}) = dir('*.mat');
        filenames.square.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
        ciFiles.(fldLabels{ii}) = fullfile(cd,filenames.square.(fldLabels{ii}));
        
        for jj = 1:length(ciFiles.(fldLabels{ii}))
            load(ciFiles.(fldLabels{ii}){jj}, 'app')
            
            %data
            outputdata.rm.(fldLabels{ii})(jj) = app.ciData.rm;
            outputdata.mtau.(fldLabels{ii})(jj) = app.ciData.mtau;
            outputdata.vsag.(fldLabels{ii})(jj) = app.ciData.vsag;
            outputdata.nRbnd.(fldLabels{ii})(jj) = app.ciData.nRbnd;
            outputdata.maxFR.(fldLabels{ii})(jj) = app.ciData.maxFR;
            outputdata.fiCurve.(fldLabels{ii}).firingRate(jj,:) = app.ciData.FIcurve(:,1);
            outputdata.spkLat.(fldLabels{ii})(jj) = app.ciData.spkLat;
            outputdata.threshold.(fldLabels{ii})(jj) = mean(app.ciData.rheoThreshold);
            outputdata.spkAmp.(fldLabels{ii})(jj) = mean(app.ciData.valRheoPeaks'-app.ciData.rheoThreshold);
            outputdata.ahpLat.(fldLabels{ii})(jj) = mean(app.ciData.rheoAHPlat);
            outputdata.ahpVal.(fldLabels{ii})(jj) = mean(app.ciData.rheoAHPval);
            outputdata.hw.(fldLabels{ii})(jj) = mean(app.ciData.halfwidth);
            outputdata.frRatio.(fldLabels{ii})(jj) = app.ciData.frRatio;
            outputdata.ampRatio.(fldLabels{ii})(jj) = app.ciData.ampRatio;
            outputdata.broadeningRatio.(fldLabels{ii})(jj) = app.ciData.broadeningRatio;
            outputdata.dAHP.(fldLabels{ii})(jj) = app.ciData.dAHP;
            
            %example traces
            if ii == 1
                exCell.(fldLabels{ii}) = inputs.ExCellEditField.Value;
            elseif ii == 2
                exCell.(fldLabels{ii}) = inputs.ExCellEditField_2.Value;
            elseif ii == 3
                exCell.(fldLabels{ii}) = inputs.ExCellEditField_3.Value;
            end
            if jj == exCell.(fldLabels{ii})
                outputdata.rWaves.(fldLabels{ii}) = app.ciData.rWaves;
                outputdata.phasePlot.(fldLabels{ii}) = app.ciData.dVdt;
                outputdata.rheo.(fldLabels{ii}) = app.ciData.rheo;
                outputdata.fiCurve.(fldLabels{ii}).currentInjection = app.ciData.FIcurve(:,2);
                
                theseFields = fields(app.ciData.ccFile);
                for ll = 1:length(fields(app.ciData.ccFile))
                    outputdata.(theseFields{ll}).(fldLabels{ii}) = app.ciData.ccFile.(theseFields{ll});
                end
                
                outputdata.rheoSweep.(fldLabels{ii}) = size(outputdata.(theseFields{end}).(fldLabels{ii}),3)...
                    - (length(outputdata.fiCurve.(fldLabels{ii}).currentInjection)...
                    - find(outputdata.rheo.(fldLabels{ii})==outputdata.fiCurve.(fldLabels{ii}).currentInjection));
            end
            clear app
        end
    end
    if inputs.RampinjectionprotocolCheckBox.Value == 1
        cd([inputs.ciDir,'/',grpLabels{ii},'/',fldrName,'/rheobase ramp']);
        contents.(fldLabels{ii}) = dir('*.mat');
        filenames.ramp.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
        rampFiles.(fldLabels{ii}) = fullfile(cd,filenames.ramp.(fldLabels{ii}));
        
        if inputs.RampinjectionprotocolCheckBox.Value == 1
            outputdata.threshold.(fldLabels{ii}) = [];
        end
        
        for jj = 1:length(rampFiles.(fldLabels{ii}))
            load(rampFiles.(fldLabels{ii}){jj}, 'app')
            outputdata.vRest.(fldLabels{ii})(jj) = mean(app.rampData.vRest);
            outputdata.threshold.(fldLabels{ii})(jj) = mean(app.rampData.spkThreshold);
            outputdata.rheobase.(fldLabels{ii})(jj) = mean(app.rampData.rheobase);
            
            if ii == 1
                exCell.(fldLabels{ii}) = inputs.ExCellEditField_4.Value;
            elseif ii == 2
                exCell.(fldLabels{ii}) = inputs.ExCellEditField_5.Value;
            elseif ii == 3
                exCell.(fldLabels{ii}) = inputs.ExCellEditField_6.Value;
            end
            if jj == exCell.(fldLabels{ii})
                outputdata.rampVoltage.(fldLabels{ii})(:,:) = app.rampData.rampData(:,1,:);
                outputdata.rampCurrent.(fldLabels{ii})(:,:) = app.rampData.rampData(:,2,:);
            end
        end
    end
end

%time vector
if inputs.SquareinjectionprotocolCheckBox.Value == 1
    t.square = samplerate^-1:samplerate^-1:length(outputdata.(theseFields{1}).(fldLabels{ii}))/samplerate;
end
if inputs.RampinjectionprotocolCheckBox.Value == 1
    t.ramp = samplerate^-1:samplerate^-1:length(outputdata.rampVoltage.(fldLabels{ii}))/samplerate;
end

%% ANALYSIS
%threshold, normality
if noGroups == 3
    for ii = 1:noGroups
        res.threshold.(fldLabels{ii}) = mean(outputdata.threshold.(fldLabels{ii})) - outputdata.threshold.(fldLabels{ii});
    end
    allRes.threshold = [res.threshold.(fldLabels{1})'; res.threshold.(fldLabels{2})'; res.threshold.(fldLabels{3})'];
    [hThreshold,~]=adtest(allRes.threshold);
elseif noGroups < 3
    for ii = 1:noGroups
        [holdH.(fldLabels{ii}),~] = adtest(outputdata.threshold.(fldLabels{ii}));
    end
    if noGroups == 1
        hThreshold = holdH.(fldLabels{1});
    elseif noGroups == 2
        hThreshold = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
    end
    clear holdH
end

%threshold, comparison tests
if noGroups == 3
    if hThreshold == 0 %if residuals normal, run anova
        [outputdata.p.thresholdAnova,~,thresholdstats]=anova1([outputdata.threshold.(fldLabels{1}) outputdata.threshold.(fldLabels{2}) outputdata.threshold.(fldLabels{3})],...
            [ones(1,length(outputdata.threshold.(fldLabels{1}))) 2.*ones(1,length(outputdata.threshold.(fldLabels{2}))) 3.*ones(1,length(outputdata.threshold.(fldLabels{3})))]);
        if outputdata.p.thresholdAnova <= 0.05
            outputdata.thresholdPostHocTukey=multcompare(thresholdstats); %done with tukey's HSD method
        end
    else %run Kruskal wallis
        [outputdata.p.thresholdKW,~,~]=kruskalwallis([outputdata.threshold.(fldLabels{1}) outputdata.threshold.(fldLabels{2}) outputdata.threshold.(fldLabels{3})],...
            [ones(1,length(outputdata.threshold.(fldLabels{1}))) 2.*ones(1,length(outputdata.threshold.(fldLabels{2}))) 3.*ones(1,length(outputdata.threshold.(fldLabels{3})))]);
        if outputdata.p.thresholdKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
            [outputdata.p.thresholdPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.threshold.(fldLabels{1}),outputdata.threshold.(fldLabels{2}));
            [outputdata.p.thresholdPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.threshold.(fldLabels{1}),outputdata.threshold.(fldLabels{3}));
            [outputdata.p.thresholdPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.threshold.(fldLabels{2}),outputdata.threshold.(fldLabels{3}));
            outputdata.pFDR.threshold=drsFDRpval([outputdata.p.thresholdPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.thresholdPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.thresholdPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
        end
    end
elseif noGroups < 3
    if hThreshold == 0 %if groups normal, run ttest
        if noGroups == 2
            [~,outputdata.p.thresholdTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.threshold.(fldLabels{1}),outputdata.threshold.(fldLabels{2}));
        end
    elseif hThreshold ~= 0 %if data not normal, run MWU test
        if noGroups == 2
            [outputdata.p.thresholdMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.threshold.(fldLabels{1}),outputdata.threshold.(fldLabels{2}));
        end
    end
end

if inputs.RampinjectionprotocolCheckBox.Value == 1
    %resting membrane voltage, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.vRest.(fldLabels{ii}) = mean(outputdata.vRest.(fldLabels{ii})) - outputdata.vRest.(fldLabels{ii});
        end
        allRes.vRest = [res.vRest.(fldLabels{1})'; res.vRest.(fldLabels{2})'; res.vRest.(fldLabels{3})'];
        [hRest,~]=adtest(allRes.vRest);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.vRest.(fldLabels{ii}));
        end
        if noGroups == 1
            hRest = holdH.(fldLabels{1});
        elseif noGroups == 2
            hRest = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %v-rest, comparison tests
    if noGroups == 3
        if hRest == 0 %if residuals normal, run anova
            [outputdata.p.restAnova,~,reststats]=anova1([outputdata.vRest.(fldLabels{1}) outputdata.vRest.(fldLabels{2}) outputdata.vRest.(fldLabels{3})],...
                [ones(1,length(outputdata.vRest.(fldLabels{1}))) 2.*ones(1,length(outputdata.vRest.(fldLabels{2}))) 3.*ones(1,length(outputdata.vRest.(fldLabels{3})))]);
            if outputdata.p.restAnova <= 0.05
                outputdata.restPostHocTukey=multcompare(reststats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.restKW,~,~]=kruskalwallis([outputdata.vRest.(fldLabels{1}) outputdata.vRest.(fldLabels{2}) outputdata.vRest.(fldLabels{3})],...
                [ones(1,length(outputdata.vRest.(fldLabels{1}))) 2.*ones(1,length(outputdata.vRest.(fldLabels{2}))) 3.*ones(1,length(outputdata.vRest.(fldLabels{3})))]);
            if outputdata.p.restKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.restPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.vRest.(fldLabels{1}),outputdata.vRest.(fldLabels{2}));
                [outputdata.p.restPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.vRest.(fldLabels{1}),outputdata.vRest.(fldLabels{3}));
                [outputdata.p.restPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.vRest.(fldLabels{2}),outputdata.vRest.(fldLabels{3}));
                outputdata.pFDR.vRest=drsFDRpval([outputdata.p.restPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.restPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.restPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hRest == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.restTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.vRest.(fldLabels{1}),outputdata.vRest.(fldLabels{2}));
            end
        elseif hRest ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.restMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.vRest.(fldLabels{1}),outputdata.vRest.(fldLabels{2}));
            end
        end
    end
    
    %rheobase, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.rheobase.(fldLabels{ii}) = mean(outputdata.rheobase.(fldLabels{ii})) - outputdata.rheobase.(fldLabels{ii});
        end
        allRes.rheobase = [res.rheobase.(fldLabels{1})'; res.rheobase.(fldLabels{2})'; res.rheobase.(fldLabels{3})'];
        [hRheobase,~]=adtest(allRes.rheobase);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.rheobase.(fldLabels{ii}));
        end
        if noGroups == 1
            hRheobase = holdH.(fldLabels{1});
        elseif noGroups == 2
            hRheobase = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %rheobase, comparison tests
    if noGroups == 3
        if hRheobase == 0 %if residuals normal, run anova
            [outputdata.p.rheobaseAnova,~,rheobasestats]=anova1([outputdata.rheobase.(fldLabels{1}) outputdata.rheobase.(fldLabels{2}) outputdata.rheobase.(fldLabels{3})],...
                [ones(1,length(outputdata.rheobase.(fldLabels{1}))) 2.*ones(1,length(outputdata.rheobase.(fldLabels{2}))) 3.*ones(1,length(outputdata.rheobase.(fldLabels{3})))]);
            if outputdata.p.rheobaseAnova <= 0.05
                outputdata.rheobasePostHocTukey=multcompare(rheobasestats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.rheobaseKW,~,~]=kruskalwallis([outputdata.rheobase.(fldLabels{1}) outputdata.rheobase.(fldLabels{2}) outputdata.rheobase.(fldLabels{3})],...
                [ones(1,length(outputdata.rheobase.(fldLabels{1}))) 2.*ones(1,length(outputdata.rheobase.(fldLabels{2}))) 3.*ones(1,length(outputdata.rheobase.(fldLabels{3})))]);
            if outputdata.p.restKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.rheobasePostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.rheobase.(fldLabels{1}),outputdata.rheobase.(fldLabels{2}));
                [outputdata.p.rheobasePostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.rheobase.(fldLabels{1}),outputdata.rheobase.(fldLabels{3}));
                [outputdata.p.rheobasePostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.rheobase.(fldLabels{2}),outputdata.rheobase.(fldLabels{3}));
                outputdata.pFDR.rheobase=drsFDRpval([outputdata.p.rheobasePostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.rheobasePostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.rheobasePostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hRheobase == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.rheobaseTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.rheobase.(fldLabels{1}),outputdata.rheobase.(fldLabels{2}));
            end
        elseif hRheobase ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.rheobaseMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.rheobase.(fldLabels{1}),outputdata.rheobase.(fldLabels{2}));
            end
        end
    end
end
if inputs.SquareinjectionprotocolCheckBox.Value == 1
    %membrane resistance, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.rm.(fldLabels{ii}) = mean(outputdata.rm.(fldLabels{ii})) - outputdata.rm.(fldLabels{ii});
        end
        allRes.rm = [res.rm.(fldLabels{1})'; res.rm.(fldLabels{2})'; res.rm.(fldLabels{3})'];
        [hRm,~]=adtest(allRes.rm);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.rm.(fldLabels{ii}));
        end
        if noGroups == 1
            hRm = holdH.(fldLabels{1});
        elseif noGroups == 2
            hRm = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %rm, comparison tests
    if noGroups == 3
        if hRm == 0 %if residuals normal, run anova
            [outputdata.p.rmAnova,~,rmstats]=anova1([outputdata.rm.(fldLabels{1}) outputdata.rm.(fldLabels{2}) outputdata.rm.(fldLabels{3})],...
                [ones(1,length(outputdata.rm.(fldLabels{1}))) 2.*ones(1,length(outputdata.rm.(fldLabels{2}))) 3.*ones(1,length(outputdata.rm.(fldLabels{3})))]);
            if outputdata.p.rmAnova <= 0.05
                outputdata.rmPostHocTukey=multcompare(rmstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.rmKW,~,~]=kruskalwallis([outputdata.rm.(fldLabels{1}) outputdata.rm.(fldLabels{2}) outputdata.rm.(fldLabels{3})],...
                [ones(1,length(outputdata.rm.(fldLabels{1}))) 2.*ones(1,length(outputdata.rm.(fldLabels{2}))) 3.*ones(1,length(outputdata.rm.(fldLabels{3})))]);
            if outputdata.p.rmKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.rmPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.rm.(fldLabels{1}),outputdata.rm.(fldLabels{2}));
                [outputdata.p.rmPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.rm.(fldLabels{1}),outputdata.rm.(fldLabels{3}));
                [outputdata.p.rmPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.rm.(fldLabels{2}),outputdata.rm.(fldLabels{3}));
                outputdata.pFDR.rm=drsFDRpval([outputdata.p.rmPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.rmPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.rmPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hRm == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.rmTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.rm.(fldLabels{1}),outputdata.rm.(fldLabels{2}));
            end
        elseif hRm ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.rmMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.rm.(fldLabels{1}),outputdata.rm.(fldLabels{2}));
            end
        end
    end
    
    %membrane decay constant, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.mtau.(fldLabels{ii}) = mean(outputdata.mtau.(fldLabels{ii})) - outputdata.mtau.(fldLabels{ii});
        end
        allRes.mtau = [res.mtau.(fldLabels{1})'; res.mtau.(fldLabels{2})'; res.mtau.(fldLabels{3})'];
        [hMemTau,~]=adtest(allRes.mtau);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.mtau.(fldLabels{ii}));
        end
        if noGroups == 1
            hMemTau = holdH.(fldLabels{1});
        elseif noGroups == 2
            hMemTau = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %membrane decay constant, comparison tests
    if noGroups == 3
        if hMemTau == 0 %if residuals normal, run anova
            [outputdata.p.mtauAnova,~,mtaustats]=anova1([outputdata.mtau.(fldLabels{1}) outputdata.mtau.(fldLabels{2}) outputdata.mtau.(fldLabels{3})],...
                [ones(1,length(outputdata.mtau.(fldLabels{1}))) 2.*ones(1,length(outputdata.mtau.(fldLabels{2}))) 3.*ones(1,length(outputdata.mtau.(fldLabels{3})))]);
            if outputdata.p.mtauAnova <= 0.05
                outputdata.mtauPostHocTukey=multcompare(mtaustats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.mtauKW,~,~]=kruskalwallis([outputdata.mtau.(fldLabels{1}) outputdata.mtau.(fldLabels{2}) outputdata.mtau.(fldLabels{3})],...
                [ones(1,length(outputdata.mtau.(fldLabels{1}))) 2.*ones(1,length(outputdata.mtau.(fldLabels{2}))) 3.*ones(1,length(outputdata.mtau.(fldLabels{3})))]);
            if outputdata.p.mtauKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.mtauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.mtau.(fldLabels{1}),outputdata.mtau.(fldLabels{2}));
                [outputdata.p.mtauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.mtau.(fldLabels{1}),outputdata.mtau.(fldLabels{3}));
                [outputdata.p.mtauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.mtau.(fldLabels{2}),outputdata.mtau.(fldLabels{3}));
                outputdata.pFDR.mtau=drsFDRpval([outputdata.p.mtauPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.mtauPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.mtauPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hMemTau == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.mtauTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.mtau.(fldLabels{1}),outputdata.mtau.(fldLabels{2}));
            end
        elseif hMemTau ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.mtauMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.mtau.(fldLabels{1}),outputdata.mtau.(fldLabels{2}));
            end
        end
    end
    
    %voltage sag, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.vsag.(fldLabels{ii}) = mean(outputdata.vsag.(fldLabels{ii})) - outputdata.vsag.(fldLabels{ii});
        end
        allRes.vsag = [res.vsag.(fldLabels{1})'; res.vsag.(fldLabels{2})'; res.vsag.(fldLabels{3})'];
        [hSag,~]=adtest(allRes.vsag);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.vsag.(fldLabels{ii}));
        end
        if noGroups == 1
            hSag = holdH.(fldLabels{1});
        elseif noGroups == 2
            hSag = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %voltage sag, comparison tests
    if noGroups == 3
        if hSag == 0 %if residuals normal, run anova
            [outputdata.p.vsagAnova,~,vsagstats]=anova1([outputdata.vsag.(fldLabels{1}) outputdata.vsag.(fldLabels{2}) outputdata.vsag.(fldLabels{3})],...
                [ones(1,length(outputdata.vsag.(fldLabels{1}))) 2.*ones(1,length(outputdata.vsag.(fldLabels{2}))) 3.*ones(1,length(outputdata.vsag.(fldLabels{3})))]);
            if outputdata.p.vsagAnova <= 0.05
                outputdata.vsagPostHocTukey=multcompare(vsagstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.vsagKW,~,~]=kruskalwallis([outputdata.vsag.(fldLabels{1}) outputdata.vsag.(fldLabels{2}) outputdata.vsag.(fldLabels{3})],...
                [ones(1,length(outputdata.vsag.(fldLabels{1}))) 2.*ones(1,length(outputdata.vsag.(fldLabels{2}))) 3.*ones(1,length(outputdata.vsag.(fldLabels{3})))]);
            if outputdata.p.vsagKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.vsagPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.vsag.(fldLabels{1}),outputdata.vsag.(fldLabels{2}));
                [outputdata.p.vsagPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.vsag.(fldLabels{1}),outputdata.vsag.(fldLabels{3}));
                [outputdata.p.vsagPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.vsag.(fldLabels{2}),outputdata.vsag.(fldLabels{3}));
                outputdata.pFDR.vsag=drsFDRpval([outputdata.p.vsagPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.vsagPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.vsagPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hSag == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.vsagTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.vsag.(fldLabels{1}),outputdata.vsag.(fldLabels{2}));
            end
        elseif hSag ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.vsagMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.vsag.(fldLabels{1}),outputdata.vsag.(fldLabels{2}));
            end
        end
    end
    
    %rebound spikes, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.nRbnd.(fldLabels{ii}) = mean(outputdata.nRbnd.(fldLabels{ii})) - outputdata.nRbnd.(fldLabels{ii});
        end
        allRes.nRbnd = [res.nRbnd.(fldLabels{1})'; res.nRbnd.(fldLabels{2})'; res.nRbnd.(fldLabels{3})'];
        [hRebound,~]=adtest(allRes.nRbnd);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.nRbnd.(fldLabels{ii}));
        end
        if noGroups == 1
            hRebound = holdH.(fldLabels{1});
        elseif noGroups == 2
            hRebound = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %rebound spikes, comparison tests
    if noGroups == 3
        if hRebound == 0 %if residuals normal, run anova
            [outputdata.p.nRbndAnova,~,nRbndstats]=anova1([outputdata.nRbnd.(fldLabels{1}) outputdata.nRbnd.(fldLabels{2}) outputdata.nRbnd.(fldLabels{3})],...
                [ones(1,length(outputdata.nRbnd.(fldLabels{1}))) 2.*ones(1,length(outputdata.nRbnd.(fldLabels{2}))) 3.*ones(1,length(outputdata.nRbnd.(fldLabels{3})))]);
            if outputdata.p.nRbndAnova <= 0.05
                outputdata.nRbndPostHocTukey=multcompare(nRbndstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.nRbndKW,~,~]=kruskalwallis([outputdata.nRbnd.(fldLabels{1}) outputdata.nRbnd.(fldLabels{2}) outputdata.nRbnd.(fldLabels{3})],...
                [ones(1,length(outputdata.nRbnd.(fldLabels{1}))) 2.*ones(1,length(outputdata.nRbnd.(fldLabels{2}))) 3.*ones(1,length(outputdata.nRbnd.(fldLabels{3})))]);
            if outputdata.p.nRbndKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.nRbndPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.nRbnd.(fldLabels{1}),outputdata.nRbnd.(fldLabels{2}));
                [outputdata.p.nRbndPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.nRbnd.(fldLabels{1}),outputdata.nRbnd.(fldLabels{3}));
                [outputdata.p.nRbndPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.nRbnd.(fldLabels{2}),outputdata.nRbnd.(fldLabels{3}));
                outputdata.pFDR.nRbnd=drsFDRpval([outputdata.p.nRbndPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.nRbndPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.nRbndPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hRebound == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.nRbndTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.nRbnd.(fldLabels{1}),outputdata.nRbnd.(fldLabels{2}));
            end
        elseif hRebound ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.nRbndMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.nRbnd.(fldLabels{1}),outputdata.nRbnd.(fldLabels{2}));
            end
        end
    end
    
    %max firing rate, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.maxFR.(fldLabels{ii}) = mean(outputdata.maxFR.(fldLabels{ii})) - outputdata.maxFR.(fldLabels{ii});
        end
        allRes.maxFR = [res.maxFR.(fldLabels{1})'; res.maxFR.(fldLabels{2})'; res.maxFR.(fldLabels{3})'];
        [hMaxFR,~]=adtest(allRes.maxFR);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.maxFR.(fldLabels{ii}));
        end
        if noGroups == 1
            hMaxFR = holdH.(fldLabels{1});
        elseif noGroups == 2
            hMaxFR = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %max firing rate, comparison tests
    if noGroups == 3
        if hMaxFR == 0 %if residuals normal, run anova
            [outputdata.p.maxFRAnova,~,maxFRstats]=anova1([outputdata.maxFR.(fldLabels{1}) outputdata.maxFR.(fldLabels{2}) outputdata.maxFR.(fldLabels{3})],...
                [ones(1,length(outputdata.maxFR.(fldLabels{1}))) 2.*ones(1,length(outputdata.maxFR.(fldLabels{2}))) 3.*ones(1,length(outputdata.maxFR.(fldLabels{3})))]);
            if outputdata.p.maxFRAnova <= 0.05
                outputdata.maxFRPostHocTukey=multcompare(maxFRstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.maxFRKW,~,~]=kruskalwallis([outputdata.maxFR.(fldLabels{1}) outputdata.maxFR.(fldLabels{2}) outputdata.maxFR.(fldLabels{3})],...
                [ones(1,length(outputdata.maxFR.(fldLabels{1}))) 2.*ones(1,length(outputdata.maxFR.(fldLabels{2}))) 3.*ones(1,length(outputdata.maxFR.(fldLabels{3})))]);
            if outputdata.p.maxFRKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.maxFRPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.maxFR.(fldLabels{1}),outputdata.maxFR.(fldLabels{2}));
                [outputdata.p.maxFRPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.maxFR.(fldLabels{1}),outputdata.maxFR.(fldLabels{3}));
                [outputdata.p.maxFRPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.maxFR.(fldLabels{2}),outputdata.maxFR.(fldLabels{3}));
                outputdata.pFDR.maxFR=drsFDRpval([outputdata.p.maxFRPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.maxFRPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.maxFRPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hMaxFR == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.maxFRTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.maxFR.(fldLabels{1}),outputdata.maxFR.(fldLabels{2}));
            end
        elseif hMaxFR ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.maxFRMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.maxFR.(fldLabels{1}),outputdata.rheobase.(fldLabels{2}));
            end
        end
    end
    
    %spike latency, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.spkLat.(fldLabels{ii}) = mean(outputdata.spkLat.(fldLabels{ii})) - outputdata.spkLat.(fldLabels{ii});
        end
        allRes.spkLat = [res.spkLat.(fldLabels{1})'; res.spkLat.(fldLabels{2})'; res.spkLat.(fldLabels{3})'];
        [hSpkLat,~]=adtest(allRes.spkLat);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.spkLat.(fldLabels{ii}));
        end
        if noGroups == 1
            hSpkLat = holdH.(fldLabels{1});
        elseif noGroups == 2
            hSpkLat = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %spike latency, comparison tests
    if noGroups == 3
        if hSpkLat == 0 %if residuals normal, run anova
            [outputdata.p.spkLatAnova,~,spkLatstats]=anova1([outputdata.spkLat.(fldLabels{1}) outputdata.spkLat.(fldLabels{2}) outputdata.spkLat.(fldLabels{3})],...
                [ones(1,length(outputdata.spkLat.(fldLabels{1}))) 2.*ones(1,length(outputdata.spkLat.(fldLabels{2}))) 3.*ones(1,length(outputdata.spkLat.(fldLabels{3})))]);
            if outputdata.p.spkLatAnova <= 0.05
                outputdata.spkLatPostHocTukey=multcompare(spkLatstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.spkLatKW,~,~]=kruskalwallis([outputdata.spkLat.(fldLabels{1}) outputdata.spkLat.(fldLabels{2}) outputdata.spkLat.(fldLabels{3})],...
                [ones(1,length(outputdata.spkLat.(fldLabels{1}))) 2.*ones(1,length(outputdata.spkLat.(fldLabels{2}))) 3.*ones(1,length(outputdata.spkLat.(fldLabels{3})))]);
            if outputdata.p.spkLatKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.spkLatPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.spkLat.(fldLabels{1}),outputdata.spkLat.(fldLabels{2}));
                [outputdata.p.spkLatPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.spkLat.(fldLabels{1}),outputdata.spkLat.(fldLabels{3}));
                [outputdata.p.spkLatPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.spkLat.(fldLabels{2}),outputdata.spkLat.(fldLabels{3}));
                outputdata.pFDR.spkLat=drsFDRpval([outputdata.p.spkLatPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.spkLatPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.spkLatPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hSpkLat == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.spkLatTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.spkLat.(fldLabels{1}),outputdata.spkLat.(fldLabels{2}));
            end
        elseif hSpkLat ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.spkLatMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.spkLat.(fldLabels{1}),outputdata.spkLat.(fldLabels{2}));
            end
        end
    end
    
    %spike amplitude, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.spkAmp.(fldLabels{ii}) = mean(outputdata.spkAmp.(fldLabels{ii})) - outputdata.spkAmp.(fldLabels{ii});
        end
        allRes.spkAmp = [res.spkAmp.(fldLabels{1})'; res.spkAmp.(fldLabels{2})'; res.spkAmp.(fldLabels{3})'];
        [hSpkAmp,~]=adtest(allRes.spkAmp);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.spkAmp.(fldLabels{ii}));
        end
        if noGroups == 1
            hSpkAmp = holdH.(fldLabels{1});
        elseif noGroups == 2
            hSpkAmp = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %spike amplitude, comparison tests
    if noGroups == 3
        if hSpkAmp == 0 %if residuals normal, run anova
            [outputdata.p.spkAmpAnova,~,spkAmpstats]=anova1([outputdata.spkAmp.(fldLabels{1}) outputdata.spkAmp.(fldLabels{2}) outputdata.spkAmp.(fldLabels{3})],...
                [ones(1,length(outputdata.spkAmp.(fldLabels{1}))) 2.*ones(1,length(outputdata.spkAmp.(fldLabels{2}))) 3.*ones(1,length(outputdata.spkAmp.(fldLabels{3})))]);
            if outputdata.p.spkAmpAnova <= 0.05
                outputdata.spkAmpPostHocTukey=multcompare(spkAmpstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.spkAmpKW,~,~]=kruskalwallis([outputdata.spkAmp.(fldLabels{1}) outputdata.spkAmp.(fldLabels{2}) outputdata.spkAmp.(fldLabels{3})],...
                [ones(1,length(outputdata.spkAmp.(fldLabels{1}))) 2.*ones(1,length(outputdata.spkAmp.(fldLabels{2}))) 3.*ones(1,length(outputdata.spkAmp.(fldLabels{3})))]);
            if outputdata.p.spkAmpKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.spkAmpPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.spkAmp.(fldLabels{1}),outputdata.spkAmp.(fldLabels{2}));
                [outputdata.p.spkAmpPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.spkAmp.(fldLabels{1}),outputdata.spkAmp.(fldLabels{3}));
                [outputdata.p.spkAmpPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.spkAmp.(fldLabels{2}),outputdata.spkAmp.(fldLabels{3}));
                outputdata.pFDR.spkAmp=drsFDRpval([outputdata.p.spkAmpPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.spkAmpPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.spkAmpPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hSpkAmp == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.spkAmpTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.spkAmp.(fldLabels{1}),outputdata.spkAmp.(fldLabels{2}));
            end
        elseif hSpkAmp ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.spkAmpMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.spkAmp.(fldLabels{1}),outputdata.spkAmp.(fldLabels{2}));
            end
        end
    end
    
    %AHP latency, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.ahpLat.(fldLabels{ii}) = mean(outputdata.ahpLat.(fldLabels{ii})) - outputdata.ahpLat.(fldLabels{ii});
        end
        allRes.ahpLat = [res.ahpLat.(fldLabels{1})'; res.ahpLat.(fldLabels{2})'; res.ahpLat.(fldLabels{3})'];
        [hLatAHP,~]=adtest(allRes.ahpLat);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.ahpLat.(fldLabels{ii}));
        end
        if noGroups == 1
            hLatAHP = holdH.(fldLabels{1});
        elseif noGroups == 2
            hLatAHP = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %AHP latency, comparison tests
    if noGroups == 3
        if hLatAHP == 0 %if residuals normal, run anova
            [outputdata.p.ahpLatAnova,~,ahpLatstats]=anova1([outputdata.ahpLat.(fldLabels{1}) outputdata.ahpLat.(fldLabels{2}) outputdata.ahpLat.(fldLabels{3})],...
                [ones(1,length(outputdata.ahpLat.(fldLabels{1}))) 2.*ones(1,length(outputdata.ahpLat.(fldLabels{2}))) 3.*ones(1,length(outputdata.ahpLat.(fldLabels{3})))]);
            if outputdata.p.ahpLatAnova <= 0.05
                outputdata.ahpLatPostHocTukey=multcompare(ahpLatstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.ahpLatKW,~,~]=kruskalwallis([outputdata.ahpLat.(fldLabels{1}) outputdata.ahpLat.(fldLabels{2}) outputdata.ahpLat.(fldLabels{3})],...
                [ones(1,length(outputdata.ahpLat.(fldLabels{1}))) 2.*ones(1,length(outputdata.ahpLat.(fldLabels{2}))) 3.*ones(1,length(outputdata.ahpLat.(fldLabels{3})))]);
            if outputdata.p.ahpLatKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.ahpLatPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.ahpLat.(fldLabels{1}),outputdata.ahpLat.(fldLabels{2}));
                [outputdata.p.ahpLatPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.ahpLat.(fldLabels{1}),outputdata.ahpLat.(fldLabels{3}));
                [outputdata.p.ahpLatPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.ahpLat.(fldLabels{2}),outputdata.ahpLat.(fldLabels{3}));
                outputdata.pFDR.ahpLat=drsFDRpval([outputdata.p.ahpLatPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.ahpLatPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.ahpLatPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hLatAHP == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.ahpLatTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.ahpLat.(fldLabels{1}),outputdata.ahpLat.(fldLabels{2}));
            end
        elseif hLatAHP ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.ahpLatMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.ahpLat.(fldLabels{1}),outputdata.ahpLat.(fldLabels{2}));
            end
        end
    end
    
    %AHP value, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.ahpVal.(fldLabels{ii}) = mean(outputdata.ahpVal.(fldLabels{ii})) - outputdata.ahpVal.(fldLabels{ii});
        end
        allRes.ahpVal = [res.ahpVal.(fldLabels{1})'; res.ahpVal.(fldLabels{2})'; res.ahpVal.(fldLabels{3})'];
        [hValAHP,~]=adtest(allRes.ahpVal);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.ahpVal.(fldLabels{ii}));
        end
        if noGroups == 1
            hValAHP = holdH.(fldLabels{1});
        elseif noGroups == 2
            hValAHP = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %AHP value, comparison tests
    if noGroups == 3
        if hValAHP == 0 %if residuals normal, run anova
            [outputdata.p.ahpValAnova,~,ahpValstats]=anova1([outputdata.ahpVal.(fldLabels{1}) outputdata.ahpVal.(fldLabels{2}) outputdata.ahpVal.(fldLabels{3})],...
                [ones(1,length(outputdata.ahpVal.(fldLabels{1}))) 2.*ones(1,length(outputdata.ahpVal.(fldLabels{2}))) 3.*ones(1,length(outputdata.ahpVal.(fldLabels{3})))]);
            if outputdata.p.ahpValAnova <= 0.05
                outputdata.ahpValPostHocTukey=multcompare(ahpValstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.ahpValKW,~,~]=kruskalwallis([outputdata.ahpVal.(fldLabels{1}) outputdata.ahpVal.(fldLabels{2}) outputdata.ahpVal.(fldLabels{3})],...
                [ones(1,length(outputdata.ahpVal.(fldLabels{1}))) 2.*ones(1,length(outputdata.ahpVal.(fldLabels{2}))) 3.*ones(1,length(outputdata.ahpVal.(fldLabels{3})))]);
            if outputdata.p.ahpValKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.ahpValPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.ahpVal.(fldLabels{1}),outputdata.ahpVal.(fldLabels{2}));
                [outputdata.p.ahpValPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.ahpVal.(fldLabels{1}),outputdata.ahpVal.(fldLabels{3}));
                [outputdata.p.ahpValPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.ahpVal.(fldLabels{2}),outputdata.ahpVal.(fldLabels{3}));
                outputdata.pFDR.ahpVal=drsFDRpval([outputdata.p.ahpValPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.ahpValPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.ahpValPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hValAHP == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.ahpValTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.ahpVal.(fldLabels{1}),outputdata.ahpVal.(fldLabels{2}));
            end
        elseif hValAHP ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.ahpValMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.ahpVal.(fldLabels{1}),outputdata.ahpVal.(fldLabels{2}));
            end
        end
    end
    
    %halfwidth, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.hw.(fldLabels{ii}) = mean(outputdata.hw.(fldLabels{ii})) - outputdata.hw.(fldLabels{ii});
        end
        allRes.hw = [res.hw.(fldLabels{1})'; res.hw.(fldLabels{2})'; res.hw.(fldLabels{3})'];
        [hHalfwidth,~]=adtest(allRes.hw);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.hw.(fldLabels{ii}));
        end
        if noGroups == 1
            hHalfwidth = holdH.(fldLabels{1});
        elseif noGroups == 2
            hHalfwidth = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %halfwidth, comparison tests
    if noGroups == 3
        if hHalfwidth == 0 %if residuals normal, run anova
            [outputdata.p.hwAnova,~,hwstats]=anova1([outputdata.hw.(fldLabels{1}) outputdata.hw.(fldLabels{2}) outputdata.hw.(fldLabels{3})],...
                [ones(1,length(outputdata.hw.(fldLabels{1}))) 2.*ones(1,length(outputdata.hw.(fldLabels{2}))) 3.*ones(1,length(outputdata.hw.(fldLabels{3})))]);
            if outputdata.p.hwAnova <= 0.05
                outputdata.hwPostHocTukey=multcompare(hwstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.hwKW,~,~]=kruskalwallis([outputdata.hw.(fldLabels{1}) outputdata.hw.(fldLabels{2}) outputdata.hw.(fldLabels{3})],...
                [ones(1,length(outputdata.hw.(fldLabels{1}))) 2.*ones(1,length(outputdata.hw.(fldLabels{2}))) 3.*ones(1,length(outputdata.hw.(fldLabels{3})))]);
            if outputdata.p.hwKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.hwPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.hw.(fldLabels{1}),outputdata.hw.(fldLabels{2}));
                [outputdata.p.hwPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.hw.(fldLabels{1}),outputdata.hw.(fldLabels{3}));
                [outputdata.p.hwPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.hw.(fldLabels{2}),outputdata.hw.(fldLabels{3}));
                outputdata.pFDR.hw=drsFDRpval([outputdata.p.hwPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.hwPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.hwPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hHalfwidth == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.hwTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.hw.(fldLabels{1}),outputdata.hw.(fldLabels{2}));
            end
        elseif hHalfwidth ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.hwMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.hw.(fldLabels{1}),outputdata.hw.(fldLabels{2}));
            end
        end
    end
    
    %Spike Frequency Accomodation, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.frRatio.(fldLabels{ii}) = mean(outputdata.frRatio.(fldLabels{ii})) - outputdata.frRatio.(fldLabels{ii});
        end
        allRes.frRatio = [res.frRatio.(fldLabels{1})'; res.frRatio.(fldLabels{2})'; res.frRatio.(fldLabels{3})'];
        [hFRA,~]=adtest(allRes.frRatio);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.frRatio.(fldLabels{ii}));
        end
        if noGroups == 1
            hFRA = holdH.(fldLabels{1});
        elseif noGroups == 2
            hFRA = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %Spike Frequency Accomodation, comparison tests
    if noGroups == 3
        if hFRA == 0 %if residuals normal, run anova
            [outputdata.p.frRatioAnova,~,frRatiostats]=anova1([outputdata.frRatio.(fldLabels{1}) outputdata.frRatio.(fldLabels{2}) outputdata.frRatio.(fldLabels{3})],...
                [ones(1,length(outputdata.frRatio.(fldLabels{1}))) 2.*ones(1,length(outputdata.frRatio.(fldLabels{2}))) 3.*ones(1,length(outputdata.frRatio.(fldLabels{3})))]);
            if outputdata.p.frRatioAnova <= 0.05
                outputdata.frRatioPostHocTukey=multcompare(frRatiostats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.frRatioKW,~,~]=kruskalwallis([outputdata.frRatio.(fldLabels{1}) outputdata.frRatio.(fldLabels{2}) outputdata.frRatio.(fldLabels{3})],...
                [ones(1,length(outputdata.frRatio.(fldLabels{1}))) 2.*ones(1,length(outputdata.frRatio.(fldLabels{2}))) 3.*ones(1,length(outputdata.frRatio.(fldLabels{3})))]);
            if outputdata.p.frRatioKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.frRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.frRatio.(fldLabels{1}),outputdata.frRatio.(fldLabels{2}));
                [outputdata.p.frRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.frRatio.(fldLabels{1}),outputdata.frRatio.(fldLabels{3}));
                [outputdata.p.frRatioPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.frRatio.(fldLabels{2}),outputdata.frRatio.(fldLabels{3}));
                outputdata.pFDR.frRatio=drsFDRpval([outputdata.p.frRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.frRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.frRatioPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hFRA == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.frRatioTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.frRatio.(fldLabels{1}),outputdata.frRatio.(fldLabels{2}));
            end
        elseif hFRA ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.frRatioMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.frRatio.(fldLabels{1}),outputdata.frRatio.(fldLabels{2}));
            end
        end
    end
    
    %Spike Amplitude Accomodation, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.ampRatio.(fldLabels{ii}) = mean(outputdata.ampRatio.(fldLabels{ii})) - outputdata.ampRatio.(fldLabels{ii});
        end
        allRes.ampRatio = [res.ampRatio.(fldLabels{1})'; res.ampRatio.(fldLabels{2})'; res.ampRatio.(fldLabels{3})'];
        [hSAA,~]=adtest(allRes.ampRatio);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.ampRatio.(fldLabels{ii}));
        end
        if noGroups == 1
            hSAA = holdH.(fldLabels{1});
        elseif noGroups == 2
            hSAA = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %Spike Amplitude Accomodation, comparison tests
    if noGroups == 3
        if hSAA == 0 %if residuals normal, run anova
            [outputdata.p.ampRatioAnova,~,ampRatiostats]=anova1([outputdata.ampRatio.(fldLabels{1}) outputdata.ampRatio.(fldLabels{2}) outputdata.ampRatio.(fldLabels{3})],...
                [ones(1,length(outputdata.ampRatio.(fldLabels{1}))) 2.*ones(1,length(outputdata.ampRatio.(fldLabels{2}))) 3.*ones(1,length(outputdata.ampRatio.(fldLabels{3})))]);
            if outputdata.p.ampRatioAnova <= 0.05
                outputdata.ampRatioPostHocTukey=multcompare(ampRatiostats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.ampRatioKW,~,~]=kruskalwallis([outputdata.ampRatio.(fldLabels{1}) outputdata.ampRatio.(fldLabels{2}) outputdata.ampRatio.(fldLabels{3})],...
                [ones(1,length(outputdata.ampRatio.(fldLabels{1}))) 2.*ones(1,length(outputdata.ampRatio.(fldLabels{2}))) 3.*ones(1,length(outputdata.ampRatio.(fldLabels{3})))]);
            if outputdata.p.ampRatioKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.ampRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.ampRatio.(fldLabels{1}),outputdata.ampRatio.(fldLabels{2}));
                [outputdata.p.ampRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.ampRatio.(fldLabels{1}),outputdata.ampRatio.(fldLabels{3}));
                [outputdata.p.ampRatioPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.ampRatio.(fldLabels{2}),outputdata.ampRatio.(fldLabels{3}));
                outputdata.pFDR.ampRatio=drsFDRpval([outputdata.p.ampRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.ampRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.ampRatioPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hSAA == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.ampRatioTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.ampRatio.(fldLabels{1}),outputdata.ampRatio.(fldLabels{2}));
            end
        elseif hSAA ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.ampRatioMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.ampRatio.(fldLabels{1}),outputdata.ahpLat.(fldLabels{2}));
            end
        end
    end
    
    %Action Potential Broadening, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.broadeningRatio.(fldLabels{ii}) = mean(outputdata.broadeningRatio.(fldLabels{ii})) - outputdata.broadeningRatio.(fldLabels{ii});
        end
        allRes.broadeningRatio = [res.broadeningRatio.(fldLabels{1})'; res.broadeningRatio.(fldLabels{2})'; res.broadeningRatio.(fldLabels{3})'];
        [hBA,~]=adtest(allRes.broadeningRatio);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.broadeningRatio.(fldLabels{ii}));
        end
        if noGroups == 1
            hBA = holdH.(fldLabels{1});
        elseif noGroups == 2
            hBA = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %Action Potential Broadening, comparison tests
    if noGroups == 3
        if hBA == 0 %if residuals normal, run anova
            [outputdata.p.broadeningRatioAnova,~,broadeningRatiostats]=anova1([outputdata.broadeningRatio.(fldLabels{1}) outputdata.broadeningRatio.(fldLabels{2}) outputdata.broadeningRatio.(fldLabels{3})],...
                [ones(1,length(outputdata.broadeningRatio.(fldLabels{1}))) 2.*ones(1,length(outputdata.broadeningRatio.(fldLabels{2}))) 3.*ones(1,length(outputdata.broadeningRatio.(fldLabels{3})))]);
            if outputdata.p.broadeningRatioAnova <= 0.05
                outputdata.broadeningRatioPostHocTukey=multcompare(broadeningRatiostats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.broadeningRatioKW,~,~]=kruskalwallis([outputdata.broadeningRatio.(fldLabels{1}) outputdata.broadeningRatio.(fldLabels{2}) outputdata.broadeningRatio.(fldLabels{3})],...
                [ones(1,length(outputdata.broadeningRatio.(fldLabels{1}))) 2.*ones(1,length(outputdata.broadeningRatio.(fldLabels{2}))) 3.*ones(1,length(outputdata.broadeningRatio.(fldLabels{3})))]);
            if outputdata.p.broadeningRatioKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.broadeningRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.broadeningRatio.(fldLabels{1}),outputdata.broadeningRatio.(fldLabels{2}));
                [outputdata.p.broadeningRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.broadeningRatio.(fldLabels{1}),outputdata.broadeningRatio.(fldLabels{3}));
                [outputdata.p.broadeningRatioPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.broadeningRatio.(fldLabels{2}),outputdata.broadeningRatio.(fldLabels{3}));
                outputdata.pFDR.broadeningRatio=drsFDRpval([outputdata.p.broadeningRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.broadeningRatioPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.broadeningRatioPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hBA == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.broadeningRatioTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.broadeningRatio.(fldLabels{1}),outputdata.broadeningRatio.(fldLabels{2}));
            end
        elseif hBA ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.broadeningRatioMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.broadeningRatio.(fldLabels{1}),outputdata.ahpLat.(fldLabels{2}));
            end
        end
    end
    
    %delta AHP, normality
    if noGroups == 3
        for ii = 1:noGroups
            res.dAHP.(fldLabels{ii}) = mean(outputdata.dAHP.(fldLabels{ii})) - outputdata.dAHP.(fldLabels{ii});
        end
        allRes.dAHP = [res.dAHP.(fldLabels{1})'; res.dAHP.(fldLabels{2})'; res.dAHP.(fldLabels{3})'];
        [hDAHP,~]=adtest(allRes.dAHP);
    elseif noGroups < 3
        for ii = 1:noGroups
            [holdH.(fldLabels{ii}),~] = adtest(outputdata.dAHP.(fldLabels{ii}));
        end
        if noGroups == 1
            hDAHP = holdH.(fldLabels{1});
        elseif noGroups == 2
            hDAHP = sum([holdH.(fldLabels{1}) holdH.(fldLabels{2})]);
        end
        clear holdH
    end
    
    %delta AHP, comparison tests
    if noGroups == 3
        if hDAHP == 0 %if residuals normal, run anova
            [outputdata.p.dAHPAnova,~,dAHPstats]=anova1([outputdata.dAHP.(fldLabels{1}) outputdata.dAHP.(fldLabels{2}) outputdata.dAHP.(fldLabels{3})],...
                [ones(1,length(outputdata.dAHP.(fldLabels{1}))) 2.*ones(1,length(outputdata.dAHP.(fldLabels{2}))) 3.*ones(1,length(outputdata.dAHP.(fldLabels{3})))]);
            if outputdata.p.dAHPAnova <= 0.05
                outputdata.dAHPPostHocTukey=multcompare(dAHPstats); %done with tukey's HSD method
            end
        else %run Kruskal wallis
            [outputdata.p.dAHPKW,~,~]=kruskalwallis([outputdata.dAHP.(fldLabels{1}) outputdata.dAHP.(fldLabels{2}) outputdata.dAHP.(fldLabels{3})],...
                [ones(1,length(outputdata.dAHP.(fldLabels{1}))) 2.*ones(1,length(outputdata.dAHP.(fldLabels{2}))) 3.*ones(1,length(outputdata.dAHP.(fldLabels{3})))]);
            if outputdata.p.dAHPKW <= 0.05 %run MWU between each data set and correct for pvals with FDR
                [outputdata.p.dAHPPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.dAHP.(fldLabels{1}),outputdata.dAHP.(fldLabels{2}));
                [outputdata.p.dAHPPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]),~]=ranksum(outputdata.dAHP.(fldLabels{1}),outputdata.dAHP.(fldLabels{3}));
                [outputdata.p.dAHPPostHocMW.([fldLabels{2} 'vs' fldLabels{3}]),~]=ranksum(outputdata.dAHP.(fldLabels{2}),outputdata.dAHP.(fldLabels{3}));
                outputdata.pFDR.dAHP=drsFDRpval([outputdata.p.dAHPPostHocMW.([fldLabels{1} 'vs' fldLabels{2}]) outputdata.p.dAHPPostHocMW.([fldLabels{1} 'vs' fldLabels{3}]) outputdata.p.dAHPPostHocMW.([fldLabels{2} 'vs' fldLabels{3}])]);
            end
        end
    elseif noGroups < 3
        if hDAHP == 0 %if groups normal, run ttest
            if noGroups == 2
                [~,outputdata.p.dAHPTtest.([fldLabels{1} 'vs' fldLabels{2}])]=ttest2(outputdata.dAHP.(fldLabels{1}),outputdata.dAHP.(fldLabels{2}));
            end
        elseif hDAHP ~= 0 %if data not normal, run MWU test
            if noGroups == 2
                [outputdata.p.dAHPMW.([fldLabels{1} 'vs' fldLabels{2}]),~]=ranksum(outputdata.dAHP.(fldLabels{1}),outputdata.ahpLat.(fldLabels{2}));
            end
        end
    end
end

%% DISPLAY DATA
close all

%trace figures
for ii = 1:noGroups
    if inputs.RampinjectionprotocolCheckBox.Value == 1
        exTraceFig.(fldLabels{ii}).ramp = figure(ii);
        exTraceFig.(fldLabels{ii}).ramp.Position = [25+(ii-1)*370 550 260 250];
        subplot(4,1,1:3)
        hold on
        noSweeps = size(outputdata.rampVoltage.(fldLabels{ii}),2);
        hiSweep = randperm(noSweeps,1);
        for jj = 1:noSweeps
            if jj ~= hiSweep
                if jj == noSweeps
                    plot(t.ramp,outputdata.rampVoltage.(fldLabels{ii})(:,noSweeps),'linewidth',1,'color',inputs.plotColorGrps(ii,:));
                    plot(t.ramp,outputdata.rampVoltage.(fldLabels{ii})(:,hiSweep),'linewidth',2,'color',.666.*inputs.plotColorGrps(ii,:))
                else
                    plot(t.ramp,outputdata.rampVoltage.(fldLabels{ii})(:,jj),'linewidth',1,'color',inputs.plotColorGrps(ii,:));
                end
            end
        end
        setAx(gca)
        xlabel('time (s)')
        ylabel('membrame voltage (mV)')
        subplot(4,1,4)
        plot(t.ramp,outputdata.rampCurrent.(fldLabels{ii})(:,hiSweep),'linewidth',2,'color',.666.*inputs.plotColorGrps(ii,:))
        setAx(gca)
        xlabel('time (s)')
        ylabel('injected current (pA)')
    end
    
    if inputs.SquareinjectionprotocolCheckBox.Value == 1
        exTraceFig.(fldLabels{ii}).square = figure(ii+noGroups);
        exTraceFig.(fldLabels{ii}).square.Position = [490+(ii-1)*265 550 260 250];
        subplot(4,1,1:3)
        hold on
        noSweeps = size(outputdata.RunTwo.(fldLabels{ii}),3);
        for jj = 2:noSweeps
            if (jj ~= outputdata.rheoSweep.(fldLabels{ii})) || (jj ~= (outputdata.rheoSweep.(fldLabels{ii})+2))
                if jj == noSweeps
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,1,noSweeps),'linewidth',1,'color',inputs.plotColorGrps(ii,:));
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,1,1),'linewidth',2,'color',.8.*inputs.plotColorGrps(ii,:))
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,1,outputdata.rheoSweep.(fldLabels{ii})+2),'linewidth',2,'color',.8.*inputs.plotColorGrps(ii,:))
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,1,outputdata.rheoSweep.(fldLabels{ii})),'linewidth',2,'color',.666.*inputs.plotColorGrps(ii,:))
                else
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,1,jj),'linewidth',1,'color',inputs.plotColorGrps(ii,:));
                end
            end
        end
        setAx(gca)
        xlim([.25 1.75])
        xlabel('time (s)')
        ylabel('membrame voltage (mV)')
        subplot(4,1,4)
        hold on
        for jj = 2:noSweeps
            if (jj ~= outputdata.rheoSweep.(fldLabels{ii})) || (jj ~= (outputdata.rheoSweep.(fldLabels{ii})+2))
                if jj == noSweeps
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,2,noSweeps),'linewidth',1,'color',inputs.plotColorGrps(ii,:));
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,2,1),'linewidth',2,'color',.8.*inputs.plotColorGrps(ii,:))
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,2,outputdata.rheoSweep.(fldLabels{ii})+2),'linewidth',2,'color',.8.*inputs.plotColorGrps(ii,:))
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,2,outputdata.rheoSweep.(fldLabels{ii})),'linewidth',2,'color',.666.*inputs.plotColorGrps(ii,:))
                else
                    plot(t.square,outputdata.RunTwo.(fldLabels{ii})(:,2,jj),'linewidth',1,'color',inputs.plotColorGrps(ii,:));
                end
            end
        end
        setAx(gca)
        xlabel('time (s)')
        ylabel('injected current (pA)')
    end
end

if inputs.RampinjectionprotocolCheckBox.Value == 1
    %vrest figure
    vrestFig=figure(noGroups.*2+1);
    vrestFig.Position = [100 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.vRest.(fldLabels{ii}))),outputdata.vRest.(fldLabels{ii}),40,lightgray,'filled')
        if hRest == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.vRest.(fldLabels{ii})),sem(outputdata.vRest.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.vRest.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hRest ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.vRest.(fldLabels{ii}),upperError.outputdata.vRest.(fldLabels{ii})] = iqrError(outputdata.vRest.(fldLabels{ii}),1);
            errorbar(ii,median(outputdata.vRest.(fldLabels{ii})),lowerError.outputdata.vRest.(fldLabels{ii}),upperError.outputdata.vRest.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.vRest.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allVrest = outputdata.vRest.(fldLabels{ii})';
        else
            allVrest = [allVrest; outputdata.vRest.(fldLabels{ii})'];
        end
    end
    ylabel('V_{rest} (mV)')
    vrestAx=gca;
    setAx(vrestAx);
    vrestAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([1.2*min(allVrest) 0])
    
    %rheobase figure
    rheobaseFig=figure(noGroups.*2+2);
    rheobaseFig.Position = [240 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.rheobase.(fldLabels{ii}))),outputdata.rheobase.(fldLabels{ii}),40,lightgray,'filled')
        if hRheobase == 0 %if data are norheobaseal, plot mean w sem
            errorbar(ii,mean(outputdata.rheobase.(fldLabels{ii})),sem(outputdata.rheobase.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.rheobase.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hRheobase ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.rheobase.(fldLabels{ii}),upperError.outputdata.rheobase.(fldLabels{ii})] = iqrError(outputdata.rheobase.(fldLabels{ii}),1);
            errorbar(ii,median(outputdata.rheobase.(fldLabels{ii})),lowerError.outputdata.rheobase.(fldLabels{ii}),upperError.outputdata.rheobase.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.rheobase.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allRheobase = outputdata.rheobase.(fldLabels{ii})';
        else
            allRheobase = [allRheobase; outputdata.rheobase.(fldLabels{ii})'];
        end
    end
    ylabel('Rheobase Current (pA; from rest)')
    rheobaseAx=gca;
    setAx(rheobaseAx);
    rheobaseAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allRheobase)])
end

%threshold figure
thresholdFig=figure(noGroups.*2+3);
thresholdFig.Position = [380 100 133 300];
hold on
for ii = 1:noGroups
    scatter(ii.*ones(1,length(outputdata.threshold.(fldLabels{ii}))),outputdata.threshold.(fldLabels{ii}),40,lightgray,'filled')
    if hThreshold == 0 %if data are nothresholdal, plot mean w sem
        errorbar(ii,mean(outputdata.threshold.(fldLabels{ii})),sem(outputdata.threshold.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,mean(outputdata.threshold.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    elseif hThreshold ~= 0 %if data non-normal, plot median w IQR
        [lowerError.outputdata.threshold.(fldLabels{ii}),upperError.outputdata.threshold.(fldLabels{ii})] = iqrError(outputdata.threshold.(fldLabels{ii}),2);
        errorbar(ii,median(outputdata.threshold.(fldLabels{ii})),lowerError.outputdata.threshold.(fldLabels{ii}),upperError.outputdata.threshold.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        scatter(ii,median(outputdata.threshold.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
    end
    if ii == 1
        allThreshold = outputdata.threshold.(fldLabels{ii})';
    else
        allThreshold = [allThreshold; outputdata.threshold.(fldLabels{ii})'];
    end
end
ylabel('Action Potential Threshold (mV)')
thresholdAx=gca;
setAx(thresholdAx);
thresholdAx.XTick=[];
xlim([0.5 noGroups+.5])
ylim([1.2*min(allThreshold) 0])


if inputs.SquareinjectionprotocolCheckBox.Value == 1
    %rm figure
    rmFig=figure(noGroups.*2+4);
    rmFig.Position = [520 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.rm.(fldLabels{ii}))),outputdata.rm.(fldLabels{ii}),40,lightgray,'filled')
        if hRm == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.rm.(fldLabels{ii})),sem(outputdata.rm.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.rm.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hRm ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.rm.(fldLabels{ii}),upperError.outputdata.rm.(fldLabels{ii})] = iqrError(outputdata.rm.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.rm.(fldLabels{ii})),lowerError.outputdata.rm.(fldLabels{ii}),upperError.outputdata.rm.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.rm.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allRm = outputdata.rm.(fldLabels{ii})';
        else
            allRm = [allRm; outputdata.rm.(fldLabels{ii})'];
        end
    end
    ylabel('Membrane Resistance (M\Omega)')
    rmAx=gca;
    setAx(rmAx);
    rmAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allRm)])
    
    %membrane decay figure
    mtauFig=figure(noGroups.*2+5);
    mtauFig.Position = [660 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.mtau.(fldLabels{ii}))),outputdata.mtau.(fldLabels{ii}),40,lightgray,'filled')
        if  hMemTau == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.mtau.(fldLabels{ii})),sem(outputdata.mtau.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.mtau.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hMemTau ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.mtau.(fldLabels{ii}),upperError.outputdata.mtau.(fldLabels{ii})] = iqrError(outputdata.mtau.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.mtau.(fldLabels{ii})),lowerError.outputdata.mtau.(fldLabels{ii}),upperError.outputdata.mtau.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.mtau.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allTau = outputdata.mtau.(fldLabels{ii})';
        else
            allTau = [allTau; outputdata.mtau.(fldLabels{ii})'];
        end
    end
    ylabel('Membrane Decay \tau (ms)')
    mtauAx=gca;
    setAx(mtauAx);
    mtauAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allTau)])
    
    %voltage sag figure
    vsagFig=figure(noGroups.*2+6);
    vsagFig.Position = [800 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.vsag.(fldLabels{ii}))),outputdata.vsag.(fldLabels{ii}),40,lightgray,'filled')
        if  hSag == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.vsag.(fldLabels{ii})),sem(outputdata.vsag.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.vsag.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hSag ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.vsag.(fldLabels{ii}),upperError.outputdata.vsag.(fldLabels{ii})] = iqrError(outputdata.vsag.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.vsag.(fldLabels{ii})),lowerError.outputdata.vsag.(fldLabels{ii}),upperError.outputdata.vsag.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.vsag.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allSag = outputdata.vsag.(fldLabels{ii})';
        else
            allSag = [allSag; outputdata.vsag.(fldLabels{ii})'];
        end
    end
    ylabel('Voltage Sag (%)')
    vsagAx=gca;
    setAx(vsagAx);
    vsagAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allSag)])
    
    %rebound spikes figure
    rbndFig=figure(noGroups.*2+7);
    rbndFig.Position = [940 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.nRbnd.(fldLabels{ii}))),outputdata.nRbnd.(fldLabels{ii}),40,lightgray,'filled')
        if  hRebound == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.nRbnd.(fldLabels{ii})),sem(outputdata.nRbnd.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.nRbnd.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hRebound ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.nRbnd.(fldLabels{ii}),upperError.outputdata.nRbnd.(fldLabels{ii})] = iqrError(outputdata.nRbnd.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.nRbnd.(fldLabels{ii})),lowerError.outputdata.nRbnd.(fldLabels{ii}),upperError.outputdata.nRbnd.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.nRbnd.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allRbnd = outputdata.nRbnd.(fldLabels{ii})';
        else
            allRbnd = [allRbnd; outputdata.nRbnd.(fldLabels{ii})'];
        end
    end
    ylabel('# Rebound Spikes)')
    rbndAx=gca;
    setAx(rbndAx);
    rbndAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allRbnd)])
    
    %max firing figure
    maxFRFig=figure(noGroups.*2+8);
    maxFRFig.Position = [1100 100 350 300];
    subplot(1,3,3)
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.maxFR.(fldLabels{ii}))),outputdata.maxFR.(fldLabels{ii}),40,lightgray,'filled')
        if  hMaxFR == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.maxFR.(fldLabels{ii})),sem(outputdata.maxFR.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.maxFR.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hMaxFR ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.maxFR.(fldLabels{ii}),upperError.outputdata.maxFR.(fldLabels{ii})] = iqrError(outputdata.maxFR.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.maxFR.(fldLabels{ii})),lowerError.outputdata.maxFR.(fldLabels{ii}),upperError.outputdata.maxFR.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.maxFR.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allFR = outputdata.maxFR.(fldLabels{ii})';
        else
            allFR = [allFR; outputdata.maxFR.(fldLabels{ii})'];
        end
    end
    ylabel('Max Firing Rate (Hz)')
    maxFRAx=gca;
    setAx(maxFRAx);
    maxFRAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allFR)])
    
    subplot(1,3,1:2)
    hold on
    for ii = 1:noGroups
        if  hMaxFR == 0 %if data are normal, plot mean w sem
            errorbar(outputdata.fiCurve.(fldLabels{ii}).currentInjection,mean(outputdata.fiCurve.(fldLabels{ii}).firingRate),sem(outputdata.fiCurve.(fldLabels{ii}).firingRate,1),...
                'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        elseif hMaxFR ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.fiCurve.(fldLabels{ii}),upperError.outputdata.fiCurve.(fldLabels{ii})] = iqrError(outputdata.fiCurve.(fldLabels{ii}).firingRate,1);
            errorbar(outputdata.fiCurve.(fldLabels{ii}).currentInjection,median(outputdata.fiCurve.(fldLabels{ii}).firingRate),lowerError.outputdata.fiCurve.(fldLabels{ii}),upperError.outputdata.fiCurve.(fldLabels{ii}),...
                'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
        end
    end
    ylabel('Firing Rate (Hz)')
    xlabel('Current Injection (pA)')
    fiAx=gca;
    setAx(fiAx);
    
    %halfwidth figure
    hwFig=figure(noGroups.*2+9);
    hwFig.Position = [100 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.hw.(fldLabels{ii}))),outputdata.hw.(fldLabels{ii}),40,lightgray,'filled')
        if hHalfwidth == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.hw.(fldLabels{ii})),sem(outputdata.hw.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.hw.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hHalfwidth ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.hw.(fldLabels{ii}),upperError.outputdata.hw.(fldLabels{ii})] = iqrError(outputdata.hw.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.hw.(fldLabels{ii})),lowerError.outputdata.hw.(fldLabels{ii}),upperError.outputdata.hw.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.hw.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allHalfwidth = outputdata.hw.(fldLabels{ii})';
        else
            allHalfwidth = [allHalfwidth; outputdata.hw.(fldLabels{ii})'];
        end
    end
    ylabel('Spike Halfwidth (ms)')
    hwAx=gca;
    setAx(hwAx);
    hwAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allHalfwidth)])
    
    %spike amp figure
    spkAmpFig=figure(noGroups.*2+10);
    spkAmpFig.Position = [240 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.spkAmp.(fldLabels{ii}))),outputdata.spkAmp.(fldLabels{ii}),40,lightgray,'filled')
        if hSpkAmp == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.spkAmp.(fldLabels{ii})),sem(outputdata.spkAmp.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.spkAmp.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hSpkAmp ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.spkAmp.(fldLabels{ii}),upperError.outputdata.spkAmp.(fldLabels{ii})] = iqrError(outputdata.spkAmp.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.spkAmp.(fldLabels{ii})),lowerError.outputdata.spkAmp.(fldLabels{ii}),upperError.outputdata.spkAmp.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.spkAmp.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allSpkAmp = outputdata.spkAmp.(fldLabels{ii})';
        else
            allSpkAmp = [allSpkAmp; outputdata.spkAmp.(fldLabels{ii})'];
        end
    end
    ylabel('Spike Amplitude (mV)')
    spkAmpAx=gca;
    setAx(spkAmpAx);
    spkAmpAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allSpkAmp)])
    
    %ahp amp figure
    ahpValFig=figure(noGroups.*2+11);
    ahpValFig.Position = [380 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.ahpVal.(fldLabels{ii}))),outputdata.ahpVal.(fldLabels{ii}),40,lightgray,'filled')
        if hValAHP == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.ahpVal.(fldLabels{ii})),sem(outputdata.ahpVal.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.ahpVal.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hValAHP ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.ahpVal.(fldLabels{ii}),upperError.outputdata.ahpVal.(fldLabels{ii})] = iqrError(outputdata.ahpVal.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.ahpVal.(fldLabels{ii})),lowerError.outputdata.ahpVal.(fldLabels{ii}),upperError.outputdata.ahpVal.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.ahpVal.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allAHPVal = outputdata.ahpVal.(fldLabels{ii})';
        else
            allAHPVal = [allAHPVal; outputdata.ahpVal.(fldLabels{ii})'];
        end
    end
    ylabel('AHP Amplitude (mV)')
    ahpValAx=gca;
    setAx(ahpValAx);
    ahpValAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allAHPVal)])
    
    %ahp latency figure
    ahpLatFig=figure(noGroups.*2+12);
    ahpLatFig.Position = [520 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.ahpLat.(fldLabels{ii}))),outputdata.ahpLat.(fldLabels{ii}),40,lightgray,'filled')
        if hLatAHP == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.ahpLat.(fldLabels{ii})),sem(outputdata.ahpLat.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.ahpLat.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hLatAHP ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.ahpLat.(fldLabels{ii}),upperError.outputdata.ahpLat.(fldLabels{ii})] = iqrError(outputdata.ahpLat.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.ahpLat.(fldLabels{ii})),lowerError.outputdata.ahpLat.(fldLabels{ii}),upperError.outputdata.ahpLat.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.ahpLat.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allAHPlat = outputdata.ahpLat.(fldLabels{ii})';
        else
            allAHPlat = [allAHPlat; outputdata.ahpLat.(fldLabels{ii})'];
        end
    end
    ylabel('AHP latency (ms)')
    ahpLatAx=gca;
    setAx(ahpLatAx);
    ahpLatAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allAHPlat)])
    
    %delta AHP
    dAHPFig=figure(noGroups.*2+13);
    dAHPFig.Position = [660 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.dAHP.(fldLabels{ii}))),outputdata.dAHP.(fldLabels{ii}),40,lightgray,'filled')
        if hDAHP == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.dAHP.(fldLabels{ii})),sem(outputdata.dAHP.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.dAHP.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hDAHP ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.dAHP.(fldLabels{ii}),upperError.outputdata.dAHP.(fldLabels{ii})] = iqrError(outputdata.dAHP.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.dAHP.(fldLabels{ii})),lowerError.outputdata.dAHP.(fldLabels{ii}),upperError.outputdata.dAHP.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.dAHP.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allDeltaAHP = outputdata.dAHP.(fldLabels{ii})';
        else
            allDeltaAHP = [allDeltaAHP; outputdata.dAHP.(fldLabels{ii})'];
        end
    end
    ylabel('\DeltaAHP (ms)')
    dAHPAx=gca;
    setAx(dAHPAx);
    dAHPAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([1.2*min(allDeltaAHP) 1.2*max(allDeltaAHP)])
    
    %latency to first AP figure
    spkLatFig=figure(noGroups.*2+14);
    spkLatFig.Position = [800 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.spkLat.(fldLabels{ii}))),outputdata.spkLat.(fldLabels{ii}),40,lightgray,'filled')
        if hSpkLat == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.spkLat.(fldLabels{ii})),sem(outputdata.spkLat.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.spkLat.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hSpkLat ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.spkLat.(fldLabels{ii}),upperError.outputdata.spkLat.(fldLabels{ii})] = iqrError(outputdata.spkLat.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.spkLat.(fldLabels{ii})),lowerError.outputdata.spkLat.(fldLabels{ii}),upperError.outputdata.spkLat.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.spkLat.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allSpkLat = outputdata.spkLat.(fldLabels{ii})';
        else
            allSpkLat = [allSpkLat; outputdata.spkLat.(fldLabels{ii})'];
        end
    end
    ylabel('Latency to first spike (ms)')
    spkLatAx=gca;
    setAx(spkLatAx);
    spkLatAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allSpkLat)])
    
    %firing rate accomadation ratio
    frRatioFig=figure(noGroups.*2+15);
    frRatioFig.Position = [940 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.frRatio.(fldLabels{ii}))),outputdata.frRatio.(fldLabels{ii}),40,lightgray,'filled')
        if hFRA == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.frRatio.(fldLabels{ii})),sem(outputdata.frRatio.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.frRatio.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hFRA ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.frRatio.(fldLabels{ii}),upperError.outputdata.frRatio.(fldLabels{ii})] = iqrError(outputdata.frRatio.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.frRatio.(fldLabels{ii})),lowerError.outputdata.frRatio.(fldLabels{ii}),upperError.outputdata.frRatio.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.frRatio.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allFRRatio = outputdata.frRatio.(fldLabels{ii})';
        else
            allFRRatio = [allFRRatio; outputdata.frRatio.(fldLabels{ii})'];
        end
    end
    ylabel('FR Accomodation Ratio')
    frRatioAx=gca;
    setAx(frRatioAx);
    frRatioAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allFRRatio)])
    
    %amplitude accomdation ratio
    ampRatioFig=figure(noGroups.*2+16);
    ampRatioFig.Position = [1080 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.ampRatio.(fldLabels{ii}))),outputdata.ampRatio.(fldLabels{ii}),40,lightgray,'filled')
        if hSAA == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.ampRatio.(fldLabels{ii})),sem(outputdata.ampRatio.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.ampRatio.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hSAA ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.ampRatio.(fldLabels{ii}),upperError.outputdata.ampRatio.(fldLabels{ii})] = iqrError(outputdata.ampRatio.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.ampRatio.(fldLabels{ii})),lowerError.outputdata.ampRatio.(fldLabels{ii}),upperError.outputdata.ampRatio.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.ampRatio.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allSARatio = outputdata.ampRatio.(fldLabels{ii})';
        else
            allSARatio = [allSARatio; outputdata.ampRatio.(fldLabels{ii})'];
        end
    end
    ylabel('Spike Amplitude Accomodation Ratio')
    ampRatioAx=gca;
    setAx(ampRatioAx);
    ampRatioAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allSARatio)])
    
    %spike broadening ratio
    broadeningRatioFig=figure(noGroups.*2+17);
    broadeningRatioFig.Position = [1220 100 133 300];
    hold on
    for ii = 1:noGroups
        scatter(ii.*ones(1,length(outputdata.broadeningRatio.(fldLabels{ii}))),outputdata.broadeningRatio.(fldLabels{ii}),40,lightgray,'filled')
        if hBA == 0 %if data are normal, plot mean w sem
            errorbar(ii,mean(outputdata.broadeningRatio.(fldLabels{ii})),sem(outputdata.broadeningRatio.(fldLabels{ii}),1),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,mean(outputdata.broadeningRatio.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        elseif hBA ~= 0 %if data non-normal, plot median w IQR
            [lowerError.outputdata.broadeningRatio.(fldLabels{ii}),upperError.outputdata.broadeningRatio.(fldLabels{ii})] = iqrError(outputdata.broadeningRatio.(fldLabels{ii}),2);
            errorbar(ii,median(outputdata.broadeningRatio.(fldLabels{ii})),lowerError.outputdata.broadeningRatio.(fldLabels{ii}),upperError.outputdata.broadeningRatio.(fldLabels{ii}),'color',inputs.plotColorGrps(ii,:),'linewidth',3,'CapSize',0)
            scatter(ii,median(outputdata.broadeningRatio.(fldLabels{ii})),125,inputs.plotColorGrps(ii,:),'filled')
        end
        if ii == 1
            allBroadeningRatio = outputdata.broadeningRatio.(fldLabels{ii})';
        else
            allBroadeningRatio = [allBroadeningRatio; outputdata.broadeningRatio.(fldLabels{ii})'];
        end
    end
    ylabel('Spike Broadening Ratio')
    broadeningRatioAx=gca;
    setAx(broadeningRatioAx);
    broadeningRatioAx.XTick=[];
    xlim([0.5 noGroups+.5])
    ylim([0 1.2*max(allBroadeningRatio)])
    
end

%% REPORT RESULTS
outputdata
outputdata.p

end
