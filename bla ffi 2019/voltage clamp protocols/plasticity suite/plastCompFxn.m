function [output] = plastCompFxn(inputs)
%for group comparisons of plasticity data
%modeled off ltpVC.m code
%created 04-21-2018, edited 10-19-18

%% INIT VARS
t.pre = -1*inputs.baseDurationminEditField.Value+.5:1:-.5;
t.post = 1.75:1:inputs.postDurationminEditField.Value+.75;
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
cd(inputs.plastDir);

for ii = 1:noGroups
    %get files,
    cd([inputs.plastDir,'/',grpLabels{ii},'/plasticity']);
    contents.(fldLabels{ii}) = dir('*.mat');
    filenames.(fldLabels{ii}) = {contents.(fldLabels{ii}).name}';
    plastFiles.(fldLabels{ii}) = fullfile(cd,filenames.(fldLabels{ii}));
    
    %init
    epscNormAmp.(fldLabels{ii}).pre = -1.*ones(length(plastFiles.(fldLabels{ii})),20); %row == cell number, column == stim intensity
    epscNormAmp.(fldLabels{ii}).post = -1.*ones(length(plastFiles.(fldLabels{ii})),120); %row == cell number, column == stim intensity
    epscAmp.(fldLabels{ii}).pre = -1.*ones(length(plastFiles.(fldLabels{ii})),20); %row == cell number, column == stim intensity
    epscAmp.(fldLabels{ii}).post = -1.*ones(length(plastFiles.(fldLabels{ii})),120); %row == cell number, column == stim intensity
    h_preAmp.(fldLabels{ii}) = -1.*ones(length(plastFiles.(fldLabels{ii})),1);
    h_postAmp.(fldLabels{ii}) = -1.*ones(length(plastFiles.(fldLabels{ii})),1);
    p_plastic_tt.(fldLabels{ii}) = -1.*ones(length(plastFiles.(fldLabels{ii})),1);
    p_plastic_mwu.(fldLabels{ii}) = -1.*ones(length(plastFiles.(fldLabels{ii})),1);
    
    for jj = 1:length(plastFiles.(fldLabels{ii}))
        %load data
        load(plastFiles.(fldLabels{ii}){jj}, 'app')
        epscAmp.(fldLabels{ii}).pre(jj,:) = app.plastData.amp.pre;
        epscAmp.(fldLabels{ii}).post(jj,:) = app.plastData.amp.post;
        epscNormAmp.(fldLabels{ii}).pre(jj,:) = app.plastData.normAmp.pre;
        epscNormAmp.(fldLabels{ii}).post(jj,:) = app.plastData.normAmp.post;
        
        %test normality, p-value of response
        thesePre = app.plastData.amp.pre;
        thesePre(isnan(thesePre)) = [];
        h_preAmp.(fldLabels{ii})(jj) = adtest(thesePre);
        h_postAmp.(fldLabels{ii})(jj) = adtest(app.plastData.lastFive);
        [~,p_plastic_tt.(fldLabels{ii})(jj)] = ttest2(thesePre,app.plastData.lastFive);
        p_plastic_mwu.(fldLabels{ii})(jj) = ranksum(thesePre,app.plastData.lastFive);
        clear app thesePre
    end
end

%init change in amplitude data
for ii = 1:noGroups
    deltaAmp.(fldLabels{ii}) = -1.*ones(length(plastFiles.(fldLabels{ii})),1);
    for jj = 1:length(plastFiles.(fldLabels{ii}))
        holder = [];
        holder = epscNormAmp.(fldLabels{ii}).post(jj,end-inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value-1:end);
        deltaAmp.(fldLabels{ii})(jj) = mean(holder(~isnan(holder)));
    end
end

%% check for normality of event distributions
%if any cell has a pre or post that is not normal, do all test as mwu
%if all pre and post are normal, do ttest
for ii = 1:noGroups
    if sum([h_preAmp.(fldLabels{ii}); h_postAmp.(fldLabels{ii})]) == 0
        output.thesePvals.(fldLabels{ii}) = p_plastic_tt.(fldLabels{ii});
        output.alphaFDR_tt.(fldLabels{ii}) = drsFDRpval(output.thesePvals.(fldLabels{ii}));
        output.sigPlastic.(fldLabels{ii}) = p_plastic_tt.(fldLabels{ii}) < output.alphaFDR_tt.(fldLabels{ii});
    else
        output.thesePvals.(fldLabels{ii}) = p_plastic_mwu.(fldLabels{ii});
        output.alphaFDR_mwu.(fldLabels{ii}) = drsFDRpval(output.thesePvals.(fldLabels{ii}));
        output.sigPlastic.(fldLabels{ii}) = p_plastic_mwu.(fldLabels{ii}) < output.alphaFDR_mwu.(fldLabels{ii});
    end
end

%% test for significance of plastic changes across groups
% ie do some conditions more likely to undergo plastic changes than others
%chi-sq tests
if noGroups == 2
    [output.pChiSq.([fldLabels{1} 'vs' fldLabels{2}]),~] = chi2test([sum(output.sigPlastic.(fldLabels{1})), length(output.sigPlastic.(fldLabels{1}))-sum(output.sigPlastic.(fldLabels{1})); sum(output.sigPlastic.(fldLabels{2})), length(output.sigPlastic.(fldLabels{2}))-sum(output.sigPlastic.(fldLabels{2}))]);
elseif noGroups == 3
    [output.pChiSq.([fldLabels{1} 'vs' fldLabels{2}]),~] = chi2test([sum(output.sigPlastic.(fldLabels{1})), length(output.sigPlastic.(fldLabels{1}))-sum(output.sigPlastic.(fldLabels{1})); sum(output.sigPlastic.(fldLabels{2})), length(output.sigPlastic.(fldLabels{2}))-sum(output.sigPlastic.(fldLabels{2}))]); %#ok<*NOPRT>
    [output.pChiSq.([fldLabels{1} 'vs' fldLabels{3}]),~] = chi2test([sum(output.sigPlastic.(fldLabels{1})), length(output.sigPlastic.(fldLabels{1}))-sum(output.sigPlastic.(fldLabels{1})); sum(output.sigPlastic.(fldLabels{3})), length(output.sigPlastic.(fldLabels{3}))-sum(output.sigPlastic.(fldLabels{3}))]);
    [output.pChiSq.([fldLabels{2} 'vs' fldLabels{3}]),~] = chi2test([sum(output.sigPlastic.(fldLabels{2})), length(output.sigPlastic.(fldLabels{2}))-sum(output.sigPlastic.(fldLabels{2})); sum(output.sigPlastic.(fldLabels{3})), length(output.sigPlastic.(fldLabels{3}))-sum(output.sigPlastic.(fldLabels{3}))]);
    %FDR correction
    output.alphaChiSq = drsFDRpval([output.pChiSq.([fldLabels{1} 'vs' fldLabels{2}]) output.pChiSq.([fldLabels{1} 'vs' fldLabels{3}]) output.pChiSq.([fldLabels{2} 'vs' fldLabels{3}])]);
end

%% directionality of plastic change
for ii = 1:noGroups
    %LTD
    ltdMod.(fldLabels{ii}) = deltaAmp.(fldLabels{ii}) < 100;
    sigLTD.(fldLabels{ii}) = output.sigPlastic.(fldLabels{ii}).*ltdMod.(fldLabels{ii});
    ltdAmp.(fldLabels{ii}) = sigLTD.(fldLabels{ii}).*deltaAmp.(fldLabels{ii});
    deleteThis = (ltdAmp.(fldLabels{ii})==0);
    ltdAmp.(fldLabels{ii})(deleteThis) = [];
    deleteThis = [];
    
    %LTP
    ltpMod.(fldLabels{ii}) = deltaAmp.(fldLabels{ii}) > 100;
    sigLTP.(fldLabels{ii}) = output.sigPlastic.(fldLabels{ii}).*ltpMod.(fldLabels{ii});
    ltpAmp.(fldLabels{ii}) = sigLTP.(fldLabels{ii}).*deltaAmp.(fldLabels{ii});
    deleteThis = (ltpAmp.(fldLabels{ii})==0);
    ltpAmp.(fldLabels{ii})(deleteThis) = [];
    deleteThis = [];
    
    %No Direction
    theseCells = (output.sigPlastic.(fldLabels{ii}) == 0);
    ndAmp.(fldLabels{ii}) = deltaAmp.(fldLabels{ii})(theseCells);
    theseCells = [];
end

%% test for equality of changes across groups
if noGroups == 2
    deltaAmp.pooled = [deltaAmp.(fldLabels{1}); deltaAmp.(fldLabels{2})];
    deltaAmp.pooledGroups = [ones(length(deltaAmp.(fldLabels{1})),1); 2.*ones(length(deltaAmp.(fldLabels{2})),1)];
elseif noGroups == 3
    deltaAmp.pooled = [deltaAmp.(fldLabels{1}); deltaAmp.(fldLabels{2}); deltaAmp.(fldLabels{3})];
    deltaAmp.pooledGroups = [ones(length(deltaAmp.(fldLabels{1})),1); 2.*ones(length(deltaAmp.(fldLabels{2})),1);...
        3.*ones(length(deltaAmp.(fldLabels{3})),1)];
end
dAmpResiduals.pooled =[];
% ie do some conditions lead to greater divergence from around 100% (no change) than others?
% % Bartlett test if normal data, Browne-Forsythe if non-normal data (see Jordan et al., 2018 Neuron)
for ii = 1:noGroups
    %test normality
    %raw values for variances
    hVar.(fldLabels{ii}) = adtest(deltaAmp.(fldLabels{ii}));
    %residuals for anova
    dAmpResiduals.pooled = [dAmpResiduals.pooled; mean(deltaAmp.(fldLabels{ii}))-deltaAmp.(fldLabels{ii})];
end
hResiduals = adtest(dAmpResiduals.pooled);
if noGroups == 2
    if sum([hVar.(fldLabels{1}) hVar.(fldLabels{2})]) == 0 %if normal...
    else %if non-normal...
    end
elseif noGroups == 3
    if sum([hVar.(fldLabels{1}) hVar.(fldLabels{2}) hVar.(fldLabels{3})]) == 0 %if normal...
        output.pBartlett = vartestn(deltaAmp.pooled,deltaAmp.pooledGroups); %variance testing
        if hResiduals == 0 %if data normal
            [output.pAnova, ~, anovaStats] = anova1(abs(100-deltaAmp.pooled),deltaAmp.pooledGroups); %delta Amplitude testing
            if output.pAnova <= .05 %run multiple comparisons
                output.anovaMC = multcompare(anovaStats);
            end
        end
    else %if non-normal...
    end
end

close all

%% cumulative distribution data
for ii = 1:noGroups
    [ProbAmp.(fldLabels{ii}), ampVal.(fldLabels{ii})] = ecdf(deltaAmp.(fldLabels{ii}));
end

%% time course data
for ii = 1:noGroups
    %pre
    %get data per min (average the 4 events for that minute)
    tAmp.(fldLabels{ii}).pre.all = -1*ones(length(plastFiles.(fldLabels{ii})),5);
    for jj = 1:length(plastFiles.(fldLabels{ii}))
        for kk = 1:inputs.AnalysisdurationminEditField.Value
            hld = [];
            hld = epscNormAmp.(fldLabels{ii}).pre(jj,kk*inputs.StimperminEditField.Value-(inputs.StimperminEditField.Value-1):kk*inputs.StimperminEditField.Value);
            if sum(isnan(hld)) > 0
                hld(isnan(hld))=[];
            end
            if isempty(hld) ~= 1
                tAmp.(fldLabels{ii}).pre.all(jj,kk) = mean(hld);
            end
        end
    end
    
    %inits
    tAmp.(fldLabels{ii}).pre.ltp = -1*ones(sum(sigLTP.(fldLabels{ii})),inputs.AnalysisdurationminEditField.Value);
    tAmp.(fldLabels{ii}).pre.ltd = -1*ones(sum(sigLTD.(fldLabels{ii})),inputs.AnalysisdurationminEditField.Value);
    tAmp.(fldLabels{ii}).pre.nd = -1*ones(sum(output.sigPlastic.(fldLabels{ii})==0),inputs.AnalysisdurationminEditField.Value);
    LTPcount = 0;
    LTDcount = 0;
    NDcount = 0;
    for jj = 1:length(plastFiles.(fldLabels{ii}))
        %ltp only
        if sigLTP.(fldLabels{ii})(jj) == 1
            tAmp.(fldLabels{ii}).pre.ltp(1+LTPcount,:) = tAmp.(fldLabels{ii}).pre.all(jj,:);
            LTPcount = LTPcount +1;
        end
        
        %ltd only
        if sigLTD.(fldLabels{ii})(jj) == 1
            tAmp.(fldLabels{ii}).pre.ltd(1+LTDcount,:) = tAmp.(fldLabels{ii}).pre.all(jj,:);
            LTDcount = LTDcount +1;
        end
        
        %no change
        if output.sigPlastic.(fldLabels{ii})(jj) == 0
            tAmp.(fldLabels{ii}).pre.nd(1+NDcount,:) = tAmp.(fldLabels{ii}).pre.all(jj,:);
            NDcount = NDcount +1;
        end
    end
    
    %get mean, SEM
    for jj = 1:inputs.AnalysisdurationminEditField.Value
        %all
        theseData = tAmp.(fldLabels{ii}).pre.all(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).pre.all(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).pre.all(jj) = sem(theseData,1);
        theseData = [];
        
        %ltp
        theseData = tAmp.(fldLabels{ii}).pre.ltp(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).pre.ltp(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).pre.ltp(jj) = sem(theseData,1);
        theseData = [];
        
        %ltd
        theseData = tAmp.(fldLabels{ii}).pre.ltd(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).pre.ltd(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).pre.ltd(jj) = sem(theseData,1);
        theseData = [];
        
        %nd
        theseData = tAmp.(fldLabels{ii}).pre.nd(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).pre.nd(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).pre.nd(jj) = sem(theseData,1);
        theseData = [];
        
    end
    
    %post
    %all
    tAmp.(fldLabels{ii}).post.all = -1*ones(length(plastFiles.(fldLabels{ii})),30);
    for jj = 1:length(plastFiles.(fldLabels{ii}))
        for kk = 1:inputs.postDurationminEditField.Value
            hld = [];
            hld = epscNormAmp.(fldLabels{ii}).post(jj,kk*inputs.StimperminEditField.Value-(inputs.StimperminEditField.Value-1):kk*inputs.StimperminEditField.Value);
            if sum(isnan(hld)) > 0
                hld(isnan(hld))=[];
            end
            if isempty(hld) ~= 1
                tAmp.(fldLabels{ii}).post.all(jj,kk) = mean(hld);
            end
        end
    end
    
    %inits
    tAmp.(fldLabels{ii}).post.ltp = -1*ones(sum(sigLTP.(fldLabels{ii})),inputs.postDurationminEditField.Value);
    tAmp.(fldLabels{ii}).post.ltd = -1*ones(sum(sigLTD.(fldLabels{ii})),inputs.postDurationminEditField.Value);
    tAmp.(fldLabels{ii}).post.nd = -1*ones(sum(output.sigPlastic.(fldLabels{ii})==0),inputs.postDurationminEditField.Value);
    LTPcount = 0;
    LTDcount = 0;
    NDcount = 0;
    for jj = 1:length(plastFiles.(fldLabels{ii}))
        %ltp only
        if sigLTP.(fldLabels{ii})(jj) == 1
            tAmp.(fldLabels{ii}).post.ltp(1+LTPcount,:) = tAmp.(fldLabels{ii}).post.all(jj,:);
            LTPcount = LTPcount +1;
        end
        
        %ltd only
        if sigLTD.(fldLabels{ii})(jj) == 1
            tAmp.(fldLabels{ii}).post.ltd(1+LTDcount,:) = tAmp.(fldLabels{ii}).post.all(jj,:);
            LTDcount = LTDcount +1;
        end
        
        %no change
        if output.sigPlastic.(fldLabels{ii})(jj) == 0
            tAmp.(fldLabels{ii}).post.nd(1+NDcount,:) = tAmp.(fldLabels{ii}).post.all(jj,:);
            NDcount = NDcount +1;
        end
    end
    
    %get mean, 95CI
    for jj = 1:inputs.postDurationminEditField.Value
        %all
        theseData = tAmp.(fldLabels{ii}).post.all(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).post.all(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).post.all(jj) = sem(theseData,1);
        theseData = [];
        
        %ltp
        theseData = tAmp.(fldLabels{ii}).post.ltp(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).post.ltp(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).post.ltp(jj) = sem(theseData,1);
        theseData = [];
        
        %ltd
        theseData = tAmp.(fldLabels{ii}).post.ltd(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).post.ltd(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).post.ltd(jj) = sem(theseData,1);
        theseData = [];
        
        %nd
        theseData = tAmp.(fldLabels{ii}).post.nd(:,jj);
        deleteThis = (theseData == -1);
        theseData(deleteThis) =[];
        deleteThis = [];
        
        mn.tAmp.(fldLabels{ii}).post.nd(jj) = mean(theseData);
        sterror.tAmp.(fldLabels{ii}).post.nd(jj) = sem(theseData,1);
        theseData = [];
        
    end
end

%% DISPLAY DATA

%response direction display fig
for ii = 1:noGroups
    pieFig.(fldLabels{ii}) = figure(1+ii);
    pieFig.(fldLabels{ii}).Position = [765+(ii-1).*230 550 225 250];
    if sum(sigLTP.(fldLabels{ii})) == 0
        if sum(sigLTD.(fldLabels{ii})) == 0
            pie([.00000001 .00000001 (length(deltaAmp.(fldLabels{ii}))-(sum(sigLTP.(fldLabels{ii}))+sum(sigLTD.(fldLabels{ii}))))],{'LTP','LTD','No plasticity'})
        elseif sum(sigLTD.(fldLabels{ii})) ~= 0
            pie([.00000001 sum(sigLTD.(fldLabels{ii})) (length(deltaAmp.(fldLabels{ii}))-(sum(sigLTP.(fldLabels{ii}))+sum(sigLTD.(fldLabels{ii}))))],{'LTP','LTD','No plasticity'})
        end
    elseif sum(sigLTP.(fldLabels{ii})) ~= 0
        if sum(sigLTD.(fldLabels{ii})) == 0
            pie([sum(sigLTP.(fldLabels{ii})) .00000001 (length(deltaAmp.(fldLabels{ii}))-(sum(sigLTP.(fldLabels{ii}))+sum(sigLTD.(fldLabels{ii}))))],{'LTP','LTD','No plasticity'})
        elseif sum(sigLTD.(fldLabels{ii})) ~= 0
            pie([sum(sigLTP.(fldLabels{ii})) sum(sigLTD.(fldLabels{ii})) (length(deltaAmp.(fldLabels{ii}))-(sum(sigLTP.(fldLabels{ii}))+sum(sigLTD.(fldLabels{ii}))))],{'LTP','LTD','No plasticity'})
        end
    end
    if ii == 1
        colormap([inputs.plotColorGrpOne; (1+(.85-max(inputs.plotColorGrpOne))./max(inputs.plotColorGrpOne)).*inputs.plotColorGrpOne; .33 .33 .33])
    elseif ii == 2
        colormap([inputs.plotColorGrpTwo; (1+(.85-max(inputs.plotColorGrpTwo))./max(inputs.plotColorGrpTwo)).*inputs.plotColorGrpTwo; .33 .33 .33])
    elseif ii == 3
        colormap([inputs.plotColorGrpThree; (1+(.85-max(inputs.plotColorGrpThree))./max(inputs.plotColorGrpThree)).*inputs.plotColorGrpThree; .33 .33 .33])
    end
    
    title(fldLabels{ii})
end

%set up colors
lightgray=[.75 .75 .75]; %light gray for individual data points
for ii = 1:noGroups
    if ii == 1
        c.(fldLabels{ii}).dark = inputs.plotColorGrpOne;
        c.(fldLabels{ii}).light = (1+(.85-max(inputs.plotColorGrpOne))./max(inputs.plotColorGrpOne)).*inputs.plotColorGrpOne;
    elseif ii == 2
        c.(fldLabels{ii}).dark = inputs.plotColorGrpTwo;
        c.(fldLabels{ii}).light = (1+(.85-max(inputs.plotColorGrpTwo))./max(inputs.plotColorGrpTwo)).*inputs.plotColorGrpTwo;
    elseif ii == 3
        c.(fldLabels{ii}).dark = inputs.plotColorGrpThree;
        c.(fldLabels{ii}).light = (1+(.85-max(inputs.plotColorGrpThree))./max(inputs.plotColorGrpThree)).*inputs.plotColorGrpThree;
    end
end

%time course fig, ltp
for ii = 1:noGroups
    %time course
    tFig.(fldLabels{ii}).ltp = figure(noGroups+1+ii);
    tFig.(fldLabels{ii}).ltp.Position = [50+(ii-1).*455 370 320 125];
    hold on
    line([t.pre(1)-1 t.post(end)+1],[100 100],'linestyle','--','color','k')
    if size(tAmp.(fldLabels{ii}).pre.ltp,1) >= 3
        scatter(t.pre,mn.tAmp.(fldLabels{ii}).pre.ltp,50,c.(fldLabels{ii}).dark,'filled')
        scatter(t.post,mn.tAmp.(fldLabels{ii}).post.ltp,50,c.(fldLabels{ii}).dark,'filled')
        errorbar(t.pre,mn.tAmp.(fldLabels{ii}).pre.ltp,sterror.tAmp.(fldLabels{ii}).pre.ltp,'color',c.(fldLabels{ii}).dark,'linewidth',1,'capsize',0,'linestyle','none')
        errorbar(t.post,mn.tAmp.(fldLabels{ii}).post.ltp,sterror.tAmp.(fldLabels{ii}).post.ltp,'color',c.(fldLabels{ii}).dark,'linewidth',1,'capsize',0,'linestyle','none')
    elseif size(tAmp.(fldLabels{ii}).pre.ltp,1) > 0
        for jj = 1:size(tAmp.(fldLabels{ii}).pre.ltp,1)
            if jj == 1
                scatter(t.pre,tAmp.(fldLabels{ii}).pre.ltp(jj,:),50,c.(fldLabels{ii}).dark,'filled')
                scatter(t.post,tAmp.(fldLabels{ii}).post.ltp(jj,:),50,c.(fldLabels{ii}).dark,'filled')
            elseif jj == 2
                scatter(t.pre,tAmp.(fldLabels{ii}).pre.ltp(jj,:),50,c.(fldLabels{ii}).light,'filled')
                scatter(t.post,tAmp.(fldLabels{ii}).post.ltp(jj,:),50,c.(fldLabels{ii}).light,'filled')
            end
        end
    end
    xlim([t.pre(1)-1 t.post(end)+1])
    ylim([0 200])
    tAx = gca;
    tAx.Box = 'off'; tAx.YColor = 'k'; tAx.XColor = 'k'; tAx.LineWidth = 1; tAx.TickDir ='out';
    tAx.XAxisLocation = 'origin';
    xlabel('time (min)')
    ylabel('normalized EPSC amplitude')
    
    %delta EPSC
    daFig.(fldLabels{ii}).ltp = figure(2*noGroups+1+ii);
    daFig.(fldLabels{ii}).ltp.Position = [371+(ii-1).*455 370 125 125];
    hold on
    line([.5 1.5],[100 100],'linestyle','--','color','k')
    if size(tAmp.(fldLabels{ii}).pre.ltp,1) >= 3
        scatter(ones(1,sum(sigLTP.(fldLabels{ii}))),ltpAmp.(fldLabels{ii}),50,c.(fldLabels{ii}).dark,'filled')
    elseif size(tAmp.(fldLabels{ii}).pre.ltp,1) > 0
        for jj = 1:size(tAmp.(fldLabels{ii}).pre.ltp,1)
            if jj == 1
                scatter(1,ltpAmp.(fldLabels{ii})(1),50,c.(fldLabels{ii}).dark,'filled')
            elseif jj == 2
                scatter(1,ltpAmp.(fldLabels{ii})(2),50,c.(fldLabels{ii}).light,'filled')
            end
        end
    end
    ylim([0 200])
    daAx = gca;
    daAx.Box = 'off'; daAx.YColor = 'k'; daAx.XColor = 'k'; daAx.LineWidth = 1; daAx.TickDir ='out'; daAx.XTick = [];
    ylabel('normalized EPSC amplitude (LTP cells), last 5 min')
end

%time course fig, ltd
for ii = 1:noGroups
    %time course
    tFig.(fldLabels{ii}).ltd = figure(3*noGroups+1+ii);
    tFig.(fldLabels{ii}).ltd.Position = [50+(ii-1).*455 180 320 125];
    hold on
    line([t.pre(1)-1 t.post(end)+1],[100 100],'linestyle','--','color','k')
    if size(tAmp.(fldLabels{ii}).pre.ltd,1) >= 3
        scatter(t.pre,mean(tAmp.(fldLabels{ii}).pre.ltd,1),50,c.(fldLabels{ii}).dark,'filled')
        scatter(t.post,mean(tAmp.(fldLabels{ii}).post.ltd,1),50,c.(fldLabels{ii}).dark,'filled')
        errorbar(t.pre,mean(tAmp.(fldLabels{ii}).pre.ltd,1),sem(tAmp.(fldLabels{ii}).pre.ltd,1),'color',c.(fldLabels{ii}).dark,'linewidth',1,'capsize',0,'linestyle','none')
        errorbar(t.post,mean(tAmp.(fldLabels{ii}).post.ltd,1),sem(tAmp.(fldLabels{ii}).post.ltd,1),'color',c.(fldLabels{ii}).dark,'linewidth',1,'capsize',0,'linestyle','none')
    elseif size(tAmp.(fldLabels{ii}).pre.ltd,1) > 0
        for jj = 1:size(tAmp.(fldLabels{ii}).pre.ltd,1)
            if jj == 1
                scatter(t.pre,tAmp.(fldLabels{ii}).pre.ltd(jj,:),50,c.(fldLabels{ii}).dark,'filled')
                scatter(t.post,tAmp.(fldLabels{ii}).post.ltd(jj,:),50,c.(fldLabels{ii}).dark,'filled')
            elseif jj == 2
                scatter(t.pre,tAmp.(fldLabels{ii}).pre.ltd(jj,:),50,c.(fldLabels{ii}).light,'filled')
                scatter(t.post,tAmp.(fldLabels{ii}).post.ltd(jj,:),50,c.(fldLabels{ii}).light,'filled')
            end
        end
    end
    xlim([t.pre(1)-1 t.post(end)+1])
    ylim([0 200])
    tAx = gca;
    tAx.Box = 'off'; tAx.YColor = 'k'; tAx.XColor = 'k'; tAx.LineWidth = 1; tAx.TickDir ='out';
    tAx.XAxisLocation = 'origin';
    xlabel('time (min)')
    ylabel('normalized EPSC amplitude')
    
    %delta EPSC
    daFig.(fldLabels{ii}).ltd = figure(4*noGroups+1+ii);
    daFig.(fldLabels{ii}).ltd.Position = [371+(ii-1).*455 180 125 125];
    hold on
    line([.5 1.5],[100 100],'linestyle','--','color','k')
    if size(tAmp.(fldLabels{ii}).pre.ltd,1) >= 3
        scatter(ones(1,sum(sigLTD.(fldLabels{ii}))),ltdAmp.(fldLabels{ii}),50,c.(fldLabels{ii}).dark,'filled')
    elseif size(tAmp.(fldLabels{ii}).pre.ltd,1) > 0
        for jj = 1:size(tAmp.(fldLabels{ii}).pre.ltd,1)
            if jj == 1
                scatter(1,ltdAmp.(fldLabels{ii})(1),50,c.(fldLabels{ii}).dark,'filled')
            elseif jj == 2
                scatter(1,ltdAmp.(fldLabels{ii})(2),50,c.(fldLabels{ii}).light,'filled')
            end
        end
    end
    ylim([0 200])
    daAx = gca;
    daAx.Box = 'off'; daAx.YColor = 'k'; daAx.XColor = 'k'; daAx.LineWidth = 1; daAx.TickDir ='out'; daAx.XTick = [];
    ylabel('normalized EPSC amplitude (LTD cells), last 5 min')
end

%time course fig, no change cells
for ii = 1:noGroups
    %time course
    tFig.(fldLabels{ii}).nd = figure(5*noGroups+1+ii);
    tFig.(fldLabels{ii}).nd.Position = [50+(ii-1).*455 0 320 125];
    hold on
    line([t.pre(1)-1 t.post(end)+1],[100 100],'linestyle','--','color','k')
    if size(tAmp.(fldLabels{ii}).pre.nd,1) >= 3
        scatter(t.pre,mean(tAmp.(fldLabels{ii}).pre.nd,1),50,c.(fldLabels{ii}).dark,'filled')
        scatter(t.post,mean(tAmp.(fldLabels{ii}).post.nd,1),50,c.(fldLabels{ii}).dark,'filled')
        errorbar(t.pre,mean(tAmp.(fldLabels{ii}).pre.nd,1),sem(tAmp.(fldLabels{ii}).pre.nd,1),'color',c.(fldLabels{ii}).dark,'linewidth',1,'capsize',0,'linestyle','none')
        errorbar(t.post,mean(tAmp.(fldLabels{ii}).post.nd,1),sem(tAmp.(fldLabels{ii}).post.nd,1),'color',c.(fldLabels{ii}).dark,'linewidth',1,'capsize',0,'linestyle','none')
    elseif size(tAmp.(fldLabels{ii}).pre.nd,1) > 0
        for jj = 1:size(tAmp.(fldLabels{ii}).pre.nd,1)
            if jj == 1
                scatter(t.pre,tAmp.(fldLabels{ii}).pre.nd(jj,:),50,c.(fldLabels{ii}).dark,'filled')
                scatter(t.post,tAmp.(fldLabels{ii}).post.nd(jj,:),50,c.(fldLabels{ii}).dark,'filled')
            elseif jj == 2
                scatter(t.pre,tAmp.(fldLabels{ii}).pre.nd(jj,:),50,c.(fldLabels{ii}).light,'filled')
                scatter(t.post,tAmp.(fldLabels{ii}).post.nd(jj,:),50,c.(fldLabels{ii}).light,'filled')
            end
        end
    end
    xlim([t.pre(1)-1 t.post(end)+1])
    ylim([0 200])
    tAx = gca;
    tAx.Box = 'off'; tAx.YColor = 'k'; tAx.XColor = 'k'; tAx.LineWidth = 1; tAx.TickDir ='out';
    tAx.XAxisLocation = 'origin';
    xlabel('time (min)')
    ylabel('normalized EPSC amplitude')
    
    %delta EPSC
    daFig.(fldLabels{ii}).nd = figure(6*noGroups+1+ii);
    daFig.(fldLabels{ii}).nd.Position = [371+(ii-1).*455 0 125 125];
    hold on
    line([.5 1.5],[100 100],'linestyle','--','color','k')
    if size(tAmp.(fldLabels{ii}).pre.nd,1) >= 3
        scatter(ones(1,sum((output.sigPlastic.(fldLabels{ii})==0))),ndAmp.(fldLabels{ii}),50,c.(fldLabels{ii}).dark,'filled')
    elseif size(tAmp.(fldLabels{ii}).pre.nd,1) > 0
        for jj = 1:size(tAmp.(fldLabels{ii}).pre.nd,1)
            if jj == 1
                scatter(1,ndAmp.(fldLabels{ii})(1),50,c.(fldLabels{ii}).dark,'filled')
            elseif jj == 2
                scatter(1,ndAmp.(fldLabels{ii})(2),50,c.(fldLabels{ii}).light,'filled')
            end
        end
    end
    ylim([0 200])
    daAx = gca;
    daAx.Box = 'off'; daAx.YColor = 'k'; daAx.XColor = 'k'; daAx.LineWidth = 1; daAx.TickDir ='out'; daAx.XTick = [];
    ylabel('normalized EPSC amplitude (no change cells), last 5 min')
end

%within group histograms
for ii = 1:noGroups
    histFig.(fldLabels{ii}) = figure(7*noGroups+1+ii);
    histFig.(fldLabels{ii}).Position = [45+(ii-1)*150 635 165 100];
    hold on
    edges = 0:10:200;
    xlim([20 180])
    ylim([0 5])
    line([100 100], [0 5],'linestyle','--','linewidth',1,'color','k')
    histogram(deltaAmp.(fldLabels{ii}),edges,'FaceColor',c.(fldLabels{ii}).dark,'EdgeColor',c.(fldLabels{ii}).dark);
    xlabel('Normalized EPSC Amplitude, last 5 min.')
    ylabel('Count')
    histAx = gca;
    histAx.Box = 'off'; histAx.YColor = 'k'; histAx.XColor = 'k'; histAx.LineWidth = 1; histAx.TickDir ='out';
    
end

%percent plastic histogram
perPlastFig = figure(26);
perPlastFig.Position = [515 535 125 200];
hold on
for ii = 1:noGroups
    bar(ii,sum(output.sigPlastic.(fldLabels{ii}))/length(output.sigPlastic.(fldLabels{ii})),'FaceColor',c.(fldLabels{ii}).dark,'EdgeColor',c.(fldLabels{ii}).dark);
end
ylabel('Proportion Plastic')
perPlastAx = gca;
setAx(perPlastAx);
perPlastAx.XTick = [];

%single deltaAmp plot to look at whole population changes following pairing 
% (use as part of data figure to show ALL data rather than just sig changes)
varFig = figure(27);
varFig.Position = [645 535 125 200];
hold on
line([.5 noGroups+.5],[100 100],'linestyle','--','color','k')
for ii = 1:noGroups
    scatter(ii.*ones(1,sum((output.sigPlastic.(fldLabels{ii})==0))),ndAmp.(fldLabels{ii}),50,lightgray,'filled')
    scatter(ii.*ones(1,sum(sigLTD.(fldLabels{ii}))),ltdAmp.(fldLabels{ii}),50,c.(fldLabels{ii}).dark,'filled')
    scatter(ii.*ones(1,sum(sigLTP.(fldLabels{ii}))),ltpAmp.(fldLabels{ii}),50,c.(fldLabels{ii}).dark,'filled')
end
xlim([0.5 noGroups+.5])
ylim([0 200])
varAx = gca;
varAx.Box = 'off'; varAx.YColor = 'k'; varAx.XColor = 'k'; varAx.LineWidth = 1; varAx.TickDir ='out'; varAx.XTick = [];
ylabel('normalized EPSC amplitude (last 5min)')

output

end
