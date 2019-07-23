function [output] = wnGrpBiophysCompFxn(inputs)
%for group comparison of plasticity data
%created 05-05-2018, edited 05-04-18

%% INITS
grpName{1} = 'responders';
grpName{2} = 'nonresponders';
params = {};

%% LOAD DATA
% Directory
cd(inputs.minStimDir);

%load data
for ii = 1:length(grpName)
    %load data
    if ii == 1
        cd([inputs.minStimDir,'/current injections/minStimResponders']);
    elseif ii == 2
        cd([inputs.minStimDir,'/current injections/minStimNonresponders']);
    end
    contents.(grpName{ii}) = dir('*.mat');
    filenames.(grpName{ii}) = {contents.(grpName{ii}).name}';
    bpFiles.(grpName{ii}) = fullfile(cd,filenames.(grpName{ii}));
    
    %fill params cell
    if inputs.RmCheckBox.Value == 1
        if isempty(params) == 1
            params = 'Rm';
        else
            params = [params, {'Rm'}];
        end
    end
    if inputs.sagCheckBox.Value == 1
        if isempty(params) == 1
            params = 'vsag';
        else
            params = [params, {'vsag'}];
        end
    end
    if inputs.maxFRCheckBox.Value == 1
        if isempty(params) == 1
            params = 'maxFR';
        else
            params = [params, {'maxFR'}];
        end
    end
    if inputs.halfwidthCheckBox.Value == 1
        if isempty(params) == 1
            params = 'hw';
        else
            params = [params, {'hw'}];
        end
    end
    if inputs.AHPmagCheckBox.Value == 1
        if isempty(params) == 1
            params = 'rAHPval';
        else
            params = [params, {'rAHPval'}];
        end
    end
    if inputs.FRadaptCheckBox.Value == 1
        if isempty(params) == 1
            params = 'arFR';
        else
            params = [params, {'arFR'}];
        end
    end
    if inputs.APbroadeningCheckBox.Value == 1
        if isempty(params) == 1
            params = 'brAP';
        else
            params = [params, {'brAP'}];
        end
    end
    if inputs.AHPlatCheckBox.Value == 1
        if isempty(params) == 1
            params = 'rAHPlat';
        else
            params = [params, {'rAHPlat'}];
        end
    end
    if inputs.AmpadaptCheckBox.Value == 1
        if isempty(params) == 1
            params = 'arAmp';
        else
            params = [params, {'arAmp'}];
        end
    end
    if inputs.APampCheckBox.Value == 1
        if isempty(params) == 1
            params = 'rAmp';
        else
            params = [params, {'rAmp'}];
        end
    end
    if inputs.APlatencyCheckBox.Value == 1
        if isempty(params) == 1
            params = 'rheo_APlat';
        else
            params = [params, {'rheo_APlat'}];
        end
    end
    if inputs.APthreshCheckBox.Value == 1
        if isempty(params) == 1
            params = 'rThresh';
        else
            params = [params, {'rThresh'}];
        end
    end
    if inputs.deltaAHPCheckBox.Value == 1
        if isempty(params) == 1
            params = 'dAHP';
        else
            params = [params, {'dAHP'}];
        end
    end
    if inputs.mTauCheckBox.Value == 1
        if isempty(params) == 1
            params = 'mtau';
        else
            params = [params, {'mtau'}];
        end
    end
    if inputs.reboundAPsCheckBox.Value == 1
        if isempty(params) == 1
            params = 'no_rspks';
        else
            params = [params, {'no_rspks'}];
        end
    end
    
    %inits
    if inputs.RmCheckBox.Value == 1
        rm.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.sagCheckBox.Value == 1
        sag.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.maxFRCheckBox.Value == 1
        maxfr.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.halfwidthCheckBox.Value == 1
        hwidth.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.AHPmagCheckBox.Value == 1
        ahpMag.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.FRadaptCheckBox.Value == 1
        frAdapt.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.APbroadeningCheckBox.Value == 1
        spkBroade.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.AHPlatCheckBox.Value == 1
        ahpLat.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.AmpadaptCheckBox.Value == 1
        ampAdapt.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.APampCheckBox.Value == 1
        spkAmp.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.APlatencyCheckBox.Value == 1
        firstSpkLat.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.APthreshCheckBox.Value == 1
        apThresh.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.deltaAHPCheckBox.Value == 1
        deltAHP.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.mTauCheckBox.Value == 1
        tau.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    if inputs.reboundAPsCheckBox.Value == 1
        noReboundSpks.(grpName{ii}) = -1.*ones(1,length(bpFiles.(grpName{ii})));
    end
    
    %actually load the data here
    for jj = 1:length(bpFiles.(grpName{ii}))
        for kk = 1:length(params)
            load(bpFiles.(grpName{ii}){jj},params{kk})
        end
        if inputs.RmCheckBox.Value == 1
            rm.(grpName{ii})(jj) = Rm.value; Rm = [];
        end
        if inputs.sagCheckBox.Value == 1
            sag.(grpName{ii})(jj) = vsag; vsag = [];
        end
        if inputs.maxFRCheckBox.Value == 1
            maxfr.(grpName{ii})(jj) = maxFR; maxFR = [];
        end
        if inputs.halfwidthCheckBox.Value == 1
            hwidth.(grpName{ii})(jj) = mean(hw); hw = [];
        end
        if inputs.AHPmagCheckBox.Value == 1
            ahpMag.(grpName{ii})(jj) = mean(rAHPval); rAHPval=[];
        end
        if inputs.FRadaptCheckBox.Value == 1
            frAdapt.(grpName{ii})(jj) = arFR; arFR=[];
        end
        if inputs.APbroadeningCheckBox.Value == 1
            spkBroade.(grpName{ii})(jj) = brAP; brAP=[];
        end
        if inputs.AHPlatCheckBox.Value == 1
            ahpLat.(grpName{ii})(jj) = mean(rAHPlat); rAHPlat=[];
        end
        if inputs.AmpadaptCheckBox.Value == 1
            ampAdapt.(grpName{ii})(jj) = arAmp; arAmp=[];
        end
        if inputs.APampCheckBox.Value == 1
            spkAmp.(grpName{ii})(jj) = mean(rAmp); rAmp=[];
        end
        if inputs.APlatencyCheckBox.Value == 1
            firstSpkLat.(grpName{ii})(jj) = rheo_APlat; rheo_APlat=[];
        end
        if inputs.APthreshCheckBox.Value == 1
            apThresh.(grpName{ii})(jj) = mean(rThresh); rThresh=[];
        end
        if inputs.deltaAHPCheckBox.Value == 1
            deltAHP.(grpName{ii})(jj) = dAHP; dAHP=[];
        end
        if inputs.mTauCheckBox.Value == 1
            tau.(grpName{ii})(jj) = mtau; mtau=[];
        end
        if inputs.reboundAPsCheckBox.Value == 1
            noReboundSpks.(grpName{ii})(jj) = no_rspks; no_rspks=[];
        end
    end
end

%% TEST DIFFERENCES
%adtest (Anderson-Darling Test) to test for normality
% only has Rm, max FR, membrane tau, halfwidth, hyperpolarization induced sag here -- add more in future if program is used for other tests.

%Rm, normality
if inputs.RmCheckBox.Value == 1
    for ii = 1:length(grpName)
        [holdH.(grpName{ii}),~] = adtest(rm.(grpName{ii}));
    end
    hRm = sum([holdH.(grpName{1}) holdH.(grpName{2})]);
end

%Rm, group differences
if inputs.RmCheckBox.Value == 1
    if hRm == 0 %if normal, do ttest
        [~,output.p.rmTtest]=ttest2(rm.(grpName{1}),rm.(grpName{2}));
    elseif hRm ~= 0 %if not normal, do MW U test
        [output.p.rmMW,~]=ranksum(rm.(grpName{1}),rm.(grpName{2}));
    end
end

%max fr, normality
if inputs.maxFRCheckBox.Value == 1
    for ii = 1:length(grpName)
        [holdH.(grpName{ii}),~] = adtest(maxfr.(grpName{ii}));
    end
    hMaxFR = sum([holdH.(grpName{1}) holdH.(grpName{2})]);
end

%max fr, group differences
if inputs.maxFRCheckBox.Value == 1
    if hMaxFR == 0 %if normal, do ttest
        [~,output.p.maxfrTtest]=ttest2(maxfr.(grpName{1}),maxfr.(grpName{2}));
    elseif hMaxFR ~= 0 %if not normal, do MW U test
        [output.p.maxfrMW,~]=ranksum(maxfr.(grpName{1}),maxfr.(grpName{2}));
    end
end

%halfwidth, normality
if inputs.halfwidthCheckBox.Value == 1
    for ii = 1:length(grpName)
        [holdH.(grpName{ii}),~] = adtest(hwidth.(grpName{ii}));
    end
    hHalfwidth = sum([holdH.(grpName{1}) holdH.(grpName{2})]);
end

%halfwidth, group differences
if inputs.halfwidthCheckBox.Value == 1
    if hHalfwidth == 0 %if normal, do ttest
        [~,output.p.hwidthTtest]=ttest2(hwidth.(grpName{1}),hwidth.(grpName{2}));
    elseif hHalfwidth ~= 0 %if not normal, do MW U test
        [output.p.hwidthMW,~]=ranksum(hwidth.(grpName{1}),hwidth.(grpName{2}));
    end
end

%hyperpolarization induced sag, normality
if inputs.sagCheckBox.Value == 1
    for ii = 1:length(grpName)
        [holdH.(grpName{ii}),~] = adtest(sag.(grpName{ii}));
    end
    hSag = sum([holdH.(grpName{1}) holdH.(grpName{2})]);
end

%hyperpolarization induced sag, group differences
if inputs.sagCheckBox.Value == 1
    if hSag == 0 %if normal, do ttest
        [~,output.p.sagTtest]=ttest2(sag.(grpName{1}),sag.(grpName{2}));
    elseif hSag ~= 0 %if not normal, do MW U test
        [output.p.sagMW,~]=ranksum(sag.(grpName{1}),sag.(grpName{2}));
    end
end

%membrane time constant, normality
if inputs.mTauCheckBox.Value == 1
    for ii = 1:length(grpName)
        [holdH.(grpName{ii}),~] = adtest(tau.(grpName{ii}));
    end
    hTau = sum([holdH.(grpName{1}) holdH.(grpName{2})]);
end

%membrane time constant, group differences
if inputs.mTauCheckBox.Value == 1
    if hTau == 0 %if normal, do ttest
        [~,output.p.tauTtest]=ttest2(tau.(grpName{1}),tau.(grpName{2}));
    elseif hTau ~= 0 %if not normal, do MW U test
        [output.p.tauMW,~]=ranksum(tau.(grpName{1}),tau.(grpName{2}));
    end
end

%% DISPLAY DATA
%colors
lightgray = [.75 .75 .75]; %light gray for individual data points
nrColor = [.33 .33 .33];

%max fr
if inputs.maxFRCheckBox.Value == 1
    maxfrFig = figure(1);
    maxfrFig.Position = [140 252.5 125 170];
    hold on
    scatter(ones(1,length(maxfr.(grpName{1}))),maxfr.(grpName{1}),40,lightgray,'filled')
    scatter(2.*ones(1,length(maxfr.(grpName{2}))),maxfr.(grpName{2}),40,lightgray,'filled')
    if hMaxFR == 0 %if data are normal, plot mean w sem
        errorbar(1,mean(maxfr.(grpName{1})),sem(maxfr.(grpName{1}),2),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,mean(maxfr.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,mean(maxfr.(grpName{2})),sem(maxfr.(grpName{2}),2),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,mean(maxfr.(grpName{2})),125,nrColor,'filled')
    elseif hMaxFR ~= 0 %if data non-normal, plot median w IQR
        [lowerError.maxfr.(grpName{1}),upperError.maxfr.(grpName{1})] = iqrError(maxfr.(grpName{1}),2);
        [lowerError.maxfr.(grpName{2}),upperError.maxfr.(grpName{2})] = iqrError(maxfr.(grpName{2}),2);
        errorbar(1,median(maxfr.(grpName{1})),lowerError.maxfr.(grpName{1}),upperError.maxfr.(grpName{1}),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,median(maxfr.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,median(maxfr.(grpName{2})),lowerError.maxfr.(grpName{2}),upperError.maxfr.(grpName{2}),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,median(maxfr.(grpName{2})),125,nrColor,'filled')
    end        
    ylabel('max firing rate (Hz)')
    maxfrAx=gca;
    setAx(maxfrAx);
    maxfrAx.XTickLabel={'resp', 'no resp'};
    xlim([0.5 2.5])
    allMaxFR = [maxfr.(grpName{1}) maxfr.(grpName{2})];
    ylim([0 1.2*max(allMaxFR)])
end

%halfwidth
if inputs.halfwidthCheckBox.Value == 1
    hwFig = figure(2);
    hwFig.Position = [140 5 125 170];
    hold on
    scatter(ones(1,length(hwidth.(grpName{1}))),hwidth.(grpName{1}),40,lightgray,'filled')
    scatter(2.*ones(1,length(hwidth.(grpName{2}))),hwidth.(grpName{2}),40,lightgray,'filled')
    if hHalfwidth == 0 %if data are normal, plot mean w sem
        errorbar(1,mean(hwidth.(grpName{1})),sem(hwidth.(grpName{1}),2),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,mean(hwidth.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,mean(hwidth.(grpName{2})),sem(hwidth.(grpName{2}),2),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,mean(hwidth.(grpName{2})),125,nrColor,'filled')
    elseif hHalfwidth ~= 0 %if data non-normal, plot median w IQR
        [lowerError.hwidth.(grpName{1}),upperError.hwidth.(grpName{1})] = iqrError(hwidth.(grpName{1}),2);
        [lowerError.hwidth.(grpName{2}),upperError.hwidth.(grpName{2})] = iqrError(hwidth.(grpName{2}),2);
        errorbar(1,median(hwidth.(grpName{1})),lowerError.hwidth.(grpName{1}),upperError.hwidth.(grpName{1}),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,median(hwidth.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,median(hwidth.(grpName{2})),lowerError.hwidth.(grpName{2}),upperError.hwidth.(grpName{2}),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,median(hwidth.(grpName{2})),125,nrColor,'filled')
    end        
    ylabel('AP halfwidth (ms)')
    hwAx=gca;
    setAx(hwAx);
    hwAx.XTickLabel={'resp', 'no resp'};
    xlim([0.5 2.5])
    allHwidth = [hwidth.(grpName{1}) hwidth.(grpName{2})];
    ylim([0 1.2*max(allHwidth)])
end

%sag
if inputs.sagCheckBox.Value == 1
    sagFig = figure(3);
    sagFig.Position = [280 252.5 125 170];
    hold on
    scatter(ones(1,length(sag.(grpName{1}))),sag.(grpName{1}),40,lightgray,'filled')
    scatter(2.*ones(1,length(sag.(grpName{2}))),sag.(grpName{2}),40,lightgray,'filled')
    if hSag == 0 %if data are normal, plot mean w sem
        errorbar(1,mean(sag.(grpName{1})),sem(sag.(grpName{1}),2),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,mean(sag.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,mean(sag.(grpName{2})),sem(sag.(grpName{2}),2),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,mean(sag.(grpName{2})),125,nrColor,'filled')
    elseif hSag ~= 0 %if data non-normal, plot median w IQR
        [lowerError.sag.(grpName{1}),upperError.sag.(grpName{1})] = iqrError(sag.(grpName{1}),2);
        [lowerError.sag.(grpName{2}),upperError.sag.(grpName{2})] = iqrError(sag.(grpName{2}),2);
        errorbar(1,median(sag.(grpName{1})),lowerError.sag.(grpName{1}),upperError.sag.(grpName{1}),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,median(sag.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,median(sag.(grpName{2})),lowerError.sag.(grpName{2}),upperError.sag.(grpName{2}),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,median(sag.(grpName{2})),125,nrColor,'filled')
    end        
    ylabel('hyperpolarization induced sag (%)')
    sagAx=gca;
    setAx(sagAx);
    sagAx.XTickLabel={'resp', 'no resp'};
    xlim([0.5 2.5])
    allSag = [sag.(grpName{1}) sag.(grpName{2})];
    ylim([0 1.2*max(allSag)])
end

%membrane resistance
if inputs.RmCheckBox.Value == 1
    rmFig = figure(4);
    rmFig.Position = [280 5 125 170];
    hold on
    scatter(ones(1,length(rm.(grpName{1}))),rm.(grpName{1}),40,lightgray,'filled')
    scatter(2.*ones(1,length(rm.(grpName{2}))),rm.(grpName{2}),40,lightgray,'filled')
    if hRm == 0 %if data are normal, plot mean w sem
        errorbar(1,mean(rm.(grpName{1})),sem(rm.(grpName{1}),2),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,mean(rm.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,mean(rm.(grpName{2})),sem(rm.(grpName{2}),2),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,mean(rm.(grpName{2})),125,nrColor,'filled')
    elseif hRm ~= 0 %if data non-normal, plot median w IQR
        [lowerError.rm.(grpName{1}),upperError.rm.(grpName{1})] = iqrError(rm.(grpName{1}),2);
        [lowerError.rm.(grpName{2}),upperError.rm.(grpName{2})] = iqrError(rm.(grpName{2}),2);
        errorbar(1,median(rm.(grpName{1})),lowerError.rm.(grpName{1}),upperError.rm.(grpName{1}),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,median(rm.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,median(rm.(grpName{2})),lowerError.rm.(grpName{2}),upperError.rm.(grpName{2}),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,median(rm.(grpName{2})),125,nrColor,'filled')
    end        
    ylabel('membrane resistance (M\Omega)')
    rmAx=gca;
    setAx(rmAx);
    rmAx.XTickLabel={'resp', 'no resp'};
    xlim([0.5 2.5])
    allRm = [rm.(grpName{1}) rm.(grpName{2})];
    ylim([0 1.2*max(allRm)])
end

%membrane decay constant
if inputs.mTauCheckBox.Value == 1
    tauFig = figure(5);
    tauFig.Position = [420 252.5 125 170];
    hold on
    scatter(ones(1,length(tau.(grpName{1}))),tau.(grpName{1}),40,lightgray,'filled')
    scatter(2.*ones(1,length(tau.(grpName{2}))),tau.(grpName{2}),40,lightgray,'filled')
    if hSag == 0 %if data are normal, plot mean w sem
        errorbar(1,mean(tau.(grpName{1})),sem(tau.(grpName{1}),2),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,mean(tau.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,mean(tau.(grpName{2})),sem(tau.(grpName{2}),2),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,mean(tau.(grpName{2})),125,nrColor,'filled')
    elseif hSag ~= 0 %if data non-normal, plot median w IQR
        [lowerError.tau.(grpName{1}),upperError.tau.(grpName{1})] = iqrError(tau.(grpName{1}),2);
        [lowerError.tau.(grpName{2}),upperError.tau.(grpName{2})] = iqrError(tau.(grpName{2}),2);
        errorbar(1,median(tau.(grpName{1})),lowerError.tau.(grpName{1}),upperError.tau.(grpName{1}),'color',inputs.plotColor,'linewidth',3,'CapSize',0)
        scatter(1,median(tau.(grpName{1})),125,inputs.plotColor,'filled')
        errorbar(2,median(tau.(grpName{2})),lowerError.tau.(grpName{2}),upperError.tau.(grpName{2}),'color',nrColor,'linewidth',3,'CapSize',0)
        scatter(2,median(tau.(grpName{2})),125,nrColor,'filled')
    end        
    ylabel('membrane decay \tau (ms)')
    tauAx=gca;
    setAx(tauAx);
    tauAx.XTickLabel={'resp', 'no resp'};
    xlim([0.5 2.5])
    allTau = [tau.(grpName{1}) tau.(grpName{2})];
    ylim([0 1.2*max(allTau)])
end
end