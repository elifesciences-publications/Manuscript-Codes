function [output] = plastFxn(inputs)
%for analysis of plasticity data
%modeled off ltpVC.m code
%created 04-20-2018, edited 04-21-18

%% INIT VARS
%sample rate
samplerate = inputs.SamplerateHzEditField.Value;

%stimulus start
stimstart = inputs.StimulusstartsecEditField.Value*samplerate;

%set baseline -501ms to -1ms before stim
blstart=stimstart-0.5*samplerate-1*(samplerate/1000); %same as getPSCData fxn
blend=blstart+0.5*samplerate;

%access resistance injection
rastart = inputs.RainjectionstartsecEditField.Value*samplerate;
raend = inputs.RainjectionendsecEditField.Value*samplerate;

%% LOAD DATA

cd(inputs.plastDir);

%stim information
stiminfo = textread('stimdur.txt'); %#ok<DTXTRD>
stimdur = stiminfo(1)/(10^6); %gives stimuulation duration in seconds
stimend = stimstart+stimdur*samplerate;

%load Plast data
%pre induction
cd([inputs.plastDir, '/pre'])
contents = dir('*.abf');
filenames = {contents.name}';
preFiles = fullfile(cd,filenames);
clear filenames
for ii = 1:length(preFiles)
    thisData = abfload(preFiles{ii},'sweeps','a');
    filtPlastData.pre(ii,:,:) = sgolayfilt(thisData(:,1,:),3,11,[],1);
    raData.pre(:,:,:,ii) = thisData;
    clear thisData
end
noSweeps.pre = size(filtPlastData.pre,3);

%post induction
cd([inputs.plastDir, '/post'])
contents = dir('*.abf');
filenames = {contents.name}';
postFiles = fullfile(cd,filenames);
clear filenames
for ii = 1:length(postFiles)
    thisData = abfload(postFiles{ii},'sweeps','a');
    filtPlastData.post(ii,:,:) = sgolayfilt(thisData(:,1,:),3,11,[],1);
    raData.post(:,:,:,ii) = thisData;
    clear thisData
end
noSweeps.post = size(filtPlastData.post,3);

%% CHECK Ra
checkRa = cat(4,raData.pre,raData.post);
[Ra, meanRa, keep] = getRa(checkRa,rastart,raend,samplerate)

if keep == 1
    %% ANALYZE DATA

    %pre induction EPSC amplitude
    %time vector for pre-induction
    t.pre = -1*inputs.baseDurationminEditField.Value:1/inputs.StimperminEditField.Value:-1/inputs.StimperminEditField.Value;
    
    %EPSC amplitude
    bsData.pre = -1.*ones(size(filtPlastData.pre,2),length(preFiles)*noSweeps.pre); %get baseline subtracted data
    count = 0;
    for ii = 1:length(preFiles)
        for jj = 1:noSweeps.pre
            count = count +1;
            bsData.pre(:,count) = filtPlastData.pre(ii,:,jj) - mean(filtPlastData.pre(ii,blstart:blend,jj));
        end
    end
    [output.amp.pre, lat.pre, rise.pre, decay.pre, tPeaks.pre] = getPSCData(bsData.pre,stimstart,stimdur,samplerate,-1);
    baseAmp = mean(output.amp.pre(~isnan(output.amp.pre)));
    output.normAmp.pre = 100.*(output.amp.pre./baseAmp);
    
    %post induction EPSP slope, Vm, Rm
    %time vector for pre-induction
    t.post = 75/60:1/inputs.StimperminEditField.Value:inputs.postDurationminEditField.Value+75/60;
    t.post(end)=[];
    
    %EPSC amplitude
    bsData.post = -1.*ones(size(filtPlastData.post,2),length(postFiles)*noSweeps.post); %get baseline subtracted data
    count = 0;
    for ii = 1:length(postFiles)
        for jj = 1:noSweeps.post
            count = count +1;
            bsData.post(:,count) = filtPlastData.post(ii,:,jj) - mean(filtPlastData.post(ii,blstart:blend,jj));
        end
    end
    [output.amp.post, lat.post, rise.post, decay.post, tPeaks.post] = getPSCData(bsData.post,stimstart,stimdur,samplerate,-1);
    output.normAmp.post = 100.*(output.amp.post./baseAmp);
    
    %Validate steady response
    output.lastFive = output.amp.post(end-(inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value-1):end);
    output.lastFive(isnan(output.lastFive)) = [];
    midFiveStart = inputs.StimperminEditField.Value*(inputs.postDurationminEditField.Value/2)-(inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value)/2+1;
    midFiveEnd = midFiveStart+inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value-1; 
    midFive = output.amp.post(midFiveStart:midFiveEnd);
    midFive(isnan(midFive)) = [];
    [PlF,xlF] = ecdf(output.lastFive);
    [PmF,xmF] = ecdf(midFive);
    h_last = adtest(output.lastFive); h_mid = adtest(midFive);
    if sum([h_last h_mid]) == 0
        [~,p_ss_tt] = ttest2(output.lastFive,midFive) %#ok<*NOPTS>
        if p_ss_tt > 0.05
            ssKeep = 1
        else
            ssKeep = 0
        end
    else
        p_ss_mwu = ranksum(output.lastFive,midFive)
        if p_ss_mwu > 0.05
            ssKeep = 1
        else
            ssKeep = 0
        end
    end
    
    %% DISPLAY DATA
    %EPSC amplitude
    plasticityFig = figure(1);
    plasticityFig.Position = [675 635 560 135];
    hold on
    line([-1+t.pre(1) t.post(end)],[baseAmp baseAmp],'color','k','linewidth',1.5,'linestyle','--')
    scatter(t.pre,output.amp.pre,50,'k','filled')
    scatter(t.post,output.amp.post,50,inputs.plotColor,'filled')
    lastFiveNorm = output.normAmp.post(end-inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value-1:end);
    lastFiveNorm(isnan(lastFiveNorm)) = [];
    text(15,2.*baseAmp,['Last 5 min.: ',num2str(mean(lastFiveNorm)),'%'],'FontWeight','bold')
    xlabel('t (min)')
    ylabel('EPSC amplitude (pA)')
    xlim([-6 31])
    pAx = gca;
    pAx.Box = 'off'; pAx.TickDir = 'out'; pAx.YColor = 'k'; pAx.XColor = 'k'; pAx.LineWidth = 1;
    
    %waveforms
    %get mean waveforms
    [output.alignedTraces.pre,output.mTrace.pre,tMaxSlope.pre, tMaxSlopeMeanTrace.pre] = meanTraceMaxRise(bsData.pre(:,~isnan(output.amp.pre)),tPeaks.pre(~isnan(output.amp.pre)),samplerate,.15,-1); %requires data in t x trial format
    thesePostData = bsData.post(:,end-inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value-1:end);
    thesePeaksPost = tPeaks.post(end-inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value-1:end);
    successSweepsLastFive = ~isnan(output.amp.post(end-inputs.StimperminEditField.Value*inputs.AnalysisdurationminEditField.Value-1:end));
    [output.alignedTraces.post,output.mTrace.post,tMaxSlope.post, tMaxSlopeMeanTrace.post] = meanTraceMaxRise(thesePostData(:,successSweepsLastFive),thesePeaksPost(successSweepsLastFive),samplerate,.15,-1); %requires data in t x trial format
    
    %plot
    traceFig = figure(2);
    traceFig.Position = [675 385 285 165];
    hold on
    plot((1/samplerate:1/samplerate:size(output.mTrace.pre,1)/samplerate).*1000,output.mTrace.pre,'color','k','linewidth',1.5)
    plot((1/samplerate:1/samplerate:size(output.mTrace.post,1)/samplerate).*1000,output.mTrace.post,'color',inputs.plotColor,'linewidth',1.5)
    xlabel('t (ms)')
    ylabel('current (pA)')
    xlim([1000*(1/samplerate) 1000*(length(output.mTrace.pre)/samplerate)])
    ylim([-1*max(output.lastFive) .15*max(output.lastFive)])
    traceAx = gca;
    traceAx.Box = 'off'; traceAx.TickDir = 'out'; traceAx.YColor = 'k'; traceAx.XColor = 'k';
    
    %steady state
    midColor = inputs.plotColor+ (.9-max(inputs.plotColor));
    ssFig = figure(3);
    ssFig.Position = [965 385 270 165];
    title('post induction steady amplitude check')
    hold on
    plot(xmF,PmF,'linewidth',1.5,'color',midColor)
    plot(xlF,PlF,'linewidth',1.5,'color',inputs.plotColor)
    xlabel('EPSC amp (pA)')
    ylabel('P')
    xlim([0 1.5*max([midFive; output.lastFive])])
    legend('middle 5min','last 5min','location','northwest')
    axes('Position',[.75 .33 .15 .5])
    hold on
    scatter(1.*ones(length(midFive),1),midFive,50,midColor,'filled')
    scatter(2.*ones(length(output.lastFive),1),output.lastFive,50,inputs.plotColor,'filled')
    ylabel('EPSC amp (pA)')
    ylim([0 max([midFive; output.lastFive])+20])
    xlim([.5 2.5])
    inAx = gca;
    inAx.XTick = [];
    
end
end