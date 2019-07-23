%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Action Potential Threshold %%
%%%%%%%%%%% for ramp protocol %%%%%%%%%%%%
%%%%%%%%%% Created: 07-08-2017 %%%%%%%%%%%
%%%%%%%%%%% Edited: 07-08-2017 %%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [apthresh]=findAP_ramp(trace,baseVm,samplerate)
no_sweeps = size(trace,2); %gives number of sweeps per file

%init vars
thresh = 25; %can also set to be a function of the sd dVm
dVmRamp = diff(trace); %takes approx derivative (lose element n = 1)
dVmRamp = dVmRamp.*(samplerate/1000); %converts to V/s
trace(1,:) = []; %remove element n = 1 from Vm data
t = samplerate^-1:samplerate^-1:size(trace,1)/samplerate;
dbaseVm = diff(baseVm); %get approx derivative of baseline V-mem to set threshold
dbaseVm = dbaseVm.*(samplerate/1000); %converts to V/s
sddVm = std(dbaseVm); %get standard deviation of the derivative of baseline V-mem
threshold = 20*ones(1,no_sweeps); %manually set to 20V/s
apthresh = -1.*ones(250,2,no_sweeps);

%find AP thresholds for each spike
for ii = 1:no_sweeps
    %inits
    oThresh = [];
    tThresh = [];
    
    %find AP thresholds
    oThresh = find(dVmRamp(:,ii)>threshold(ii)); %find points above the threshold of dVm during ramp
    if isempty(oThresh) ~= 1
        %inits per sweep
        noSpkSweep = 1;
        for jj = 1:length(oThresh)
            if jj == 1
                %get point before and after threshold
                p0 = oThresh(jj)-1; %after
                p1 = oThresh(jj)-2; %before
                afterVal = abs(20-dVmRamp(p0,ii));
                beforeVal = abs(20-dVmRamp(p1,ii));
                if afterVal <= beforeVal %if value after threshold is closer to threshold...
                    tThresh(noSpkSweep) = p0;
                else %if value before threshold is closer to threshold...
                    tThresh(noSpkSweep) = p1;
                end
                noSpkSweep = noSpkSweep+1;
            else
                if oThresh(jj) > .001*samplerate + tThresh(noSpkSweep-1) %if time over threshold is at least 1ms past the last AP (refractory period)
                    %get point before and after threshold
                    p0 = oThresh(jj)-1; %after
                    p1 = oThresh(jj)-2; %before
                    afterVal = abs(20-dVmRamp(p0,ii));
                    beforeVal = abs(20-dVmRamp(p1,ii));
                    if afterVal <= beforeVal %if value after threshold is closer to threshold...
                        tThresh(noSpkSweep) = p0;
                    else %if value before threshold is closer to threshold...
                        tThresh(noSpkSweep) = p1;
                    end
                    noSpkSweep = noSpkSweep +1;
                end
            end
        end
        noSpkSweep = noSpkSweep -1
        
        %Plots AP traces and thresholding
        threshFig = figure(99);
        threshFig.Position = [700 425 475 275];
        subplot(1,2,1)
        hold on
        plot(t,trace(:,ii))
        scatter(tThresh/samplerate,trace(tThresh,ii))
        ylabel('Membrane Potential (mV)')
        xlabel('Time (sec)')
        tAx = gca;
        tAx.Box ='off'; tAx.YColor = 'k'; tAx.XColor = 'k'; tAx.TickDir = 'out'; tAx.LineWidth = 1;
        subplot(1,2,2)
        hold on
        plot(t,dVmRamp(:,ii))
        line([t(1) t(end)],[threshold(ii) threshold(ii)],'color','k','linewidth',1.25);
        ylabel('dV/dt')
        xlabel('Time (sec)')
        tAx = gca;
        tAx.Box ='off'; tAx.YColor = 'k'; tAx.XColor = 'k'; tAx.TickDir = 'out'; tAx.LineWidth = 1;
        if ii ~= no_sweeps
            clf
        else
            close(threshFig)
        end
        
        %put into matrix
        apthresh(1:noSpkSweep,1,ii) = trace(tThresh,ii);
        apthresh(1:noSpkSweep,2,ii) = tThresh;
    end
end