function [output] = rampinjFxn(inputs)
%for analysis of biophysics differences bw responding and nonresponding neurons in minstim data
%created 02-11-19 %modified 02-11-18
close all

%% INIT VARS
%sample rate
samplerate=inputs.SamplerateHzEditField.Value;

%stimulus start in ms
injstart = inputs.RampstartmsEditField.Value*(samplerate/1000);

%baseline
blstart=injstart-inputs.BaselinedurationmsEditField_2.Value*(samplerate/1000)-1*(samplerate/1000);
blend=injstart-1*(samplerate/1000);

%% LOAD DATA
cd(inputs.ciDir);

%data
contents = dir('*.abf');
filenames = {contents.name}';
ciFiles = fullfile(cd,filenames);

output.rampData = abfload(ciFiles{1},'sweeps','a');
vTrace(:,:) = output.rampData(:,1,:);
iTrace(:,:) = output.rampData(:,2,:);
t = (1/samplerate:1/samplerate:size(output.rampData,1)/samplerate)';

%% ANALYZE DATA
noSweeps = size(output.rampData,3); %gives number of sweeps per file

%init vars
dVmRamp = diff(vTrace); %takes approx derivative (lose element n = 1)
dVmRamp = dVmRamp.*(samplerate/1000); %converts to V/s
vTrace(1,:) = []; %remove element n = 1 from Vm data
t(1,:) = []; iTrace(1,:) = [];
threshold = 20*ones(1,noSweeps); %manually set to 20V/s
apthresh = -1.*ones(250,2,noSweeps);

%find AP thresholds for each spike
for ii = 1:noSweeps
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
        plot(t,vTrace(:,ii))
        scatter(tThresh/samplerate,vTrace(tThresh,ii))
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
        if ii ~= noSweeps
            clf
        else
            close(threshFig)
        end
        
        %put into matrix
        apthresh(1:noSpkSweep,1,ii) = vTrace(tThresh,ii);
        apthresh(1:noSpkSweep,2,ii) = tThresh;
    end
    
    %get rheobase
    output.rheobase(ii) = iTrace(apthresh(1,2,ii),ii); %in pA
    
    %get threshold
    output.spkThreshold(ii) = apthresh(1,1,ii); %in mV, use over threshold from square pulses
    
    %get Vrest
    output.vRest(ii) = mean(output.rampData(blstart:blend,1,ii)); %in mV
end

%plot
hiSweep = randperm(noSweeps,1);
rheoFig = figure(1);
rheoFig.Position = [570 285 775 450];
subplot(4,1,1:3)
hold on
for ii = 1:noSweeps
    if ii ~= hiSweep
        if ii == noSweeps
            plot(t,vTrace(:,noSweeps),'linewidth',1,'color',[.67 .67 .67]);
            plot(t,vTrace(:,hiSweep),'linewidth',2,'color',inputs.plotColor)
        else
            plot(t,vTrace(:,ii),'linewidth',1,'color',[.67 .67 .67]);
        end
    end
end
setAx(gca)
xlabel('time (s)')
ylabel('membrame voltage (mV)')
subplot(4,1,4)
plot(t,iTrace(:,hiSweep),'linewidth',2,'color',inputs.plotColor)
setAx(gca)
xlabel('time (s)')
ylabel('injected current (pA)')

output %#ok<*NOPRT>
end
