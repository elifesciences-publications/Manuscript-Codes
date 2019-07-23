function [waveforms] = findAPwaveforms(data,centerTimes,sweep,sample_rate)
%gets trace for each AP waveform for plotting purposes

waveforms=zeros(length(centerTimes),61);
for ii = 1:length(centerTimes)
    waveforms(ii,:)=data(centerTimes(ii)-0.002*sample_rate:centerTimes(ii)+0.004*sample_rate,1,sweep);
end
end