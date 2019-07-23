function [rms] = getRMS(data)
%used to determine standard deviation of electrical noise for detecting synaptic events
%make sure data is baseline subtracted
rms = sqrt(mean(data));

end