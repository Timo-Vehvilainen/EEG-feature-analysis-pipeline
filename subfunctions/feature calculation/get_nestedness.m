function [PLV] = get_nestedness(signal, highcomp, fs)
% CALCULATE NESTEDNESS
% input:   signal     = signal from one channel
%          highcomp   = frequency band for higher components: [3 8] [8-15] [15-30]
%          fs         = sampling frequency
% output:  PLV        = nestedness: PLV value between the low component and 
%                       filtered envelope of high component

lowcomp = [0.2 0.6]; 

% filter out 'low' (nesting) component 

[b,a] = butter(5,lowcomp(2)/(fs/2),'low'); 
s_low = filtfilt(b,a,signal'); 
[b,a] = butter(5,lowcomp(1)/(fs/2),'high');
s_low= filtfilt(b,a,s_low);

% filter out the and 'high' (nested) component
[b,a] = butter(7,highcomp(2)/(fs/2),'low'); 
s_high = filtfilt(b,a,signal'); 
[b,a] = butter(7,highcomp(1)/(fs/2),'high');
s_high= filtfilt(b,a,s_high);

% Take envelope for higher components and filter it with the same filter 
% that you used for lower component ([0.2 - 0.6] Hz).
[yupper,ylower] = envelope(s_high); 

[b,a] = butter(5,lowcomp(2)/fs,'low'); 
env_high= filtfilt(b,a,yupper); 
[b,a] = butter(5,lowcomp(1)/fs,'high');
env_high= filtfilt(b,a,env_high);

% Convert both: low component and filtered envelope of high component 
% into complex form (Hilbert transform).
h_high = hilbert(env_high); 
h_low = hilbert(s_low); 

% Compute PLV between them.
PLV = get_PLV(h_high, h_low); 

end


