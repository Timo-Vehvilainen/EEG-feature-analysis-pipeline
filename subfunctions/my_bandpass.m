function [filtered_signal] = my_bandpass(signal, thresholds, sampling_rate)
%MY_BANDPASS A simple wrapper for a basic band-passing function with some
%preprocessing such as mirroring. The bandpass filter is run as two pairs
%of forward-reverse filters, one highpass filter (5th order Butterworth) and 
%one lowpass filter (7th order Butterworth)
%   PARAMETERS:
%               signal:
%                       a vector of time_samples from a signal
%               thresholds:
%                       a two-element vector of the form:
%                           [highpass_threshold, lowpass_threshold] (in Hz)
%               sampling_rate:
%                       the sampling rate for the signal
%OUTPUT is the filtered version of the original signal

L = length(signal);
% Add edges before filtering and Hilbert (to fix edge effects)
  s1_big = [flipud(signal); signal; flipud(signal)];
% s1_big = signal;

% HPF
  [b, a] = butter(5, thresholds(1)/(sampling_rate/2), 'high');
   s1_big_f1 = filtfilt(b, a, s1_big);
   
% LPF
  [b, a] = butter(7, thresholds(2)/(sampling_rate/2), 'low');  
   s1_big_f1_f2 = filtfilt(b, a, s1_big_f1);
  
% cut edges
  filtered_signal = s1_big_f1_f2(L+1:L+L, :);
% filtered_signal = s1_big_f1_f2;
end

