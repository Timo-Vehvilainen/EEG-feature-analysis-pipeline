%Timo Vehviläinen August 2019

function [line_length] = get_line_length(epoch, varargin)
%GET_LINE_LENGTH Calculates the line length of a signal, or the running sum
%of the absolute difference between consecutive samples and normalizing it.
%
%PARAMETERS:
%       epoch:
%           this is the signal to be analysed, as a vector of time samples
%OPTIONAL PARAMETERS:
%       window_size_samples:
%           the size of the windows to be used for calculating line length.
%           Defaults to 250 samples.
%       overlap_samples:
%           the number of samples that consecutive windows are supposed to
%           overlap. Defaults to 30.
%   

%parse input
p = inputParser;
addParameter(p, 'window_size_samples', 250, @isnumeric);
addParameter(p, 'overlap_samples', 30, @isnumeric);
parse(p, varargin{:});
S = struct2cell(p.Results); 
[overlap_samples, window_size_samples] = deal(S{:});

%calculate intermediate variabels
window_interval = window_size_samples - overlap_samples;
array_size = floor((numel(epoch) - window_interval) / overlap_samples);

%initilize
line_length = zeros(array_size, 1);
j = 1;
i = 1;

%go through the epoch in window-sized chunks
while j < (numel(epoch) - window_size_samples)
    window = epoch(j:(j+window_size_samples));
    
    %for each window, compute the total sum of absolute differences between
    %consecutive samples
    L_i = 0;
    for z = 1:(numel(window)-1)
        L_i = nansum([L_i, abs(window(z) - window(z+1))]);
    end
    line_length(i) = L_i;
    
    %move to the next window
    i = i + 1;
    j = j + window_interval;
end
%normalize by dividing line lengths by the total sum of the line lengths
line_length = line_length./nansum(line_length);
end

