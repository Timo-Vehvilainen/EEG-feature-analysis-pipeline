function wpli = get_wPLI(X, Y, debias, freq_band, fs)
% WEIGHTED PLI
% input:   X        = complex form signal from first channel
%          Y        = complex form signal from second channel
%          debias   = [0]: biased, [1]: debiased
% output:  wpli     = wPLI value between the signals

% NOTE: there are 2 variations of wPLI: biased and debiased. 
% The last one is better. So, put '1' for input 'debias' variable.

    L = length(X);

  % cross-spectral density
    Pxy = cpsd(X, Y, L, 1, [], []);

  % compute wpli
    Pxy = imag(Pxy);         % make everything imaginary
    outsum = nansum(Pxy, 1); % compute the sum; this is 1 x size(2:end)
    outsumW = nansum(abs(Pxy), 1); % normalization of the WPLI

    if debias == 1
       outssq = nansum(Pxy .^ 2, 1);
       wpli = (outsum .^ 2 - outssq) ./ (outsumW .^ 2 - outssq); % do the pairwise

    else
       wpli = outsum ./ outsumW; % estimator of E(Im(X))/E(|Im(X)|)

    end

end 

