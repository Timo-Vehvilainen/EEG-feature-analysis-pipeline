
function PLV = get_PLV(X, Y)
% CALCULATE PLV
% input:   X   = complex form signal from first channel
%          Y   = complex form signal from second channel
% output:  PLV = PLV value between components

% NOTE:    signals must be in complex format (Hilbert) 

    X = X(:);
    Y = Y(:);
    PLV = abs(sum(exp(1i * (angle(X .* conj(Y))))) ./ length(X));

end % end
