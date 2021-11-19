function [h, p, t, df] = ttestch(Xm, Xs, Xc, n, alpha)

% Crawford & Howell's modified t-test
% Crawford JR, Garthwaite PH. Neuropsychology 2005 May;19(3):318-31)
% One tailed test
%
% Xm, Xs, Xc: test X mean, std, single case value
% n: sample size of normative sample
% alpha: significance value
%
%
% Author: maarten.schrooten@uzleuven.be

% Formula 1
t = ( Xc - Xm ) / ( Xs * sqrt( (n + 1) / n ) );

% p value (http://www.statsci.org/matlab/statbox.html tp.m)
df = n - 1;
tails = 1;
p = (1 - 0.5 * ( 1 + betainc( t^2 / (df + t^2), 0.5, 0.5 * df ) ) ) * tails;

% Hypothesis rejected?
if p < alpha, h = 1; else h = 0; end

end
