function [h, p, t, df] = rsdt(Xm, Xs, Xc, Ym, Ys, Yc, n, r, alpha)

% Revised Standardized Difference Test
% Crawford JR, Garthwaite PH. Neuropsychology 2005 May;19(3):318-31)
% Two tailed test
%
% Xm, Xs, Xc: test X mean, std, single case value
% Ym, Ys, Yc: test Y mean, std, single case value
% n: sample size of normative sample
% r: correlation between X & Y in normative sample
% alpha: significance value
%
%
% Author: maarten.schrooten@uzleuven.be

% Z scores
Xz = ( Xc - Xm ) / Xs;
Yz = ( Yc - Ym ) / Ys;

% Degrees of freedom
df = n - 1;
df2 = df^2;
r2 = r^2;

% Compute Psi = formula 6
tails = 2;
% tinv from Alois Schloegl <a.schloegl@ieee.org> tinv.m $Revision: 1.1
y = abs( (sign(alpha/tails - 1/2).*sqrt(df./betainv(2*min(alpha/tails, 1-alpha/tails), df/2, 1/2) - df)) );
y2 = y^2;
psi = abs( Xz - Yz ) ...
      / ...
      sqrt( ( ( n + 1 ) / n ) * ...
            (   ( 2 - 2 * r ) + ...
              + ( 2 * (1 - r2) / (n - 1) ) ...
              + ( (5 + y2) * (1 - r2) / ( 2 * df2 ) ) ...
              + ( r * ( 1 + y2 ) * (1 - r2) / ( 2 * df2 ) ) ...
            ) ...
          );

% Hypothesis proven?
if psi > y, h = 1; else h = 0; end

% t score
a = (1 + r) * (1 - r2);
b = (1 - r) * ( 4*df2 + 4*(1+r)*df + (1+r)*(5+r) );
c = -2 * (Xz - Yz)^2 * ( n * df^2 / (n + 1) );
t = sqrt( ( -b + sqrt( b^2 - 4*a*c ) ) / ( 2 * a ) ); % formula 7

% p value (http://www.statsci.org/matlab/statbox.html tp.m)
p = (1 - 0.5 * ( 1 + betainc( t^2 / (df + t^2), 0.5, 0.5 * df ) ) ) * tails; % two-tailed

end
