function [v,ess] = bartlett_mv(Z,taper_method,M)

switch taper_method
  case 'none'
    taper = 0;
  case 'tukey'
    taper = 1;
  case 'parzen'
    taper = 2;
  case 'bartlett'
    taper = 3;
end

if nargin < 2
  taper = 0;
end

% Z is a NxMx2

N = size(Z,1); % Sample size
dim = size(Z,2); % Number of correlations
v = zeros(dim);

% Select cut-off for ACF/tapering
if nargin < 4
%   M = ceil(2*sqrt(N));
  M = N-2;
end

% Set up tapering
ks = abs([-M+1:0, 1:M-1]);

lambda_k = ones(2*M-1,1);
if taper == 1 % Tukey
  lambda_k(:) = 0.5 .* (1 + cos(pi.*ks./M) );
elseif taper == 2 % Parzen
  co = ceil(M/2);
  lambda_k(ks<=co) = 1-6.*(ks(ks<=co)./M).^2+6.*(ks(ks<=co)./M).^3;
  lambda_k(ks>co) = 2.*(1-ks(ks>co)./M).^3;
elseif taper == 3 % Bartlett
  lambda_k(:) = 1 - ks./M;
end

ac_lags = 1:N-1;
is_ac_tap = ismember(ac_lags,ks);
ac_lambda_k = lambda_k(M+1:end);

acfs = zeros(length(ac_lags),dim);
for i = 1:dim
  A = zscore(Z(:,i,1));
  r_ii = autocorr(A,M-1);
  
  r_ii = r_ii(2:end);
  r_ii(is_ac_tap) = r_ii(is_ac_tap) .* ac_lambda_k;
  r_ii(~is_ac_tap) = 0;
  
  acfs(:,i) = r_ii;
end

for i = 1:dim  
  vi = zeros(dim,1);
  acfi = acfs(:,i);
  for j = i+1:dim
    us = 1:N-1;
    vi(j) = (1 + 2*sum((N-us)'.*acfi.*acfs(:,j)) ./ N) ./ N;
  end
  v(:,i) = vi;
end

v = v + v';
ess = 1 + 1 ./ v;
ess(eye(dim)==1) = N;