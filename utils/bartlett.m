function [v,ess] = bartlett(Z,taper_method,full,M)

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

if nargin < 3
  full = false;
  if nargin < 2
    taper = 0;
  end
end

% Z is a NxMx2

N = size(Z,1); % Sample size
dim = size(Z,2); % Number of correlations
assert(size(Z,3) == 2);
v = zeros(dim);

% Select cut-off for ACF/tapering
if nargin < 4
%   M = ceil(2*sqrt(N));
  M = N-2;
end

% Set up tapering
ks = abs([-M+1:0, 1:M-1]);

lambda_k = ones(2*M-1,1);
biasing = 'unbiased';
if taper == 1 % Tukey
  lambda_k(:) = 0.5 .* (1 + cos(pi.*ks./M) );
elseif taper == 2 % Parzen
  co = ceil(M/2);
  lambda_k(ks<=co) = 1-6.*(ks(ks<=co)./M).^2+6.*(ks(ks<=co)./M).^3;
  lambda_k(ks>co) = 2.*(1-ks(ks>co)./M).^3;
elseif taper == 3 % Bartlett
  lambda_k(:) = 1 - ks./M;
end

if taper > 0
  biasing = 'biased';
end

lags = [N-1:-1:0 1:N-1];
is_tap = ismember(abs(lags),ks);

ac_lags = 1:N-1;
is_ac_tap = ismember(ac_lags,ks);
ac_lambda_k = lambda_k(M+1:end);

a = 1; b = 2;
c = 3; d = 4;
sig = zeros(4,1);
xc = zeros(2*N-1,4,4);
for i = 1:dim
  A = zscore(Z(:,i,1));
  B = zscore(Z(:,i,2));
  
  if full
    sig(a) = std(A);
    sig(b) = std(B);
  
    xc(:,a,b) = xcorr(A,B,biasing)./sig(a)./sig(b);
    for j = 1:dim
        C = Z(:,j,1);
        D = Z(:,j,2);
        sig(c) = std(C);
        sig(d) = std(D);

        xc(:,a,c) = xcorr(A,C,biasing)./sig(a)./sig(c);
        xc(:,a,d) = xcorr(A,D,biasing)./sig(a)./sig(d);

        xc(:,b,d) = xcorr(B,D,biasing)./sig(b)./sig(d);
        xc(:,b,c) = xcorr(B,C,biasing)./sig(b)./sig(c);
        xc(:,c,d) = xcorr(C,D,biasing)./sig(c)./sig(d);

        for a_i = 1:4
          for c_i = 1:4
            xc(is_tap,a_i,c_i) = xc(is_tap,a_i,c_i) .* lambda_k;
            xc(~is_tap,a_i,c_i) = 0;
          end
        end

        rho_ab = xc(N,a,b);
        rho_cd = xc(N,c,d);

        v(i,j) = sum( xc(:,a,c).*xc(:,b,d) + xc(:,a,d).*xc(:,b,c) ...
                  - rho_ab .* ( xc(:,a,c).*xc(:,a,d) + xc(:,b,c).*xc(:,b,d) ) ...
                  - rho_cd .* ( xc(:,a,c).*xc(:,b,c) + xc(:,a,d).*xc(:,b,d) ) ...
                  + 0.5 .* rho_ab .* rho_cd .* ( xc(:,a,c).^2 + xc(:,b,c).^2 + xc(:,a,d).^2 + xc(:,b,d).^2 ) );
    end
  else
    r_ii = autocorr(A,M-1); r_ii = r_ii(2:end);
    r_jj = autocorr(B,M-1); r_jj = r_jj(2:end);
    
    % Taper the results
    r_ii(is_ac_tap) = r_ii(is_ac_tap) .* ac_lambda_k;
    r_ii(~is_ac_tap) = 0;
    r_jj(is_ac_tap) = r_jj(is_ac_tap) .* ac_lambda_k;
    r_jj(~is_ac_tap) = 0;
    
    u = 1:N-1;
    v(i,i) = (1 + 2*sum((N-u)'.*r_ii.*r_jj) ./ N) ./ N;
  end
end

ess = 1 + 1 ./ v;