function [v,V] = bartlett_correction_mv(Z,taper,full)

if nargin < 3
  full = false;
  if nargin < 2
    taper = 0;
  end
end

% Z is a 2xMxN

N = size(Z,1); % Sample size
dim = size(Z,2); % Number of correlations
assert(size(Z,3) == 2);
v = zeros(dim);

% Set up tapering
M = ceil(2*sqrt(N));
ks = abs([-M+1:0, 1:M-1]);
ac_ks = 1:M-1;

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
  A = Z(:,i,1);
  B = Z(:,i,2);
  
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
    acA = autocorr(A,N-1); acA = acA(2:end);
    acB = autocorr(B,N-1); acB = acB(2:end);
    
    % Taper the results
    acA(is_ac_tap) = acA(is_ac_tap) .* ac_lambda_k;
    acA(~is_ac_tap) = 0;
    acB(is_ac_tap) = acB(is_ac_tap) .* ac_lambda_k;
    acB(~is_ac_tap) = 0;
    
    v(i,i) =  1 + 2*sum( (N-1:-1:1)' .*acA.*acB )/N;
  end
end

V = v / N;