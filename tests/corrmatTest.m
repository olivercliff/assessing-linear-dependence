addpath(genpath('..'))

% Generate random samples
T = 100;
D = 3;
X = randn(T,D);

[rho_MATLAB,pval_MATLAB] = corr(X);

[rho,pval] = corrb(X,'test','finite');

fprintf('Correlation computed from MATLAB in-built: %s (p = %s)\n', mat2str(rho_MATLAB), mat2str(pval_MATLAB));
fprintf('Correlated computed from this toolkit: %s (p = %s)\n', mat2str(rho), mat2str(pval));