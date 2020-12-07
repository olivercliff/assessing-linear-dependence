dim_X = 1;
dim_Y = 1;
dim_W = 0;

seed = 'shuffle';

% Dimension of Z
D = dim_X + dim_Y + dim_W;

% Partitions
p_X = 1:dim_X;
p_Y = dim_X+1:dim_X+dim_Y;
p_W = dim_X+dim_Y+1:D;

% Seed the RNG
rng(seed);

% Generate random samples
Z = randn(D,T);

X = Z(p_X,:)';
Y = Z(p_Y,:)';
W = Z(p_W,:)';

[rho_MATLAB, pval_MATLAB] = partialcorr(X,Y,W);
[rho, pval] = pcorr(X,Y,W,'test','asymptotic');
[rho1, pval1] = pcorr(X,Y,W,'test','modified','varianceEstimator','none');

  %% Print results

fprintf('---\n');
fprintf('Correlation computed from MATLAB in-built: %.5g (p = %.5g)\n', rho_MATLAB, pval_MATLAB);
fprintf('Correlated computed from this toolkit: %.5g (p = %.5g)\n', rho, pval);
fprintf('Correlated computed from this toolkit (option 2): %.5g (p = %.5g)\n', rho1, pval1);