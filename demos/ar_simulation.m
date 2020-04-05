clear
close all

if ~exist('granger_causality.m','file')
  addpath('..');
end

%% Reproduce a figure from the paper? Else choose params

figure_opts = {'2a','2b';
               '3a','3b';
               '4a','4b';
               '5a','5b';
               '6a','6b';
               '7a','7b';
               '8a','8b'};

which_figure = '4a';

fig_id = strcmp(figure_opts,which_figure);

verbose = true;

if any(fig_id(:))
  [fig,subfig] = find(fig_id);
  config = get_configuration(fig,subfig);
else
  config.T = 2^9; % Dataset length
  config.R = 100; % Number of runs/trials
  config.S = 1000; % Number of samples for MC distribution
  config.dim_X = 1; % set X dimension
  config.dim_Y = 1; % set Y dimension
  config.dim_W = 0; % set W dimension

  config.to_filter = 2; % None=0, FIR=1, IIR=2
  config.filter_order = 8; % 8th order
  config.ar = true; % Include the AR parameters (or use IID, false)

  config.causal = false; % Non-causal (null model, should we be testing true negatives or false positives?)

  config.is_granger = true; % Granger causality (not CMI)
  config.embedding = [-1 -1]; % Optimal (set 0 for p = q = 1, and >0 for a specific embedding on both)

  config.alpha = 0.05; % significance level

  config.seed = now; % RNG seed
  
  if ~config.is_granger && any(config.embedding > 0)
    warning('Embedding level %d set but variables will not be embeded (computing CMI, not GC).\n', config.embedding);
  end
end

config.R = 100;


%% Set up filters and simulator

if config.is_granger
  compute_measure = @(X,Y,W) granger_causality(X,Y,W,config.embedding,'none',false);
else
  compute_measure = @(X,Y,W) conditional_mutual_information(X,Y,W,'none',false);
end

if config.to_filter == 1
  % FIR filter
  a_coeff = fir1(config.filter_order, 0.5);
  b_coeff = 1;
elseif config.to_filter == 2
  % Butterworth (IIR) filter
  [a_coeff, b_coeff] = butter(config.filter_order, 0.5);
end

% Is the original signal autoregressive or spectrally white?
if config.ar
  phi_X = 0.3 .* eye(config.dim_X);
  phi_Y = -0.8 .* eye(config.dim_Y);
  phi_W = 0.4 .* eye(config.dim_W);
else
  phi_X = zeros(config.dim_X);
  phi_Y = zeros(config.dim_Y);
  phi_W = zeros(config.dim_W);
end

% Should we include a causal influence from Y to X?
if config.causal
  phi_XY = 0.03 .* eye(config.dim_X,config.dim_Y);
else
  phi_XY = zeros(config.dim_X,config.dim_Y);
end

% Autoregression parameters
phi_XW = zeros(config.dim_X, config.dim_W);

phi_YX = zeros(config.dim_Y, config.dim_X);
phi_YW = zeros(config.dim_Y, config.dim_W);

phi_WX = zeros(config.dim_W, config.dim_X);
phi_WY = zeros(config.dim_W, config.dim_Y);

Phi = [phi_X, phi_XY, phi_XW;
       phi_YX, phi_Y, phi_YW;
       phi_WX, phi_WY, phi_W];

% Innovation covariance
M = config.dim_X + config.dim_Y + config.dim_W;
Sigma = eye(M);

% Partitions for X,Y,W
p_X = 1:config.dim_X;
p_Y = config.dim_X+1:config.dim_X+config.dim_Y;
p_W = config.dim_X+config.dim_Y+1:M;
  
%% Run simulations

% Pre-allocate measure (GC or MI value)..
measure = zeros(config.R,1);

% ..and p-values
pvals_LR = zeros(config.R,1); % LR test
pvals_E = zeros(config.R,1); % exact test

rng(config.seed);

fprintf('Running simulations...');
if ~verbose
  fprintf(' (turn on verbose to see progress)');
else
  fprintf(' (turn off verbose to suppress output)');
end
fprintf('\n');

% Run sims
for r = 1:config.R
  Z = Sigma*randn(M,config.T);
  for t = 2:config.T
    Z(:,t) = Phi*Z(:,t-1) + Z(:,t);
  end

  % Partition the dataset (Z) into multiple time series (X,Y, and conditional W)
  X = Z(p_X,:)';
  Y = Z(p_Y,:)';
  W = Z(p_W,:)';

  % Filter the data to induce higher autocorrelation (if opted)
  if config.to_filter > 0
    X = filter(a_coeff,b_coeff,X,[],1);
    Y = filter(a_coeff,b_coeff,Y,[],1);
    W = filter(a_coeff,b_coeff,W,[],1);
  end

  % Compute measure (GC or MI)
  [measure(r),stats] = compute_measure(X,Y,W);
  
  % Generate p-values
  pvals_LR(r) = compute_significance(measure(r),stats,'lr');
  pvals_E(r) = compute_significance(measure(r),stats,'exact');
  
  if verbose
    if mod(r,10) == 0
      fprintf('Completed run %d/%d.\n', r, config.R);
    end
  end
end
  
%% Plot results
col_LR = [1 0 0];
col_E = [0 0 0];

figure;
hold on;
plot([0 1], [0 1], 'k--');
ph1 = plot(sort(pvals_LR),linspace(0,1,config.R), '-', 'color', col_LR, 'linewidth', 1);
ph2 = plot(sort(pvals_E),linspace(0,1,config.R), '-', 'color', col_E, 'linewidth', 1);
legend([ph1 ph2], 'LR test', 'Exact test','location', 'best');

fprintf('LR test FPR at %d%% significance: %.3g\n', config.alpha*100, mean(pvals_LR <= config.alpha) );
fprintf('Exact test FPR at %d%% significance: %.3g\n', config.alpha*100, mean(pvals_E <= config.alpha) );