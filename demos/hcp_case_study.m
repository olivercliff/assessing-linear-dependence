
clear
% close all

if ~exist('granger_causality.m','file')
  addpath('..');
end

%% Load HCP data

fprintf('Loading data...\n');
load data/hcp_rsfMRI.mat
fprintf('Done.\n');

% Normalise (otherwise we have global effects in trends, etc.)
fprintf('Normalising and detrending...\n');
for i = 1:size(dat,3)
  s_bold = dat(:,:,i)';
  s_bold = detrend(s_bold);
  s_bold = zscore(s_bold);
  dat(:,:,i) = s_bold';
end

%% Reproduce which HCP experiment from the paper?

% 1 (MI), 2 (MI with MV dim 2), 3 (GC, optimal embedding), 4 (GC, high
% embedding), 5 (GC with MV dim 2)
which_test = 4;

switch which_test
  case 1
    is_granger = false; % MI
    dims = 1; % 1D
  case 2
    is_granger = false; % MI
    dims = 2; % 2D
  case 3
    is_granger = true; % GC
    embedding = [-1 -1]; % auto-embedding
    dims = 1; % 1D
  case 4
    is_granger = true; % GC
    embedding = [100 100]; % Over optimistic (pesimistic?) embedding
    dims = 1; % 1D
  case 5
    is_granger = true; % GC
    embedding = [-1 1]; % auto-embedding
    dims = 2; % 1D
end

%% Default params

seed = 1; % Used for the paper
generate_tikz = true;

R = 1000; % Number of trials
S = 5000; % MC sample size
alpha = 0.05; % Significance level
seq = 201:1000; % Sequence of the data to take (cut off first+last 200)

to_filter = 2; % Filter the BOLD data?

% Bias towards only selecting the highest ACF in the subject (increases the FPR slightly)
only_high_ac_regions = false;

verbose = true;

%% Prelim stuff

fs = 0.72; % Sample time (s)

if to_filter == 1
%   fc = 0.2; % Low-pass (Hz)
  fc = [0.01 0.08]; % Passband (Hz)
  fo = 8; % Order
  Wn = fc./(fs/2);
  
  a_coeff = fir1(fo,Wn);
  b_coeff = 1;
elseif to_filter == 2
  fc = [0.01 0.08]; % Passband (Hz)
  fo = 4; % Order
  Wn = fc./(fs/2);
  
  [a_coeff, b_coeff] = butter(fo,Wn);
else
  warning('No filtering selected');
end

D = size(dat,1);
T = size(dat,2);
M = size(dat,3);

%% Run experiments

% Pre-allocate measure (GC or MI value)..
measure = zeros(R,1);

% ..and p-values
pvals_LR = zeros(R,1); % log-likelihood ratio test p-values
pvals_E = zeros(R,1); % exact test p-values

rng(seed);

% Run experiments
fprintf('Running experiments...\n');
for r = 1:R
  
  % Randomly sample subjects (without replacement) and..
  ss = zeros(dims,2);  
  ss(:) = randsample(M,dims*2);
  
  % .. regions (without replacement)
  ds = zeros(dims,2);
  ds(:) = randsample(D,dims*2);
  
  % Multivariate time series
  X = dat(ds(:,1),:,ss(:,1));
  Y = dat(ds(:,2),:,ss(:,2));
  
  W = []; % No conditional

  % Digitally filter the series?
  if to_filter > 0
    X = filter(a_coeff,b_coeff,X,[],2);
    Y = filter(a_coeff,b_coeff,Y,[],2);
  end
  % Chop off samples to avoid boundary effects from filtering, etc.
  X = X(:,seq)';
  Y = Y(:,seq)';
  
  % Compute our measures
  if is_granger
    [measure(r),dist] = granger_causality(X,Y,W,embedding,'none',false);
  else
    [measure(r),dist] = conditional_mutual_information(X,Y,W,'none',false);
  end
  
  pvals_LR(r) = compute_significance(measure(r), dist, 'lr');
  pvals_E(r) = compute_significance(measure(r), dist, 'exact',S);
  
  if verbose
    if mod(r,10) == 0
      fprintf('Completed run %d/%d.\n', r, R);
    end
  end
end

%% Plot results
col_LR = [1 0 0];
col_E = [0 0 0];
  
figure;
hold on;
plot([0 1], [0 1], 'k:');
ph1 = plot(sort(pvals_LR),linspace(0,1,R), '-', 'color', col_LR, 'linewidth', 1);
ph2 = plot(sort(pvals_E),linspace(0,1,R), '-', 'color', col_E, 'linewidth', 1);
legend([ph1 ph2],'LR test', 'Exact test','location', 'best');

fprintf('LR test FPR at %.2f: %.3g\n', alpha, mean(pvals_LR <= alpha) );
fprintf('Exact test FPR at %.2f: %.3g\n', alpha, mean(pvals_E <= alpha) );

if generate_tikz
  if ~exist('tikz','dir')
    mkdir('tikz');
  end
  tikz_fname = sprintf('tikz/FPR_hcp_test-%d.tikz', which_test);
  matlab2tikz( 'filename', tikz_fname, 'noSize', true );
end