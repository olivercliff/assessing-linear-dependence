function numericalEvaluation(figure_or_config,which_exp,output_file,plot_results,parallel)

% NUMERICAL_EVALUATION   Runs simulations of all numerical evalutions from our paper.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Oliver M. Cliff <oliver.m.cliff@gmail.com>
%
% If you use this code for your research, please cite the following paper:
%
% Oliver M. Cliff, Leonardo Novelli, Ben D Fulcher, James M. Shine,
% Joseph T. Lizier, "Exact Inference of Linear Dependence for Multiple
% Autocorrelated Time Series," arXiv preprint arXiv:2003.03887 (2020).
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if ~exist('mvgc.m','file')
  addpath('..');
  addpath('../utils/');
end

if nargin < 5
  parallel = false;
  if nargin < 4
    plot_results = true;
    if nargin < 3
      output_file = [];
      if nargin < 2
        which_exp = [];
      end
    end
  end
end

%% Reproduce a figure from the paper? Else, choose the params below

verbose = true;

if ischar(figure_or_config)
  figure_opts = {'1a','1b';
                 '2a','2b';
                 '3a','3b';
                 '4a','4b';
                 '5a','5b';
                 '6a','6b';
                 '7a','7b';
                 '8a','8b';
                 '9a','9b';
                 '10a','10b';
                 '11a','11b'};

  fig_id = strcmp(figure_opts,figure_or_config);

  if any(fig_id(:))
    % Is the input a valid figure (1a-7b)...?
    [fig,subfig] = find(fig_id);
    config = getConfiguration(fig,subfig,which_exp);
  end
  
  fprintf('Performing numerical simulations from Fig. %s.\n', figure_or_config);
elseif isstruct(figure_or_config)
  config = figure_or_config;
  fprintf('User-supplied configuration.\n');
else
  fprintf('Unknown input.\n');
end

fprintf('Using the following params:\n');
disp(config);

%% Set up filters and simulator

univariate = config.dim_X == 1 && config.dim_Y == 1;

if config.whiten && ~univariate
  fprintf('The pre-whitening option is only available for univariate time series. Exiting.\n');
  return;
end

if config.is_pc
  computeMeasure = @(X,Y,W,varargin) pcorr(X,Y,W,varargin{:});
else
  if config.is_granger
    computeMeasure = @(X,Y,W,varargin) mvgc(X,Y,W,'p',config.p,'q',config.q,varargin{:});
  else
    computeMeasure = @(X,Y,W,varargin) mvmi(X,Y,W,varargin{:});
  end
end

if config.to_filter && config.filter_order > 0
  if config.to_filter == 1
    % FIR filter
    b_coeff = fir1(config.filter_order, 0.5);
    a_coeff = 1;
  elseif config.to_filter == 2
    % Butterworth (IIR) filter
    [b_coeff, a_coeff] = butter(config.filter_order, 0.5);
  end
else
  b_coeff = nan;
  a_coeff = nan;
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
  phi_XY = 0.2104 .* eye(config.dim_X,config.dim_Y);
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

fprintf('Using the following tests:\n');

fprintf('\t1. The modified F-test (ours)\n');
if config.whiten
  pvals_E = nan(config.R,1); % exact test (no prewhitening)
  pvals_F = nan(config.R,4); % 4 different approaches to prewhitening
  pw_ar_orders = nan(config.R,1);
  pw_arma_orders = nan(config.R,2);
  pw_timer = nan(config.R,4);
  fprintf('\t2. The F-test after prewhitening\n');
else
  fprintf('\t2. The asymptotic LR (chi-square) test\n');
  pvals_E = nan(config.R,1); % exact test
  pvals_chi2 = nan(config.R,1); % LR test
  eta_mean = nan(config.R,1);
  
  % Save the F-test results if univariate
  if univariate
    pvals_F = nan(config.R,1); % F-test
    fprintf('\t3. The F-test\n');
  end
end

timer = nan(config.R,1);

rng(config.seed);

fprintf('Running simulations...');
if ~verbose
  fprintf(' (turn on verbose to see progress)');
else
  fprintf(' (turn off verbose to suppress output)');
end
fprintf('\n');

if ~parallel
  parforargs = 0;
else
  parforargs = inf;
end

% Run sims
% parfor (r = 1:config.R, parforargs)
for r = 1:config.R
  
  if config.to_filter == 3 && univariate && config.dim_W == 0
    phi = generateARMAParams(config.filter_order,0);
    mdl = arima('AR',phi,'constant',0,'variance',1);
    X = simulate(mdl,config.T);
    Y = simulate(mdl,config.T);
    W = [];
  else
    Z = Sigma*randn(M,config.T);
    if config.ar
      for t = 2:config.T
        Z(:,t) = Phi*Z(:,t-1) + Z(:,t);
      end
    elseif config.lm
      Z(:,1) = rand(M,1);
      for t = 2:config.T
        Z(:,t) = config.r.*Z(:,t-1).*(1-Z(:,t-1));
      end
    end

    % Partition the dataset (Z) into multiple time series (X,Y, and conditional W)
    X = Z(p_X,:)';
    Y = Z(p_Y,:)';
    W = Z(p_W,:)';

    % Filter the data to induce higher autocorrelation (if opted)
    if config.to_filter > 0 && config.filter_order > 0
      X = filter(b_coeff,a_coeff,X,[],1);
      Y = filter(b_coeff,a_coeff,Y,[],1);
      W = filter(b_coeff,a_coeff,W,[],1);
    end
  end
  
  % Normalise the data (N.B. a column is a single time series)
  X = zscore(detrend(X));
  Y = zscore(detrend(Y));
  if ~isempty(W)
    W = zscore(detrend(W));
  end
  
  % Our modified F-test
  tic;
  [measure,pvals_E(r),~,stats] = computeMeasure(X,Y,W,'test','modified','taperMethod','tukey');
  timer(r) = toc;
  eta_mean(r) = mean(diag(stats.N_e));
  
  if config.whiten > 0
    
    try
      tic;
      [X_pw_ar1,Y_pw_ar1,W_pw_ar1] = prewhitenAR(X,Y,W,1,1,config.whiten_both);
      pw_timer(r,1) = toc;
      [~,pvals_F(r,1)] = computeMeasure(X_pw_ar1,Y_pw_ar1,W_pw_ar1,...
                                        'test','finite');
    catch
      warning('Run %i failed to learn AR(1) model\n',r);
    end
    
    try
      tic;
      [X_pw_arma11,Y_pw_arma11,W_pw_arma11] = prewhitenARMA(X,Y,W,1,1,true);
      pw_timer(r,2) = toc;
      [~,pvals_F(r,2)] = computeMeasure(X_pw_arma11,Y_pw_arma11,W_pw_arma11,...
                                        'test','finite');
    catch
      warning('Run %i failed to learn ARMA(1,1) model\n',r);
    end
    
    try
      tic;
      [X_pw_arp,Y_pw_arp,W_pw_arp,pw_ar_orders(r)] = prewhitenAR(X,Y,W,[],[],config.whiten_both);
      pw_timer(r,3) = toc;
      [~,pvals_F(r,3)] = computeMeasure(X_pw_arp,Y_pw_arp,W_pw_arp,...
                                        'test','finite');
    catch
      warning('Run %i failed to learn AR(p) model\n',r);
    end
    
    try
      tic;
      [X_pw_armapq,Y_pw_armapq,W_pw_armapq,pw_arma_orders(r,:)] = prewhitenARMA(X,Y,W,config.arma_p_max,config.arma_q_max,false);
      pw_timer(r,4) = toc;
      [~,pvals_F(r,4)] = computeMeasure(X_pw_armapq,Y_pw_armapq,W_pw_armapq,...
                                        'test','finite');
    catch
      warning('Run %i failed to learn ARMA(p,q) model\n',r);
    end
  else
    pvals_chi2(r) = significance(measure,stats,'test','asymptotic');
    
    if univariate
      % F-test
      pvals_F(r) = significance(measure,stats,'test','finite');
    end
  end
  
  if verbose
    if mod(r,10) == 0
      fprintf('Completed run %d/%d. Size (at 5%%):\n', r, config.R);
      fprintf('\tModified F-test: %.3g\n', mean(pvals_E(1:r) <= 0.05 ));
      if config.whiten
        fprintf('\tF-test (AR(1) whitened): %.3g\n', mean(pvals_F(1:r,1) <= 0.05 ));
        fprintf('\tF-test (ARMA(1,1) whitened): %.3g\n', mean(pvals_F(1:r,2) <= 0.05 ));
        fprintf('\tF-test (AR(p) whitened): %.3g\n', mean(pvals_F(1:r,3) <= 0.05 ));
        fprintf('\tF-test (ARMA(p,q) whitened): %.3g\n', mean(pvals_F(1:r,4) <= 0.05 ));
      else
        fprintf('\tChi^2-test: %.3g\n', mean(pvals_chi2(1:r) <= 0.05 ));
        if univariate
          fprintf('\tF-test: %.3g\n', mean(pvals_F(1:r) <= 0.05 ));
        end
      end
    end
  end
end

if ~isempty(output_file)
  if config.whiten
    save(output_file,'config','pvals_E','pvals_F','pw_ar_orders','pw_arma_orders','eta_mean');
  else
    if univariate
      save(output_file,'config','pvals_E','pvals_chi2','pvals_F','eta_mean');
    else
      save(output_file,'config','pvals_E','pvals_chi2');
    end
  end
end
  
if plot_results
  col_LR = [1 0 0];
  col_E = [0 0 0];

  figure;
  hold on;
  plot([0 1], [0 1], 'k--');
  ph(1) = plot(sort(pvals_E),linspace(0,1,config.R), '-', 'color', col_E, 'linewidth', 1);
  ph(2) = plot(sort(pvals_chi2),linspace(0,1,config.R), '--', 'color', col_LR, 'linewidth', 1);
  
  labels{1} = 'Modified $\Lambda$-test';
  labels{2} = '$\chi^2$-test';

  fprintf('Modified Lambda-test FPR at %d%% significance: %.3g\n', config.alpha*100, mean(pvals_E <= config.alpha) );
  fprintf('chi2-test FPR at %d%% significance: %.3g\n', config.alpha*100, mean(pvals_chi2 <= config.alpha) );
  if univariate
    ph(3) = plot(sort(pvals_F),linspace(0,1,config.R), '-', 'color', col_LR, 'linewidth', 1);
    
    labels{3} = '$F$-test';
    
    fprintf('F-test FPR at %d%% significance: %.3g\n', config.alpha*100, mean(pvals_F <= config.alpha) );
    if config.whiten
      ph(4) = plot(sort(pvals_chi2_pw),linspace(0,1,config.R), '-.', 'color', col_LR, 'linewidth', 1);
      ph(5) = plot(sort(pvals_F_pw),linspace(0,1,config.R), ':', 'color', col_LR, 'linewidth', 1);
      
      labels{4} = '$\chi^2$-test (whitened)';
      labels{5} = '$F$-test (whitened)';
      
      fprintf('Pre-whitened Chi-square test FPR at %d%% significance: %.3g\n', config.alpha*100, mean(pvals_chi2_pw <= config.alpha) );
      fprintf('Pre-whitened F-test FPR at %d%% significance: %.3g\n', config.alpha*100, mean(pvals_F_pw <= config.alpha) );
    end
  end
  legend(ph,labels,'interpreter','latex', 'location', 'best');
end