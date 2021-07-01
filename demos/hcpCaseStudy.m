% HCP_CASE_STUDY   Runs the Human Connectome Project case study from the paper.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Oliver M. Cliff <oliver.m.cliff@gmail.com>,
%
% If you use this code for your research, please cite the following paper:
%
% Oliver M. Cliff, Leonardo Novelli, Ben D Fulcher, James M. Shine,
% Joseph T. Lizier, "Assessing the significance of directed and multivariate
% measures of linear dependence between time series," Phys. Rev. Research 3,
% 013145 (2021).
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

clear
close all

addpath(genpath('..'));

prewhiten = false;

% Download the HCP data (if it's not already there)
datafile = './data/hcp_rsfMRI.mat';

if ~exist(datafile,'file')
  article_url = '<a href = "https://figshare.com/articles/hcp_rsfMRI_mat/12084636">https://figshare.com/articles/hcp_rsfMRI_mat/12084636</a>';
  download_url = 'https://ndownloader.figshare.com/files/22217589';
  
  fprintf('HCP data not found in the expected location: %s.\n', datafile);
  
  fprintf('You can download it here: %s\n', article_url);
  prompt = 'Do you want me to download it now? Y/N [Y]: ';
  reply = input(prompt,'s');
  
  if isempty(reply)
    reply = 'Y';
  end
  
  if upper(reply) ~= 'Y'
    fprintf('OK, exiting program now.\n');
    return
  end
  
  if ~exist('data','dir')
    mkdir('data');
  end
  
  fprintf('OK, getting it now. This may take a while...\n');
  
  % Download the hcp_rsfMRI.mat file from figshare
  out_fname = websave(datafile,download_url);
  fprintf(['Done. Human Connectome Project rsfMRI data ',...
            'downloaded from figshare to:\n%s\n'],out_fname);

end

%% Load HCP data

fprintf('Loading data...\n');
load(datafile);
fprintf('Done.\n');

% Normalise (otherwise we have global effects in trends, etc.)
fprintf('Normalising and detrending...\n');
for i = 1:size(dat,3)
  s_bold = dat(:,:,i)';
  s_bold = detrend(s_bold);
  s_bold = zscore(s_bold);
  dat(:,:,i) = s_bold';
end

% Reproduce which HCP experiment from the paper?
% - 1: MI
% - 2: MI with MV dim 2
% - 3: GC, optimal embedding
% - 4: GC, high embedding
% - 5: GC with MV dim 2 (not in paper)
which_test = 4;

% Do you want to use the F-test or the asymptotic LR (chi^2) test 
f_test = false;

switch which_test
  case 1
    is_granger = false; % MI
    dims = 1; % 1D
  case 2
    is_granger = false; % MI
    dims = 2; % 2D
  case 3
    is_granger = true; % GC
    config.p = 'auto';
    config.q = 'auto';
    dims = 1; % 1D
  case 4
    is_granger = true; % GC
    config.p = '100';
    config.q = '100'; % Overly optimistic (pessimistic?) embedding
    dims = 1; % 1D
  case 5
    is_granger = true; % GC
    config.p = 'auto';
    config.q = 'auto';
    dims = 2; % 1D
  case 6
    is_granger = true; % GC
    config.p = 'auto';
    config.q = '1'; % auto-embedding
    dims = 1; % 1D
end

if f_test && dims > 1
  warning(['No known finite-sample distribution for ',...
            'F-tests with multivariate models. Using LR test instead.']);
  f_test = true;
end

%% Default params

seed = 1; % Kind of irrelevant since parfor is used
generate_tikz = false;

R = 1000; % Number of trials
surrogates = 5000; % MC sample size
alpha = 0.05; % Significance level
seq = 201:1000; % Sequence of the data to take (cut off first+last 200)

to_filter = 2; % Filter the BOLD data? (1 = FIR, 2 = IIR)

verbose = true;

%% Prelim stuff

fs = 1/0.72; % Sample time (s)

if to_filter == 1
%   fc = 0.2; % Low-pass (Hz)
  fc = [0.02 0.08]; % Passband (Hz)
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
measure_pw = zeros(R,1);

% ..and p-values
pvals_LR = zeros(R,1); % log-likelihood ratio test p-values
pvals_LR_pw = zeros(R,1); % log-likelihood ratio test p-values
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
  
  if prewhiten && dims == 1
    try
      [X_pw,Y_pw,W_pw] = prewhitenAR(X,Y,W);
    catch
      X_pw = X; Y_pw = Y; W_pw = W;
    end
  end
  
  % Compute our measures
  if is_granger
    [measure(r),pvals_E(r)] = mvgc(X,Y,W,'p',config.p,'q',config.q,'test','modified','surrogates',surrogates);
    [~,pvals_LR(r)] = mvgc(X,Y,W,'p',config.p,'q',config.q,'test','asymptotic');
    
    if prewhiten && dims == 1
      [measure_pw(r),pvals_LR_pw(r)] = mvgc(X_pw,Y_pw,W_pw,'p',config.p,'q',config.q,'test','asymptotic');
    end
  else
    [measure(r),pvals_E(r)] = mvmi(X,Y,W,'test','modified','surrogates',surrogates);
    [measure(r),pvals_LR(r)] = mvmi(X,Y,W,'test','asymptotic');
    
    if prewhiten && dims == 1
      [measure_pw(r),pvals_LR_pw(r)] = mvgc(X_pw,Y_pw,W_pw,'test','asymptotic');
    end
  end
  
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
if prewhiten && dims == 1
  ph3 = plot(sort(pvals_LR_pw),linspace(0,1,R), '--', 'color', col_LR, 'linewidth', 1);
end

if f_test
  legend([ph1 ph2],'Standard F-test', 'Modified $F$-test','location', 'best','interpreter','latex');
else
  legend([ph1 ph2],'$\chi^2$-test', 'Modified $F$-test','location', 'best','interpreter','latex');
end

fprintf('Exact test FPR at %.2f: %.4g\n', alpha, mean(pvals_E <= alpha) );
fprintf('LR test FPR at %.2f: %.4g\n', alpha, mean(pvals_LR <= alpha) );
if prewhiten && dims == 1
  fprintf('LR test (after prewhitening) FPR at %.2f: %.4g\n', alpha, mean(pvals_LR_pw <= alpha) );
end

if generate_tikz
  if ~exist('tikz','dir')
    mkdir('tikz');
  end
  tikz_fname = sprintf('tikz/FPR_hcp_test-%d.tikz', which_test);
  matlab2tikz( 'filename', tikz_fname, 'noSize', true );
end