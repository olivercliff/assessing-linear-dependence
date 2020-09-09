
% Choose settings from the paper (or '' for your own custom configuration)
% myfig = '8b';
myfig = '10a';
exp = 2;

if isempty(myfig)
  config.T = 2^5; % Dataset length
  config.R = 20; % Number of runs/trials
  config.S = 1000; % Number of samples for MC distribution
  config.dim_X = 1; % set X dimension
  config.dim_Y = 1; % set Y dimension
  config.dim_W = 0; % set W dimension

  config.to_filter = 3; % None=0, MA/FIR=1, ARMA/IIR=2, AR=3
  config.filter_order = 10; % 8th order
  config.ar = false; % Include the AR parameters (or use IID, false)

  config.causal = false; % Non-causal (null model, should we be testing true negatives or false positives?)

  config.is_granger = false; % Granger causality (not CMI)
  config.p = 'auto'; % Optimal embedding, set to numeric for a specific embedding
  config.q = 'auto';
  
  config.whiten = true; % Granger causality (not CMI)

  config.is_pc = false;

  config.alpha = 0.05; % significance level

  config.seed = now; % RNG seed
  
  numericalEvaluation(config);
else
%   numericalEvaluation(myfig,exp,['./results/fig-' myfig '_exp-' num2str(exp) '.mat'],false,false);
  numericalEvaluation(myfig,exp,[],true,false);
end