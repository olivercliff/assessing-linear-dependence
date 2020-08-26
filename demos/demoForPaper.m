
% Choose settings from the paper (or '' for your own custom configuration)
myfig = '8b';
exp = 8;

if isempty(myfig)
  config.T = 2^9; % Dataset length
  config.R = 1000; % Number of runs/trials
  config.S = 1000; % Number of samples for MC distribution
  config.dim_X = 1; % set X dimension
  config.dim_Y = 1; % set Y dimension
  config.dim_W = 0; % set W dimension

  config.to_filter = 2; % None=0, FIR=1, IIR=2
  config.filter_order = 8; % 8th order
  config.ar = true; % Include the AR parameters (or use IID, false)

  config.causal = false; % Non-causal (null model, should we be testing true negatives or false positives?)

  config.is_granger = true; % Granger causality (not CMI)
  config.p = 'auto'; % Optimal embedding, set to numeric for a specific embedding
  config.q = 'auto';

  config.is_pc = false;

  config.alpha = 0.05; % significance level

  config.seed = now; % RNG seed
  numericalEvaluation(config);
else
%   numericalEvaluation(myfig,exp,['./results/fig-' myfig '_exp-' num2str(exp) '.mat'],false,false);
  numericalEvaluation(myfig,exp,[],true,true);
end

