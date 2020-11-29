% Choose settings from the paper (or '' for your own custom configuration)
myfig = '10a';
exp = 1;

addpath ../utils/

if isempty(myfig)
  config.T = 2^11; % Dataset length
  config.R = 1000; % Number of runs/trials
  config.S = 1000; % Number of samples for MC distribution
  config.dim_X = 1; % set X dimension
  config.dim_Y = 1; % set Y dimension
  config.dim_W = 0; % set W dimension

  config.to_filter = 0; % None=0, MA/FIR=1, ARMA/IIR=2
  config.filter_order = 0; % FIR/IIR Filter order
  config.ar = false; % Include the AR parameters (or use IID, false)
  config.lm = true;
  config.r = 3.65;

  config.causal = false; % Non-causal (null model, should we be testing true negatives or false positives?)

  config.is_granger = false; % Granger causality (not CMI)
  config.p = 'auto'; % Optimal embedding, set to numeric for a specific embedding
  config.q = 'auto';
  
  config.whiten = false; % Granger causality (not CMI)
  config.arma_p_max = 5;
  config.arma_q_max = 5;

  config.is_pc = false;

  config.alpha = 0.05; % significance level

  config.seed = now; % RNG seed
  
  numericalEvaluation(config);
else
%   numericalEvaluation(myfig,exp,['./results/fig-' myfig '_exp-' num2str(exp) '.mat'],false,false);
  numericalEvaluation(myfig,exp,[],true,false);
end