function config = get_configuration(fig,subfig)
  
  % start on fig 2
  fig = fig + 1;

  if subfig == 1
    config.to_filter = 1; % FIR
  else
    config.to_filter = 2; % IIR
  end
  config.filter_order = 8; % 8th order

  config.T = 2^9; % Dataset length
  config.R = 1000; % Number of runs/trials
  config.S = 1000; % Number of samples for MC distribution
  config.alpha = 0.05; % significance level
  
  config.ar = true; % Include the AR parameters (or use IID, false)
  config.causal = false; % Non-causal (null model, should we be testing true negatives or false positives?)
  
  config.seed = now; % RNG seed (we didn't use one for the numerical sims)
  
  switch fig
    case 2
      config.is_granger = false; % GC or MI?
      config.embedding = nan;
      
      config.dim_X = 1; % set X dimension
      config.dim_Y = 1; % set Y dimension
      config.dim_W = 0; % set W dimension
    case 3
      config.is_granger = false;
      config.embedding = nan;
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 1;
    case 4
      config.is_granger = false;
      config.embedding = nan;
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 100;
    case 5
      config.is_granger = false;
      config.embedding = nan;
      
      config.dim_X = 3;
      config.dim_Y = 3;
      config.dim_W = 0;
    case 6
      config.is_granger = true;
      config.embedding = [-1 -1];
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 0;
    case 7
      config.is_granger = true;
      config.embedding = [-1 20];
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 0;
    case 8
      config.is_granger = true;
      config.embedding = [-1 -1];
      
      config.dim_X = 3;
      config.dim_Y = 3;
      config.dim_W = 0;
  end
end