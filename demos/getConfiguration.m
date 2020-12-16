function config = getConfiguration(fig,subfig,exp)
  
  % start on fig 2
  fig = fig + 1;

  if subfig == 1
    config.to_filter = 1; % FIR
  else
    config.to_filter = 2; % IIR
  end

  config.T = 2^9; % Dataset length
  config.R = 1000; % Number of runs/trials
  config.S = 1000; % Number of samples for MC distribution
  config.alpha = 0.05; % significance level
  
  config.lm = false; % Include the AR parameters (or use IID, false)
  config.ar = true; % Include the AR parameters (or use IID, false)
  config.causal = false; % Non-causal (null model, should we be testing true negatives or false positives?)
  
  config.seed = 'shuffle'; % RNG seed (we didn't use one for the numerical sims)
  
  config.use_arma = false;
  
  filter_orders = 0:4:32;
  dims = 1:5;
  conditionals = 0:20:200;
  qs = [1,20,40:20:100];
  
  switch fig
    case 2
      config.is_granger = false; % GC or MI?
      
      config.dim_X = 1; % set X dimension
      config.dim_Y = 1; % set Y dimension
      config.dim_W = 0; % set W dimension
      
      config.whiten = false;
      
      if isempty(exp)
        config.filter_order = 8;
      else
        config.filter_order = filter_orders(exp);
      end
    case 3
      config.is_granger = false;
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 1;
      
      config.whiten = false;
      
      if isempty(exp)
        config.filter_order = 8;
      else
        config.filter_order = filter_orders(exp);
      end
    case 4
      config.is_granger = false;
      
      config.dim_X = 1;
      config.dim_Y = 1;
      
      config.filter_order = 8;
      
      config.whiten = false;
      
      if isempty(exp)
        config.dim_W = 100;
      else
        config.dim_W = conditionals(exp);
      end
    case 5
      config.is_granger = false;
      
      config.dim_W = 0;
      
      config.filter_order = 8;
      
      config.whiten = false;
      
      if isempty(exp)
        config.dim_X = 3;
        config.dim_Y = 3;
      else
        config.dim_X = dims(exp);
        config.dim_Y = dims(exp);
      end
    case 6
      config.is_granger = true;
      
      config.p = 'auto';
      config.q = 'auto';
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 0;
      
      config.whiten = false;
      
      if isempty(exp)
        config.filter_order = 8;
      else
        config.filter_order = filter_orders(exp);
      end
    case 7
      config.is_granger = true;
      
      config.p = 'auto';
      
      config.filter_order = 8;
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 0;
      
      config.whiten = false;
      
      if isempty(exp)
        config.q = '20';
      else
        config.q = num2str(qs(exp));
      end
    case 8
      config.is_granger = true;
      config.p = 'auto';
      config.q = '1';
      config.filter_order = 8;
      config.whiten = false;
      
      config.dim_X = 3;
      config.dim_Y = 3;
      config.dim_W = 0;
      
      if isempty(exp)
        config.dim_X = 3;
        config.dim_Y = 3;
      else
        config.dim_X = dims(exp);
        config.dim_Y = dims(exp);
      end
    case 9
      config.is_granger = false;
      config.filter_order = 8;
      config.whiten = false;
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 0;
      
      if isempty(exp)
        config.T = 2^9;
      else
        config.T = 2^(5+2*exp);
      end
    case 10
      config.is_granger = true;
      
      config.filter_order = 8;
      
      config.whiten = false;
      
      config.p = 'auto';
      config.q = 'auto';
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 0;
      
      if isempty(exp)
        config.T = 2^9;
      else
        config.T = 2^(5+2*exp);
      end
    case 11
      config.is_granger = false; % GC or MI?
      
      config.dim_X = 1; % set X dimension
      config.dim_Y = 1; % set Y dimension
      config.dim_W = 0; % set W dimension
      
      config.whiten = true;
      config.whiten_both = false;
      
      if isempty(exp)
        config.filter_order = 8;
      else
        config.filter_order = filter_orders(exp);
      end
    case 12
      config.is_granger = true;
      
      config.p = 'auto';
      config.q = 'auto';
      
      config.dim_X = 1;
      config.dim_Y = 1;
      config.dim_W = 0;
      
      config.whiten = true;
      config.whiten_both = false;
      
      if isempty(exp)
        config.filter_order = 8;
      else
        config.filter_order = filter_orders(exp);
      end
  end
  
  config.is_pc = false;
end