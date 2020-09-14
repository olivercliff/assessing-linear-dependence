function [X_tilde,Y_tilde,W_tilde,orders] = prewhitenARMA(X,Y,W,p_max,q_max,fixed_order)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the ARMA(p,q) model of X)

if nargin < 6
  fixed_order = false;
end

% Get the model order
if nargin < 4 || isempty(p_max)
  p_max = 1;
end
if nargin < 5 || isempty(q_max)
  q_max = 1;
end

if ~fixed_order
  % Cache the log-likelihood of the estimated models
  logL = zeros(1+p_max,1+q_max);
  mdl = cell(1+p_max,1+q_max);
  for p = 0:p_max
    for q = 0:q_max
      % Try two ARMA fitting procedures (interior-point is more stable than sqp).
      % if both fail then the try-catch outside should set this run to NaN)
      mdl{p+1,q+1} = arima(p,0,q);
      mdl{p+1,q+1}.Constant = 0;
      try
        [mdl{p+1,q+1},~,logL(p+1,q+1)] = estimate(mdl{p+1,q+1},X,'Display','off');
      catch
        opts = optimset('fmincon');
        opts.Algorithm = 'interior-point'; % Try the interior-point rather than sqp algo
        [mdl{p+1,q+1},~,logL(p+1,q+1)] = estimate(mdl{p+1,q+1},X,'options',opts,'Display','off');
      end
    end
  end

  nparams = (1:p_max+1)'*(1:q_max+1);
  [~,bic] = aicbic(logL(:), nparams(:), length(X)*ones((p_max+1)*(q_max+1),1));

  % Use AIC
  [~,opt_id] = min(bic);
  [p_opt,q_opt] = ind2sub([1+p_max, 1+q_max],opt_id);

  est_mdl = mdl{p_opt,q_opt};

  orders = [p_opt-1,q_opt-1];
else
  mdl = arima(p_max,0,q_max);
  mdl.Constant = 0;
  
  try
    est_mdl = estimate(mdl,X,'Display','off');
  catch
    opts = optimset('fmincon');
    opts.Algorithm = 'interior-point';
    est_mdl = estimate(mdl,X,'options',opts,'Display','off');
  end
end

X_tilde = infer(est_mdl,X);
Y_tilde = infer(est_mdl,Y);

if nargin > 2 || isempty(W)
  W_tilde = zeros(size(X_tilde,1),size(W,2));
  for i = 1:size(W,2)
    W_tilde(:,i) = infer(est_mdl,W(:,i));
  end
else
  W_tilde = [];
end