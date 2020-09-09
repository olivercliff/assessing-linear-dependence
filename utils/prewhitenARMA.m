function [X_tilde,Y_tilde,W_tilde,orders] = prewhitenARMA(X,Y,W,p,q)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the ARMA(p,q) model of X)

% Get the model order
if p == 0 && q == 0
  X_tilde = X;
  Y_tilde = Y;
  W_tilde = W;
  orders = [0, 0];
  return;
end

% Try two ARMA fitting procedures (if both fail then the try-catch outside should set this run to NaN)
try
  mdl = arima(p,0,q);
  mdl = estimate(mdl,X,'Display','off');
catch
  opts = optimset('fmincon');
  opts.Algorithm = 'interior-point';
  mdl = estimate(mdl,X,'options',opts,'Display','off');
end

orders = [p,q];

X_tilde = infer(mdl,X);
Y_tilde = infer(mdl,Y);

if nargin > 2 || isempty(W)
  W_tilde = zeros(size(X_tilde,1),size(W,2));
  for i = 1:size(W,2)
    W_tilde(:,i) = infer(mdl,W(:,i));
  end
else
  W_tilde = [];
end