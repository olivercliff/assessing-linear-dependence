function [X_tilde,Y_tilde,W_tilde,p] = prewhiten_arma(X,Y,W,p,q)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
while q > 1
  try
    mdl = arima(p,0,q);
    mdl = estimate(mdl,X);
    break;
  catch
    q = q-1;
    p = p-1;
  end
end

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