function [X_tilde,Y_tilde,W_tilde,orders] = prewhitenARMA(X,Y,W,p,q)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
while q > 0 || p > 0
  try
    opts = optimset('fmincon');
    opts.Algorithm = 'interior-point';
    mdl = arima(p,0,q);
    mdl = estimate(mdl,X,'options',opts,'Display','off');
    break;
  catch
    warning('Reducing filter order.\n');
    q = q-1;
    p = p-1;
  end
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