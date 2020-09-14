function [X_tilde,Y_tilde,W_tilde,p] = prewhitenAR(X,Y,W,p)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
if nargin < 4
  [p,pacf] = order(X);
  Mdl = arima('AR',pacf);
else
  Mdl = arima(p,0,0);
end

% infer AR parameters
Mdl = estimate(Mdl,X,'Display','off');

% Take residuals
X_tilde = infer(Mdl,X);
Y_tilde = infer(Mdl,Y);

if nargin > 2 || isempty(W)
  W_tilde = zeros(size(X_tilde,1),size(W,2));
  for i = 1:size(W,2)
    W_tilde(:,i) = infer(Mdl,W(:,i));
  end
else
  W_tilde = [];
end