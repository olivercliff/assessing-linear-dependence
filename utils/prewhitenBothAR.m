function [X_tilde,Y_tilde,W_tilde,pX,pY] = prewhitenBothAR(X,Y,W,pX,pY)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
if nargin < 4
  [pX,pacfX] = order(X);
  [pY,pacfY] = order(Y);
  MdlX = arima('AR',pacfX);
  MdlY = arima('AR',pacfY);
else
  MdlX = arima(pX,0,0);
  MdlY = arima(pY,0,0);
end

% infer AR parameters
MdlX = estimate(MdlX,X,'Display','off');
MdlY = estimate(MdlY,Y,'Display','off');

% Take residuals
X_tilde = infer(MdlX,X);
Y_tilde = infer(MdlY,Y);

if nargin > 2 || isempty(W)
  W_tilde = zeros(size(X_tilde,1),size(W,2));
  for i = 1:size(W,2)
    W_tilde(:,i) = infer(MdlX,W(:,i));
  end
else
  W_tilde = [];
end