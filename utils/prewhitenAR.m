function [X_tilde,Y_tilde,W_tilde,pX,pY] = prewhitenAR(X,Y,W,pX,pY,pw_both)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
if pw_both
  if nargin < 4 || isempty(pX)
    pX = order(X);
    pY = order(Y);
  end
  MdlX = ar(iddata(X),pX,'approach','burg');
  MdlY = ar(iddata(Y),pY,'approach','burg');
else
  if nargin < 4 || isempty(pX)
    pX = order(X);
  end
  MdlX = ar(iddata(X),pX,'approach','burg');
  MdlY = MdlX;
end

% Take residuals
X_tilde = X-predict(MdlX,iddata(X),1).OutputData;
Y_tilde = Y-predict(MdlY,iddata(Y),1).OutputData;

if nargin > 2 || isempty(W)
  W_tilde = zeros(size(X_tilde,1),size(W,2));
  for i = 1:size(W,2)
    W_tilde(:,i) = infer(MdlX,W(:,i));
  end
else
  W_tilde = [];
end