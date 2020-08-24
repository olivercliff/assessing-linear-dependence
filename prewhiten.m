function [X_tilde,Y_tilde,W_tilde,p] = prewhiten(X,Y,W,M)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
if nargin > 3
  p = order(X,M);
else
  p = order(X);
end

% Embed time series and infer AR parameters
[Xf,Yp,Xp] = embed(X,Y,p,p);
Yf = Y(p+1:end);

pi_B = Xf \ [ones(size(Xp,1),1), Xp];

X_tilde = Xf - [ones(size(Xp,1),1),Xp]*pi_B';
Y_tilde = Yf - [ones(size(Xp,1),1),Yp]*pi_B';

if nargin > 2 || isempty(W)
  W_tilde = W(p+1:end,:);
else
  W_tilde = [];
end