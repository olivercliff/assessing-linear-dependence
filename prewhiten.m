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

XpC = [ones(size(Xp,1),1), Xp];
YpC = [ones(size(Xp,1),1), Yp];

pi_B = XpC\Xf;

X_tilde = Xf - XpC*pi_B;
Y_tilde = Yf - YpC*pi_B;

% pi_B = parcorr(X,'numLags',p);
% 
% X_tilde = Xf;
% Y_tilde = Yf;
% for i = 1:p
%   X_tilde = X_tilde - Xp(:,i) * pi_B(i+1);
%   Y_tilde = Y_tilde - Yp(:,i) * pi_B(i+1);
% end

if nargin > 2 || isempty(W)
  W_tilde = W(p+1:end,:);
else
  W_tilde = [];
end