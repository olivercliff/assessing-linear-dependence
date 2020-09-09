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
[Xf,Xp,Yp] = embed(X,p,Y,p);
Yf = Y(p+1:end);

XpC = [ones(size(Xp,1),1), Xp];
YpC = [ones(size(Xp,1),1), Yp];

pi_B = XpC\Xf;

X_tilde = Xf - XpC*pi_B;
Y_tilde = Yf - YpC*pi_B;

if nargin > 2 || isempty(W)
  W_tilde = zeros(size(X_tilde,1),size(W,2));
  for i = 1:size(W,2)
    [Wf,Wp] = embed(W(:,i),p);
    WpC = [ones(size(Xp,1),1), Wp];
    W_tilde(:,i) = Wf - WpC*pi_B;
  end
else
  W_tilde = [];
end