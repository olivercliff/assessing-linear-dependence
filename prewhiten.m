function [X_tilde,Y_tilde,W_tilde,p] = prewhiten(X,Y,W)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
p = order(X);

% Embed time series and infer AR parameters
[Xf,Yp,Xp] = embed(X,Y,p,p);
Yf = Y(p+1:end);

pi_B = Xf \ Xp;

X_tilde = Xf - Xp*pi_B';
Y_tilde = Yf - Yp*pi_B';
W_tilde = W(p+1:end,:);