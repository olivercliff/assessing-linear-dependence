function [X_tilde,Y_tilde,ar_order] = prewhitenAR(X,Y,ar_order,use_aic)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
if nargin < 3 || isempty(ar_order)
  
  if ~use_aic
    % Uses Burg's method/partial autocorrelation and chooses first value <
    % 1.96/sqrt(T-1)
    ar_order = order(X);
  else
    max_order = length(X)-1;

    Mdls = cell(max_order,1);
    aics = inf(max_order,1);
    for p = 1:max_order
      try
        Mdls{p} = ar(iddata(X),p,'approach','burg');
        aics(p) = aic(Mdls{p});
      catch
        break;
      end
    end
    [~,minaic] = min(aics);
    ar_order = minaic;
    Mdl = Mdls{ar_order};
  end
end

if ~exist('Mdl','var')
  Mdl = ar(iddata(X),ar_order,'approach','burg');
end

% Take residuals
X_tilde = resid(iddata(X),Mdl).OutputData;
Y_tilde = resid(iddata(Y),Mdl).OutputData;