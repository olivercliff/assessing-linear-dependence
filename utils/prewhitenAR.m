function [X_tilde,Y_tilde,ar_order] = prewhitenAR(X,Y,ar_order,method)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the AR(p) model of X)

% Get the model order
if nargin < 4
  method = 'AICc';
end

if nargin < 3 || isempty(ar_order)
  
  if strcmp(method,'Burgs')
    % Uses Burg's method/partial autocorrelation and chooses first value <
    % 1.96/sqrt(T-1)
    ar_order = order(X);
  elseif strcmp(method,'AICc') || strcmp(method,'BIC')
    
    max_order = round(length(X)/2);

    Mdls = cell(max_order,1);
    aics = inf(max_order,1);
    XI = iddata(X);
    for p = 1:max_order
      try
        Mdls{p} = ar(XI,p);
        aics(p) = aic(Mdls{p},method);
      catch
        break;
      end
    end
    [~,ar_order] = min(aics);
    Mdl = Mdls{ar_order};
  else
    error('Unknown method: %s',method);
  end
end

if ~exist('Mdl','var')
  Mdl = ar(iddata(X),ar_order);
end

% Take residuals
X_tilde = resid(iddata(X),Mdl).OutputData;
Y_tilde = resid(iddata(Y),Mdl).OutputData;