function [X_tilde,Y_tilde,arma_orders] = prewhitenARMA(X,Y,ar_order,ma_order,aic_or_bic)
% Takes in two vectors, X and Y, (optionally a third, W) and outputs the
% pre-whitened time series (based on the ARMA(p,q) model of X)

if nargin < 6
  aic_or_bic = 'aic';
end

if nargin < 4 || isempty(ar_order) || isempty(ma_order)
  ar_max = length(X) - 1;
  ma_max = ar_max;
  
  logL = zeros(1+ar_max,1+ma_max);
  Mdls = cell(1+ar_max,1+ma_max);
  for p = 0:ar_max
    for q = 0:ma_max
      
      Mdls{p+1,q+1} = arima(p,0,q);
      Mdls{p+1,q+1}.Constant = 0;
      try
        [Mdls{p+1,q+1},~,logL(p+1,q+1)] = estimate(Mdls{p+1,q+1},X,'Display','off');
      catch
        Mdls = Mdls(1:p-1,1:q-1);
        logL = logL(1:p-1,1:q-1);
        break;
      end
    end
  end

  nparams = (1:p)'*(1:q);
  [aic,bic] = aicbic( logL(:), nparams(:), repmat(length(X),p*q,1) );

  if strcmp(aic_or_bic,'bic')
    [~,opt_id] = min(bic);
  else
    [~,opt_id] = min(aic);
  end
  [ar_order,ma_order] = ind2sub(size(nparams),opt_id);

  Mdl = Mdls{ar_order,ma_order};

  arma_orders = [ar_order-1,ma_order-1];
else
  Mdls = arima(ar_order,0,ma_order);
  Mdls.Constant = 0;
  
  try
    Mdl = estimate(Mdls,X,'Display','off');
  catch
    opts = optimset('fmincon');
    opts.Algorithm = 'interior-point';
    Mdl = estimate(Mdls,X,'options',opts,'Display','off');
  end
end

X_tilde = infer(Mdl,X);
Y_tilde = infer(Mdl,Y);