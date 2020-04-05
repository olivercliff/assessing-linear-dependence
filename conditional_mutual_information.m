function [cmi,stats] = conditional_mutual_information(X,Y,W,taper_method_str,mv_bartlett_correction)

if nargin < 5
  mv_bartlett_correction = false;
end

if nargin < 4
  taper = 0;
else
  switch taper_method_str
    case 'none'
      taper = 0;
    case 'tukey'
      taper = 1;
    case 'parzen'
      taper = 2;
    case 'bartlett'
      taper = 3;
  end
end

[pr,nu,c] = compute_partial_correlations(X,Y,W,taper,mv_bartlett_correction);

% Outputs for computing the significance (variance estimation and number of
% condtiionals)
stats.nu = nu;
stats.c = c;
stats.mv = mv_bartlett_correction;
stats.dof = length(c);
stats.N_o = length(X);
stats.to_cmi = @(x) x;

stats.cmi = true;

cmi = -0.5*sum(log(1-pr.^2));

end