function [pvalue,dist] = compute_significance(estimate,dist,test,S,correction_method_str)

% Log-likelihood ratio test
if strcmp(test,'lr')
  
  % Statistic is 2 * nested log ratio * number of samples (removed the
  % order of autoregression)
  stat = 2 * dist.N_o * dist.to_cmi(estimate);
  
  % Get p-value from survival function of chi-squared dist
  pvalue = 1 - chi2cdf( stat, dist.dof );
  
% (Our) exact test
elseif strcmp(test,'exact')
  
  if nargin < 4
    S = 1000;
  end
  if nargin < 5
    correction = 1;
  else
    switch correction_method_str
      case 'none'
        correction = 0;
      case 'bartlett'
        correction = 1;
      case 'roy'
        correction = 2;
    end
  end

  if correction > 1 && ~dist.mv
    warning('Roy''s correction only available if CMI is computed with full multivariate setting. Use, e.g., conditional_mutual_information(X,Y,W,''none'',true)\n');
  end

  % Initial effective sample size (remove the order of autoregression)
  dist.N_eff = dist.N_o;
  
  % Bartlett-corrected effective sample size
  if correction > 0
    dist.N_eff = dist.N_eff./diag(dist.nu);
  end
  dist.N_eff = dist.N_eff - dist.c;

  if any(dist.N_eff < 1)
    sum_lt1 = sum(dist.N_eff < 1);
    warning('F-statistics with effective sample size less than one: %f.\n', sum_lt1);
  end
  
  % 2nd input parameter to F-distribution (or only param to Student's t)
  dist.d_2 = dist.N_eff - 2;
  
  % Monte carlo sample the t-distributed random variables
  t_rvs = zeros(S,dist.dof);
  for i = 1:dist.dof
    t_rvs(:,i) = trnd(dist.d_2(i),[S,1]);
  end
  
  % Compute p-value
  if dist.cmi
    % Conditional mutual information (sums of log-F dist. RVs)
    f_rvs = t_rvs.^2;
    logf_rvs = log(f_rvs./dist.d_2'+1);
    dist_exact = sum(logf_rvs,2);
    stat = 2*dist.to_cmi(estimate);
  else
    % Partial correlation (sums of Student's t dist. RVs)
    dist_exact = sum(t_rvs,2);
    stat = estimate;
  end
  
  % Proportion of surrogates less than statistic
  pvalue = mean( stat <= dist_exact );
end