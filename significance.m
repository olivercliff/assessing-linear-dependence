function [pval,stat] = significance(estimate,dist,test,S,correction_method_str)
%SIGNIFICANCE Significance testing for linear-dependence measures.
%   PVAL = SIGNIFICANCE(ESTIMATE,STATS) returns PVAL, the p-value for the
%   scalar ESTIMATE and statistics STATS using our exact hypothesis tests.
%   Both ESTIMATE and STATS should be computed from either the function <a href="matlab:help mvgc">mvgc</a>
%   <a href="matlab:help mvmi">mvmi</a> in this package.
%
%   PVAL = SIGNIFICANCE(...,TEST) uses either the 'exact' test or the 'lr' (likelihood
%   ratio) test.
%
%   PVAL = SIGNIFICANCE(...,S) if using the 'exact' test, specify the number of
%   surrogates used in the Monte Carlo sampling.
%
%   PVAL = SIGNIFICANCE(...,CORRECTION) optional string input specifying the
%   method used for bartlett correction:
%     - 'none' do not Bartlett correct.
%     - 'bartlett' Bartlett correction as per paper.
%     - 'roy' use full Multivariate Bartlett formula as per Roy (1989). Not currently
%             implemented.
%
%   [PVAL,STAT] = SIGNIFICANCE(...) also returns STAT, the statistic
%   compared against the null distribution.
%
%   See also <a href="matlab:help mvmi">mvmi</a>, <a href="matlab:help mvgc">mvgc</a>, <a href="matlab:help pcd">pcd</a>

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Oliver M. Cliff <oliver.m.cliff@gmail.com>,
%
% If you use this code for your research, please cite the following paper:
%
% Oliver M. Cliff, Leonardo Novelli, Ben D Fulcher, James M. Shine,
% Joseph T. Lizier, "Exact Inference of Linear Dependence for Multiple
% Autocorrelated Time Series," arXiv preprint arXiv:2003.03887 (2020).
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 3
  test = 'exact';
end

if strcmp(test,'lr')
  
  % Statistic is 2 * nested log ratio * number of samples (removed the
  % order of autoregression)
  stat = 2 * dist.N_o * dist.to_cmi(estimate);
  
  % Get p-value from survival function of chi-squared dist
  pval = 1 - chi2cdf( stat, dist.dof );
  
% (Our) exact test
elseif strcmp(test,'exact')
  
  if nargin < 4
    S = 5000;
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
    dist.N_eff = dist.N_eff./diag(dist.eta);
  end
  dist.N_eff = dist.N_eff - dist.cs;

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
  pval = mean( stat <= dist_exact );
end