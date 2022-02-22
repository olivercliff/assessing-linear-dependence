function [pval,dist] = significance(estimate,stats,varargin)
%SIGNIFICANCE Significance testing for linear-dependence measures.
%   PVAL = SIGNIFICANCE(ESTIMATE,STATS) returns PVAL, the p-value for the
%   scalar ESTIMATE and statistics STATS using our exact hypothesis tests.
%   Both ESTIMATE and STATS should be computed from either the function <a href="matlab:help mvgc">mvgc</a>
%   <a href="matlab:help mvmi">mvmi</a> in this package.
%
%   [...] = SIGNIFICANCE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values.  Valid parameters are the
%   following:
%
%         Parameter                   Value
%          'test'                     'modified' (the default) uses our
%                                     modified lambda-test,
%                                     'finite' uses the typical two-tail
%                                     t-test and 'asymptotic uses the
%                                     chi-squared test.
%          'surrogates'               Integer denoting the number of
%                                     surrogates used in generating the
%                                     numerical null distributions.
%          'varianceEstimator'        'bartlett' (default) uses Bartlett's
%                                     formula assuming no
%                                     cross-correlations, 'roy' makes no
%                                     assumptions about cross-correlations.
%                                     (N.B. this needs a measure to be run
%                                     with param 'multivariateBartlett' set
%                                     to true.)
%          'RV'                       Which random variables to use
%                                     (F-distributed or beta-distributed)
%
%   See also <a href="matlab:help mvmi">mvmi</a>, <a href="matlab:help mvgc">mvgc</a>, <a href="matlab:help pcd">pcd</a>

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Oliver M. Cliff <oliver.m.cliff@gmail.com>,
%
% If you use this code for your research, please cite the following paper:
%
% Oliver M. Cliff, Leonardo Novelli, Ben D Fulcher, James M. Shine,
% Joseph T. Lizier, "Assessing the significance of directed and multivariate
% measures of linear dependence between time series," Phys. Rev. Research 3,
% 013145 (2021).
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

parser = inputParser;

addRequired(parser,'estimate',@isnumeric);
addRequired(parser,'stat',@isstruct);
addOptional(parser,'RV','B',@(x) any(validatestring(x,{'F','B'})));

parser = parseParameters(parser,estimate,stats,varargin{:});

S = parser.Results.surrogates;


if (~stats.cmi && ~strcmp(parser.Results.test,'modified')) || strcmp(parser.Results.test,'asymptotic')
  
  if stats.cmi
    % LR statistic is 2 * nested log ratio * numer of samples (removed the
    % order of autoregression)
    estimate = 2*stats.to_cmi(estimate);

    % Get p-value from quantile function of chi-squared dist
    pval = 1-chi2cdf(estimate*stats.N_o,stats.dof);

    if nargout > 1
      dist = chi2inv(linspace(0,1,S),stats.dof)./stats.N_o;
    end
  else
    t = estimate.*sqrt((stats.N_o-2)./(1-estimate.^2));

    % Get p-value from quantile function of F-dist
    pval = 1-fcdf(t.^2,1,stats.N_o-2);
    
    if nargout > 1
      dist = tinv(linspace(0,1,S),stats.N_o-2)./sqrt(stats.N_o-2);
    end
  end
  
% (Our) exact test
elseif strcmp(parser.Results.test,'finite')
  
  % Monte carlo sample the F/t-distributed random variables
  stats.d_2 = stats.N_o - stats.cs - 2;
  
  l_rvs = zeros(S,stats.dof);
  for i = 1:stats.dof
    l_rvs(:,i) = trnd(stats.d_2(i),[S,1]).^2 ./ stats.d_2(i);
  end
  
  dist = log(l_rvs + 1 );
  dist = sum(dist,2);
    
  estimate = 2*stats.to_cmi(estimate);
    
    % Proportion of surrogates less than statistic
  pval = mean(estimate<=dist);
    
  if nargout > 1
    dist = sort(dist);
  end
      
elseif strcmp(parser.Results.test,'modified')
  
  switch parser.Results.varianceEstimator
    case 'none'
      correction = 0;
    case 'bartlett'
      correction = 1;
    case 'roy'
      correction = 2;
  end

  if correction > 1 && ~stats.mv
    warning('Roy''s correction only available if CMI is computed with full multivariate setting. Use, e.g., MVMI(X,Y,W,''multivariateBartlett'',true)\n');
  end

  % Initial effective sample size (remove the order of autoregression)
  
  
  % Bartlett-corrected effective sample size
  mv_correlation = false;
  if correction > 0
%     stat.d_2 = stat.N_o ./ diag(stat.var_r);
    if isfield(stats,'correlation_matrix') && stats.correlation_matrix
      stats.d_2 = stats.N_e;
      mv_correlation = true;
    else
      stats.d_2 = diag(stats.N_e);
    end
  else
    stats.d_2 = repmat(stats.N_o,size(estimate));
  end
  stats.d_2 = stats.d_2 - stats.cs - 2;

  if any(stats.d_2 < 1)
    sum_lt1 = sum(stats.d_2 < 1);
    warning('%d Effective DOF < 1 (%s).\n',...
              sum_lt1, mat2str(stats.d_2(stats.d_2 < 1)));
  end
  
  % Compute p-value
  if stats.cmi
    
    l_rvs = zeros(S,stats.dof);
    if parser.Results.RV == 'F'
      % Monte carlo sample the F/t-distributed random variables
      l_rvs = zeros(S,stats.dof);
      for i = 1:stats.dof
        t_rvs = trnd(stats.d_2(i),[S,1]);
        l_rvs(:,i) = stats.d_2(i) ./ (stats.d_2(i) + t_rvs.^2);
      end
    else
      % Monte carlo sample the beta-distributed random variables
      for i = 1:stats.dof
        l_rvs(:,i) = betarnd(stats.d_2(i)/2,1/2,[S,1]);
      end
    end
    
    dist = prod(l_rvs,2);
    
    estimate = exp(-2*stats.to_cmi(estimate));
    
    % Proportion of surrogates less than statistic
    pval = 1-mean(estimate<=dist);
    
    if nargout > 1
      dist = sort(dist);
    end
  else
    % Sum the partial correlations
    r = estimate;
    nu = stats.d_2;
    
    % Compute t-stat...
    t = r.*sqrt((nu)./(1-r.^2));
    
    % ...and get p-value analytically
    pval = 1-fcdf(t.^2,1,nu);
    
    if numel(nu) == 1 && nargout > 1
      if mv_correlation
        M = length(nu);
        dist = zeros(S,M);
        for i = 1:M
          dist(:,i) = tinv(linspace(0,1,S),nu(i))./sqrt(nu(i));
        end
      else
        dist = tinv(linspace(0,1,S),stats.d_2)./sqrt(stats.d_2);
      end
    end
  end
end
