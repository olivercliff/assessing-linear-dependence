function [pval,dist] = significance_mv(estimate,stats,varargin)
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

if ~strcmp(parser.Results.test,'modified')) || strcmp(parser.Results.test,'asymptotic')
  t = estimate.*sqrt((stats.N_o-2)./(1-estimate.^2));

  % Get p-value from quantile function of F-dist
  pval = 1-fcdf(t^2,1,stats.N_o-2);

  if nargout > 1
    dist = tinv(linspace(0,1,S),stats.N_o-2)./sqrt(stats.N_o-2);
  end
elseif strcmp(parser.Results.test,'modified')
end