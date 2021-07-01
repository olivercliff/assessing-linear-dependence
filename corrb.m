function [coeff,pval,stats] = corrb(X,varargin)
% Pearson Correlation for vector AR processes.
%   R = CORRB(X) returns the scalar estimate of the correlation
%   between each pair of columns in the N-by-P matrix X.
%
%   [R,PVAL] = CORRB(...) also returns PVAL, the p-value for testing the
%   hypothesis of no correlation
%
%   [...] = CORRB(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%
%         Parameter                   Value
%          'test'                     'modified' (the default) uses our
%                                     modified t-test,
%                                     'finite' or 'asymptotic'
%                                     uses the two-tailed t-test.
%          'surrogates'               Integer denoting the number of
%                                     surrogates used in generating the
%                                     numerical null distributions.
%          'varianceEstimator'        'bartlett' (default) uses Bartlett's
%                                     formula assuming no
%                                     cross-correlations, 'roy' makes no
%                                     assumptions about cross-correlations.
%          'taperMethod'              'none' (default) to compute
%                                     sample autocorrelations without
%                                     tapering, 'tukey' to use the Tukey
%                                     windowing, 'parzen' for Parzen
%                                     windows, or 'bartlett' to use
%                                     Barttlett's correction. 
%
%   Example:
%     % Compute the sample correlation between X and Y and obtain
%     % both Student's t-test p-value and the exact test p-value.
%     % (these values should be similar)
%     X = randn(100,2);
%     [R,PVAL] = CORRB(X)
%     [R,PVAL] = CORRB(X,'test','modified')

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

addRequired(parser,'X',@ismatrix);

parser = parseParameters(parser,X,varargin{:});

coeff = corr(X);

M = length(coeff);
T = size(X,1);
if nargout > 1
   
  [eta,ess] = bartlett_mv(X,parser.Results.taperMethod);

  % Outputs for computing the significance (variance estimation and number of
  % condtiionals)
  stats.N_e = ess;
  stats.eta = eta;
  stats.cs = zeros(M,1);
  stats.mv = parser.Results.multivariateBartlett;
  stats.dof = length(stats.cs);
  stats.N_o = T-1;
  stats.cmi = false;
  stats.correlation_matrix = true;

  pval = significance(coeff,stats,varargin{:});
end