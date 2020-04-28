function [coeff,pval,dist,stats] = pcorr(X,Y,varargin)
%PC Correlation (partial or Pearson) for vector AR processes.
%   PR = PCORR(X,Y) returns the scalar estimate of the correlation
%   between the N-by-1 vectors X and Y.
%
%   PR = PCORR(X,Y,W,...) returns the scalar estimate of partial
%   correlation between X and Y conditioned on the N-by-C matrix W.
%
%   [PR,PVAL] = PCORR(...) also returns PVAL, the p-value for testing the
%   hypothesis of no correlation
%
%   [...] = PCORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%
%         Parameter                   Value
%          'test'                     'exact' (the default) uses a
%                                     Bartlett-corrected Student's t-test,
%                                     'standard' uses the typical two-tail
%                                     t-test.
%          'surrogates'               Numeric denoting the number of
%                                     surrogates used in generating the
%                                     exact null distributions.
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
%          'multivariateBartlett'     False (default) to assume all pairs
%                                     of correlations are independent, and
%                                     true to Bartlett correct for full
%                                     covariance matrix.
%
%   Example:
%     % Compute the sample correlation between X and Y and obtain
%     % both Student's t-test p-value and the exact test p-value.
%     % (these values should be similar)
%     X = randn(100,1);
%     Y = randn(100,1);
%     [PR,PVAL] = PCORR(X,Y)
%     [PR,PVAL] = PCORR(X,Y,'test','standard')

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

parser = inputParser;

addRequired(parser,'X',@isvector);
addRequired(parser,'Y',@isvector);
addOptional(parser,'W',[],@(x) (isnumeric(x) && ismatrix(x)));

parser = parseParameters(parser,X,Y,varargin{:});

params = varargin;
if ~contains(parser.UsingDefaults,'W')
  params = varargin(2:end);
end

[pr,eta,cs] = pcd(X,Y,parser.Results.W,params{:});

% Outputs for computing the significance (variance estimation and number of
% condtiionals)
stats.eta = eta;
stats.cs = cs;
stats.mv = parser.Results.multivariateBartlett;
stats.dof = length(cs);
stats.N_o = length(X);
stats.cmi = false;

coeff = sum(pr);

if nargout > 1
  
  sig = @(coeff,stats) significance(coeff,stats,params{:});
  
  if nargout > 2
    [pval,dist] = sig(coeff,stats);
  else
    pval = sig(coeff,stats);
  end
end