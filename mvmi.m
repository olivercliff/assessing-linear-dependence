function [cmi,pval,dist,stats] = mvmi(X,Y,varargin)
%MVMI Mutual information (conditional or unconditional) for vector AR processes.
%   I = MVMI(X,Y) returns the scalar estimate of the linear-Gaussian
%   Mutual information between the N-by-K matrix X and the N-by-L
%   matrix Y. Columns of X and Y correspond to time series and rows
%   correspond to time indices.
%
%   I = MVMI(X,Y,W,...) returns the scalar estimate of linear-Gaussian
%   Mutual information between X and Y conditioned on the N-by-C matrix W.
%
%   [I,PVAL] = MVMI(...) also returns PVAL, the p-value for testing the
%   hypothesis of no correlation
%
%   [...] = MVMI(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%
%         Parameter                   Value
%          'test'                     'modified' (the default) uses our
%                                     modified lambda-test,
%                                     'finite' uses the F-test, and
%                                     'asymptotic' uses the chi-square
%                                     test
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
%          'multivariateBartlett'     False (default) to assume all pairs
%                                     of correlations are independent, and
%                                     true to Bartlett correct for full
%                                     covariance matrix.
%
%   Example:
%     % Compute mutual information between X and Y and obtain
%     % both Student's t-test p-value and the exact test p-value.
%     % (these values should be similar)
%     X = randn(100,3);
%     Y = randn(100,5);
%     [I,PVAL] = MVMI(X,Y)
%     [I,PVAL] = MVMI(X,Y,'test','finite')
%
%   See also <a href="matlab:help mvgc">mvgc</a>, <a href="matlab:help pcd">pc</a>

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

isnumericMatrix = @(x) (isnumeric(x) && ismatrix(x));

addRequired(parser,'X',isnumericMatrix);
addRequired(parser,'Y',isnumericMatrix);
addOptional(parser,'W',[],isnumericMatrix);

parser = parseParameters(parser,X,Y,varargin{:});
params = varargin;
if ~contains(parser.UsingDefaults,'W')
  params = varargin(2:end);
end

[pr,resids,cs] = pcd(X,Y,parser.Results.W,params{:});

cmi = -0.5*sum(log(1-pr.^2));

if nargout > 1
  
  [var_r,N_e] = bartlett(resids,parser.Results.taperMethod,parser.Results.multivariateBartlett);
  
  % Outputs for computing the significance (variance estimation and number of
  % condtiionals)
  stats.pr = pr;
  stats.var_r = var_r;
  stats.N_e = N_e;
  stats.cs = cs;
  stats.mv = false;
  stats.dof = length(cs);
  stats.N_o = length(X);
  stats.to_cmi = @(x) x;
  stats.cmi = true;
  
  sig = @(cmi,stats) significance(cmi,stats,params{:});
  if nargout > 2
    [pval,dist] = sig(cmi,stats);
  else
    pval = sig(cmi,stats);
  end
end

end