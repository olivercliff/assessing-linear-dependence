function [F,pval,dist,stats] = mvgc(X,Y,varargin)
%MVGC Multivariate Granger causality (conditional or unconditional).
%   F = MVGC(X,Y) returns the scalar estimate of the Granger causality
%   from the N-by-L matrix Y to the N-by-K matrix X. Columns of X and Y
%   correspond to time series and rows correspond to time indices.
%
%   F = MVGC(X,Y,W,...) returns the scalar estimate of Granger causality
%   between X and Y conditioned on the N-by-C matrix W.
%
%   [F,PVAL] = MVGC(...) also returns PVAL, the p-value for testing the
%   hypothesis of no dependence
%
%   [...] = MVGC(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%
%         Parameter                   Value
%          'p'                        'auto' (the default) uses automatic
%                                     embedding, otherwise input the
%                                     desired history length as a string
%          'q'                        'auto' (the default) uses automatic
%                                     embedding, otherwise input the
%                                     desired history length as a string
%          'test'                     'modified' (the default) uses our
%                                     modified lambda-test,
%                                     'finite' uses the F-test, and
%                                     'asymptotic' uses the chi-square
%                                     test.
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
%     % Compute the Granger causality from Y to X and obtain the LR test p-value
%     % and the exact test p-value.
%     X = randn(100,5);
%     Y = randn(100,3);
%     [F,pval_exact] = MVGC(X,Y);
%     [F,pval_LR] = MVGC(X,Y,'test','finite');
%
%   See also <a href="matlab:help order">order</a>, <a href="matlab:help mvmi">mvmi</a>, <a href="matlab:help pcd">pcd</a>

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

isnumericMatrix = @(x) (isnumeric(x) && ismatrix(x));

addRequired(parser,'X',isnumericMatrix);
addRequired(parser,'Y',isnumericMatrix);
addOptional(parser,'W',[],isnumericMatrix);

parser = parseParameters(parser,X,Y,varargin{:});
params = varargin;
if ~contains(parser.UsingDefaults,'W')
  params = varargin(2:end);
end

% Embedding (i.e., history length)
if strcmp(parser.Results.p,'auto')
  p = order(X);
else 
  p = str2double(parser.Results.p);
end

if strcmp(parser.Results.q,'auto')
  q = order(Y);
else 
  q = str2double(parser.Results.q);
end

% Embed the vectors for input to CMI calculator
[Xf,Xp,Yp,Wp] = embed(X,p,Y,q,parser.Results.W);

% Add any conditional matrix
if isempty(Wp)
  XpW = Xp;
else
  XpW = [Xp, Wp];
end

% Calculate CMI (also returning structure for computing significance)
if nargout > 1
  if nargout > 2
    [cmi,pval,dist,stats] = mvmi(Xf,Yp,XpW,params{:});
    stats.p = p;
    stats.q = q;
    stats.to_cmi = @(x) 0.5*x;
  else
    [cmi,pval] = mvmi(Xf,Yp,XpW,params{:});
  end
  
else
  cmi = mvmi(Xf,Yp,XpW,params{:});
end
F = 2*cmi;

end