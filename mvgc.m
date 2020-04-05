function [F,stats] = mvgc(X,Y,W,embedding,taper_method_str)
%MVGC Multivariate Granger causality (conditional or unconditional).
%   F = MVGC(X,Y) returns the scalar estimate of the Granger causality
%   from the N-by-L matrix Y to the N-by-K matrix X. Columns of X and Y
%   correspond to time series and rows correspond to time indices.
%
%   F = MVGC(X,Y,W,...) returns the scalar estimate of Granger causality
%   between X and Y conditioned on the N-by-C matrix W.
%
%   [F,STATS] = MVGC(...) also returns STATS, a structure of statistics that
%   are used as input to <a href="matlab:help significance">significance</a> function.
%
%   [...] = MVGC(...,EMBEDDING) optional 2 element vector for setting the
%   embeddings (history lengths) of both X and Y. If EMBEDDING(i) < 0, the
%   embedding is optimally chosen using the <a href="matlab:help order">order</a> function. If EMBEDDING(i) > 0,
%   then a fixed embedding length of EMBEDDING(i) is used. Default is
%   EMBEDDING = [-1 -1].
%
%   [...] = MVGC(...,TAPER) optional string specifying the taper method used
%   for inferring the autocorrelation function for Bartlett''s formula using one
%   of the following techniques:
%     - 'none' (default)
%     - 'tukey'
%     - 'parzen'
%     - 'bartlett'
%
%   Example:
%     % Compute the Granger causality from Y to X and obtain the LR test p-value
%     % and the exact test p-value.
%     X = randn(100,5);
%     Y = randn(100,3);
%     [F,stats] = MVGC(X,Y);
%     pval_LR = significance(F,stats,'lr');
%     pval_exact = significance(F,stats,'exact');
%
%   See also <a href="matlab:help significance">significance</a>, <a href="matlab:help order">order</a>, <a href="matlab:help mvmi">mvmi</a>, <a href="matlab:help pcd">pcd</a>

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
  W = [];
end

if nargin < 4
  embedding = [-1 -1]; % optimal embedding
end
if nargin < 5
  taper_method_str = 'none';
end

% Embedding (i.e., history length)
if embedding(1) < 0
  %   <0: optimal embedding,
  p = order(X);
  
elseif embedding(1) == 0
  %   ==0: 1st-order AR embedding (p=q=1),
  p = 1;
else
  %   >0: set predictee and predictor embedding to value (p=q=embedding),
  p = embedding(1);
end

if embedding(2) < 0
  q = order(Y);
elseif embedding(2) == 0
  q = 1;
else
  q = embedding(2);
end

% Embed the vectors for input to CMI calculator
[Xf,Yp,Xp,Wp] = embed(X,Y,p,q,W);

% Add any conditional matrix
if isempty(Wp)
  XpW = Xp;
else
  XpW = [Xp, Wp];
end

% Calculate CMI (also returning structure for computing significance)
[cmi,stats] = mvmi(Xf,Yp,XpW,taper_method_str);
F = 2*cmi;

% How to convert GC to CMI for significance calcs
stats.to_cmi = @(x) 0.5*x;

stats.p = p;
stats.q = q;

end