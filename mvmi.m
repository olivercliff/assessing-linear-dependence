function [cmi,stats] = mvmi(X,Y,W,taper_method_str)
%MVMI Mutual information (conditional or unconditional) for vector AR processes.
%   I = MVMI(X,Y) returns the scalar estimate of the linear-Gaussian
%   Mutual information between the N-by-K matrix X and the N-by-L
%   matrix Y. Columns of X and Y correspond to time series and rows
%   correspond to time indices.
%
%   I = MVMI(X,Y,W,...) returns the scalar estimate of linear-Gaussian
%   Mutual information between X and Y conditioned on the N-by-C matrix W.
%
%   [I,STATS] = MVMI(...) also returns STATS, a structure of statistics that
%   are used as input to <a href="matlab:help significance">significance</a> function, e.g.,
%
%       [I,STATS] = MVMI(X,Y);
%       pval = significance(F,STATS);      
%
%   [...] = MVMI(...,TAPER) optional string specifying the taper method used
%   for inferring the autocorrelation function for Bartlett''s formula using one
%   of the following techniques:
%     - 'none' (default)
%     - 'tukey'
%     - 'parzen'
%     - 'bartlett'
%
%   Example:
%     % Compute the linear-Gaussian mutual information between X and Y and obtain
%     % both the LR test p-value and the exact test p-value.
%     X = randn(100,5);
%     Y = randn(100,3);
%     [I,stats] = MVMI(X,Y);
%     pval_LR = significance(I,stats,'lr');
%     pval_exact = significance(I,stats,'exact');
%
%   See also <a href="matlab:help significance">significance</a>, <a href="matlab:help mvgc">mvgc</a>, <a href="matlab:help pcd">pcd</a>

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
  taper_method_str = 'none';
end

[pr,eta,cs] = pcd(X,Y,W,taper_method_str);

% Outputs for computing the significance (variance estimation and number of
% condtiionals)
stats.eta = eta;
stats.cs = cs;
stats.mv = false;
stats.dof = length(cs);
stats.N_o = length(X);
stats.to_cmi = @(x) x;

stats.cmi = true;

cmi = -0.5*sum(log(1-pr.^2));

end