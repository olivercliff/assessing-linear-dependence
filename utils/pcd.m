function [pr,resids,cs] = pcd(X,Y,varargin)
%PCD Partial correlation decomposition for vector AR processes.
%   R = PCD(X,Y) returns the KL-by-1 vector of the partial correlation
%   decomposition between the N-by-K matrix X and the N-by-L matrix
%   Y, where KL=K*L. Columns of X and Y correspond to time series
%   and rows correspond to time indices.
%
%   R = PCD(X,Y,W,...) returns the partial correlation decomposition
%   between X and Y conditioned on the N-by-C matrix W.
%
%   [R,RESIDS,C] = PCD(...) also returns the N-by-KL matrix RESIDS, the
%   residual vectors for each independent partial correlation in R.
%
%   [R,RESIDS,CS] = PCD(...) also returns the KL-by-1 vector CS, the vector
%   specifying the dimension of the multivariate condtiional process
%   going into each computation of R.
%
%   See also <a href="matlab:help mvmi">mvmi</a>, <a href="matlab:help mvgc">mvgc</a>, <a href="matlab:help significance">significance</a>

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

addRequired(parser,'X',@isnumeric);
addRequired(parser,'Y',@isnumeric);
addOptional(parser,'W',[],@isnumeric);

parser = parseParameters(parser,X,Y,varargin{:});

T = size(X,1); % Dataset length
k = size(X,2); % Dimension of X
l = size(Y,2); % Dimension of Y

resids = zeros(T,k*l,2);
pr = zeros(k*l,1); % Partial correlations
cs = zeros(k*l,1); % Number of conditionals

% Concomittant variable
W = parser.Results.W;

for j = 1:k
  X_j = X(:,j);
  for i = 1:l
    ij = ((j-1)*l)+i;
    
    Y_i = Y(:,i);
    
    XjYi = [X_j Y_i];
    C_ij = [ones(T,1), W, Y(:,1:i-1), X(:,1:j-1)];

    % Compute and store residuals
    eXYjZ = XjYi - C_ij*(C_ij\XjYi);
    resids(:,ij,:) = eXYjZ;

    % Compute parcorr from a correlation b/w residuals
    pr(ij) = corr(eXYjZ(:,1),eXYjZ(:,2));
    
    % Number of conditions (exclude intercept)
    cs(ij) = size(C_ij,2)-1;
  end
end

if parser.Results.verifyStationary
  acf = autocorr(resids(:,end,1),T-1);
  if mean(acf(2:end) > 1.96/sqrt(T)) > 0.05
    warning('Residuals for X process are non-stationary');
  end
  acf = autocorr(resids(:,end,2),T-1);
  if mean(acf(2:end) > 1.96/sqrt(T)) > 0.05
    warning('Residuals for Y process are non-stationary');
  end
end