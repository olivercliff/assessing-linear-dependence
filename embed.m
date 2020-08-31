function [Xf,Xp,Yp,Wf] = embed(X,p,Y,q,W)
%EMBED Embed the past of vector AR process.
%   [Xf,Yp,Xp] = EMBED(X,Y,P,Q) returns the K-by-N matrix Xf of future values of X,
%   the LQ-by-N matrix Yp of past values of Y, and the KP-by-N matrix Xp of past
%   values of X. The inputs are the N-by-L matrix Y, the N-by-K matrix X, and
%   their embedding (history) lengths P and Q as scalars. Columns of X and Y
%   correspond to time series and rows correspond to time indices.
%
%   [...,Wp] = EMBED(...,W) returns the N-by-C matrix Wp of past values of
%   the conditoinal N-by-C matrix W.
%
%   Example:
%     % Fix history lengths p and q to 10, then compute the transfer entropy
%     % from Y to X with this embedding. Confirm this is half the granger causality
%     X = randn(100,5);
%     Y = randn(100,3);
%     p = 10;
%     q = 10;
%     [Xf,Xp,Yp] = embed(X,p,Y,q);
%     te = mvmi(Xf,Yp,Xp)
%     gc = mvgc(X,Y,'p',p,'q',q)
%
%   See also <a href="matlab:help mvgc">mvgc</a>, <a href="matlab:help mvmi">mvmi</a>, <a href="matlab:help order">order</a>

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

  N = size(X,1);
  dim_X = size(X,2);

  Xf = X(p+1:end,:);
  Xp = zeros(N-p,dim_X*p);
  for i = 1:p
    seq = i:N-p+i-1;
    dims = (i-1)*dim_X+1:i*dim_X;
    Xp(:,dims) = X(seq,:);
  end
  
  if nargin > 2
    dim_Y = size(Y,2);
    
    Yp = zeros(N-q,dim_Y*q);
    for i = 1:q
      seq = i:N-q+i-1;
      dims = (i-1)*dim_Y+1:i*dim_Y;
      Yp(:,dims) = Y(seq,:);
    end

    pq_max = max([p,q]);
    Xf = Xf(1-p+pq_max:end,:);
    Xp = Xp(1-p+pq_max:end,:);
    Yp = Yp(1-q+pq_max:end,:);
    
    if nargin > 4
      % Ensure W is the right length
      Wf = W(pq_max+1:end,:);
    end
  end
end