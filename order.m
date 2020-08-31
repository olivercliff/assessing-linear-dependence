function p = order(Z,numLags)
%EMBED Embed the past of vector AR process.
%   P = ORDER(X) returns the optimal AR order P of the the N-by-K matrix X
%   using the first partial autocorrelation that is less than 1.96/sqrt(N) or
%   round(N/3), whichever comes first.
%
%   Example:
%     % Infer optimal history lengths p and q, then compute the transfer entropy
%     % from Y to X with this embedding. Confirm this is half the granger causality
%     X = randn(100,5);
%     Y = randn(100,3);
%     p = order(X);
%     q = order(Y);
%     [Xf,Yp,Xp] = embed(X,Y,p,q);
%     te = mvmi(Xf,Yp,Xp)
%     gc = mvgc(X,Y,[p q])
%
%   See also <a href="matlab:help mvgc">mvgc</a>, <a href="matlab:help mvmi">mvmi</a>, <a href="matlab:help embed">embed</a>

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

T = size(Z,1);
D = size(Z,2);

if nargin == 1
  numLags = round(T/3);
end

ps = zeros(D,1);
for i = 1:D
  X = Z(:,i);
  alpha = [1;zeros(numLags,1)];
  
  cp = [];
  for k = 1:numLags
    lagX = lagmatrix(X,1:k+1);
    
    [Q,R] = qr([ones((length(X)-k),1)  lagX(k+1:end,1:k)],0);
    b = R\(Q'*X(k+1:end));
    alpha(k+1) = b(end);
    
    if abs(alpha(k+1)) <= 1.96/sqrt(T-k)
      cp = k;
      break;
    end
  end
  if isempty(cp)
    ps(i) = numLags;
  else
    ps(i) = cp;
  end
end
p = max(ps);

if p == 0
  p = 1;
end