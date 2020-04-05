function [pr,nu,c] = compute_partial_correlations(X,Y,W,taper,mv_correction)

if nargin < 4
  taper = 0;
end
if nargin < 5
  mv_correction = false;
end

T = size(X,1); % Dataset length
k = size(X,2); % Dimension of X
l = size(Y,2); % Dimension of Y

all_resids = zeros(T,k*l,2);
pr = zeros(k*l,1); % Partial correlations
c = zeros(k*l,1); % Number of conditionals

for j = 1:k
  X_j = X(:,j);
  for i = 1:l
    ij = ((j-1)*l)+i;
    
    Y_i = Y(:,i);
    C_ij = [W,Y(:,1:i-1),X(:,1:j-1)];

    % Compute and store residuals
    eXjZ = X_j - C_ij*(C_ij\X_j);
    eYiZ = Y_i - C_ij*(C_ij\Y_i);
    
    all_resids(:,ij,1) = eXjZ;
    all_resids(:,ij,2) = eYiZ;

    % Parcorrs are correlation b/w residuals
    pr(ij) = corr(eXjZ,eYiZ);
    
    % Number of conditions
    c(ij) = size(C_ij,2);
  end
end

% Compute Bartlett corrections here so we don't have to pass back the residuals
nu = bartlett_correction_mv(all_resids,taper,mv_correction);