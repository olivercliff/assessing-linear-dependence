function [gc,stats] = granger_causality(X,Y,W,embedding,taper_method_str,mv_bartlett_correction)

% Compute granger causality from predictee Y to predictor X

if nargin < 6
  mv_bartlett_correction = false;
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
  p = compute_order(X);
  
elseif embedding(1) == 0
  %   ==0: 1st-order AR embedding (p=q=1),
  p = 1;
else
  %   >0: set predictee and predictor embedding to value (p=q=embedding),
  p = embedding(1);
end

if embedding(2) < 0
  q = compute_order(Y);
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

% Calculate CMI (returning dist structure for computing significance)
[cmi,stats] = conditional_mutual_information(Xf,Yp,XpW,taper_method_str,mv_bartlett_correction);
gc = 2*cmi;

% How to convert GC to CMI for significance calcs
stats.to_cmi = @(x) 0.5*x;

stats.p = p;
stats.q = q;

end