function p = compute_order(X)

T = size(X,1);
D = size(X,2);
M = round(T/3);

ps = zeros(D,1);
for i = 1:D
  cp = find(abs(parcorr(X(:,i),M)) < 1.96./sqrt(T),1);
  if isempty(cp)
    ps(i) = M;
  else
    ps(i) = cp - 2;
  end
end
p = max(ps);

if p == 0
  p = 1;
end