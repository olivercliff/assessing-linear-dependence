function [phi,theta,r_p,r_q] = generateARMAParams(p,q)

if p > 0
  r_p = zeros(p,1);
  y = zeros(p,p);
  for k = 1:p
    r_p(k) = betarnd(floor(0.5*(k+1)),floor(0.5*k)+1)*2-1;
    for i = 1:k-1
      y(i,k) = y(i,k-1) - r_p(k) * y(k-i,k-1);
    end
    y(k,k) = r_p(k);
  end
  phi = y(:,p);
else
  r_p = [];
  phi = [];
end

if q > 0
  r_q = zeros(q,1);
  y = zeros(q,q);
  for k = 1:q
    r_q(k) = betarnd(floor(0.5*(k+1)),floor(0.5*k)+1)*2-1;
    for i = 1:k-1
      y(i,k) = y(i,k-1) - r_q(k) * y(k-i,k-1);
    end
    y(k,k) = r_q(k);
  end
  theta = y(:,q);
else
  r_q = [];
  theta = [];
end