function [Xf,Yp,Xp,Wp] = embed(X,Y,p,q,W)
  % Embed the Y and X vectors

  N = size(X,1);
  dim_X = size(X,2);
  dim_Y = size(Y,2);

  Xf = X(p+1:end,:);
  Xp = zeros(N-p,dim_X*p);
  for i = 1:p
    seq = i:N-p+i-1;
    dims = (i-1)*dim_X+1:i*dim_X;
    Xp(:,dims) = X(seq,:);
  end

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

  % Ensure W is the right length
  Wp = W(pq_max+1:end,:);
end