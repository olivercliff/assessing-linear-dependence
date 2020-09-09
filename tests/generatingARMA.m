addpath ../utils/

runs = 50;
N = 2^9;
max_order = 15;

bart_ar = zeros(runs,max_order);
bart_arma = zeros(runs,max_order);

for p = 1:max_order
  for r = 1:runs
    failed = true;
    while failed
      phi = generateARMAParams(p,0);
      try
        mdl = arima('AR',phi,'constant',0,'variance',1);
        failed = false;
      catch
      end
    end
    X = simulate(mdl,N);
    Y = simulate(mdl,N);
    bart_ar(r,p) = sum(autocorr(X,50).*autocorr(Y,50));
  end
  
  fprintf('Completed order %i/%i\n',p,max_order)
end

figure; plot(max(bart_ar,[],1))
figure; plot(mean(bart_ar,1))