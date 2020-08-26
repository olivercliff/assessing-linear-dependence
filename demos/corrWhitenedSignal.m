clear
close all

R = 1000;
N = 2^5;

addpath('..');
addpath('../utils/');

orders = [1, 5:5:20];
O = length(orders);

pval_pw_standard = zeros(O,R);
pval_standard = zeros(O,R);
pval_pw_exact = zeros(O,R);
pval_exact = zeros(O,R);

ess_pw = zeros(O,R);
ess = zeros(O,R);

for o = 1:O
  
  order = orders(o);
     
  cpval_pw_standard = zeros(R,1);
  cpval_standard = zeros(R,1);
  cpval_pw_exact = zeros(R,1);
  cpval_exact = zeros(R,1);

  cess = zeros(R,1);
  cess_pw = zeros(R,1);

  parfor r = 1:R
    
    not_set = true;
    fprintf('Attempted coefficients for X''s VAR(%d): ', order);
    while not_set
      not_set = false;
      ar_X = num2cell(rand(order,1).*2-1);
      fprintf('%s, ',mat2str([ar_X{:}],2));
      try
        Mdl_X = arima('Constant',0.0,'AR',ar_X,'Variance',1);
      catch
        not_set = true;
      end
    end
    fprintf('\n[%d/%d X] Successfully selected X''s coefficients: %s.\n',r,R,mat2str([ar_X{:}],2));

    not_set = true;
    fprintf('Attempted coefficients for Y''s VAR(%d): ', order);
    while not_set
      not_set = false;
      ar_Y = num2cell(rand(order,1).*2-1);
      fprintf( '%s, ',mat2str([ar_Y{:}],2));
      try
        Mdl_Y = arima('Constant',0.0,'AR',ar_Y,'Variance',1);
      catch
        not_set = true;
      end
    end
    fprintf('\n[%d/%d Y] Successfully selected Y''s coefficients: %s.\n',r,R,mat2str([ar_Y{:}],2));
    
    X = simulate(Mdl_X,N);
    Y = simulate(Mdl_Y,N);

    X = zscore(detrend(X));
    Y = zscore(detrend(Y));

    [r_XY,cpval_standard(r),~,stat] = pcorr(X,Y,'test','standard','taperMethod','none');
    cess(r) = stat.N_e;

    [X_pw,Y_pw] = prewhiten(X,Y,[],order*2);

    [r_XY_pw,cpval_pw_standard(r),~,stat_pw] = pcorr(X_pw,Y_pw,'test','standard','taperMethod','none');
    cess_pw(r) = stat_pw.N_e;
  end
  
  ess(o,:) = cess;
  ess_pw(o,:) = cess_pw;
    
  pval_pw_standard(o,:) = cpval_pw_standard;
  pval_standard(o,:) = cpval_standard;
end

%% Plots
figure;
errorbar(mean(ess,2), var(ess,[],2) ./ sqrt(R)); xlabel('Order'); ylabel('ESS');
figure;
errorbar(mean(ess_pw,2), var(ess_pw,[],2) ./ sqrt(R)); xlabel('Order'); ylabel('ESS after whitening');

figure;
errorbar(mean(pval_standard <= 0.05,2), var(pval_standard <= 0.05,[],2) ./ sqrt(R)); xlabel('Order'); ylabel('FPR');

figure;
errorbar(mean(pval_pw_standard <= 0.05,2), var(pval_pw_standard <= 0.05,[],2) ./ sqrt(R)); xlabel('Order'); ylabel('FPR after whitening');