% Test to use the Bartlett correction for inferring significance of
% cross-correlations

% Select to either run the test without a causal influence (causal=false)
% or with a causal influence between X and Y at a specified lag
% (causal=true). The parametres "ar" and "to_filter" will introduce higher
% levels of autocorrelation (see numerical_evaluation.m)

%% Configure

causal = false; % Causal influence from Y to X
lag = 10; % Which lag is this influence at?

ar = true; % Are processes autoregressive?
to_filter = true; % Are they filtered (increase autocorrelation)

T = 1000; % Number of samples
M = ceil(T/5); % Max lag to compute xcorr for (ensure M < T-3)

alpha = 0.05; % Significance level

if ~causal
  out = input('Run null test (no causality)? y/n [y]? ','s');
  if out == 'n'
    return
  end
else
  out = input(sprintf('Run alternative test (causal influence at lag %d)? y/n [y]? ',lag),'s');
  if out == 'n'
    return
  end
end

%% Setup

% Is the original signal autoregressive or spectrally white?
if ar
  phi_X = 0.4;
  phi_Y = -0.8;
else
  phi_X = 0;
  phi_Y = 0;
end

% Should we include a causal influence from Y to X?
if causal
  phi_XY = 0.25 .* eye(dim_X,dim_Y);
else
  phi_XY = 0;
end

% Autoregression matrices
auto_Phi = [phi_X, 0;
            0, phi_Y];
          
lag_Phi = [0, phi_XY;
          0, 0];

% Innovation covariance
m = 2;
Sigma = eye(2);

% Simulate (AR) process
Z = Sigma*randn(2,T);
for t = 2:T
  Z(:,t) = auto_Phi*Z(:,t-1) + Z(:,t);
  if t > lag
    Z(:,t) = lag_Phi*Z(:,t-lag) + Z(:,t);
  end
end

% Partition the dataset (Z) into X and Y
X = Z(1,:)';
Y = Z(2,:)';

% Filter AR process
if to_filter
  [a_coeff, b_coeff] = butter(20, 0.5);
  X = filter(a_coeff,b_coeff,X,[],1);
  Y = filter(a_coeff,b_coeff,Y,[],1);
end

% Setup our outputs
C = zeros(2*M-1,1);
lags = zeros(2*M-1,1);
pval = zeros(2*M-1,1);
quants = zeros(2*M-1,2);

% Setup Matlab's outputs
Cm = zeros(2*M-1,1);
pvalm = zeros(2*M-1,1);
quantsm = zeros(2*M-1,2);

%% Compute cross-correlations

S = 1000;

lag_seq = 1:2*M-1;
for u = lag_seq
  
  m = u-M;
  if m > 0
    xids = (1:T-m)+m;
    yids = (1:T-m);
  else
    yids = (1:T+m)-m;
    xids = (1:T+m);
  end
  
  cX = X(xids,:);
  cY = Y(yids,:);
  
  % Get our correlation implementation
  [C(u),pval(u),cdist,stat] = pcorr(cX,cY,'surrogates',S);
  lags(u) = m;
  quants(u,1) = cdist(round(0.975*S));
  quants(u,2) = cdist(round(0.025*S));
  
  [~,cdistm] = significance(C(u),stat,'surrogates',S,'test','standard');
  quantsm(u,1) = cdistm(round(0.975*S));
  quantsm(u,2) = cdistm(round(0.025*S));
  
  % Get MATLAB's correlation implementation
  [Cm(u),pvalm(u)] = corr(cX,cY);
  
  fprintf('[%i/%i] lag %d (%d), r=%.2f, pval=%.2f\n', u, length(lag_seq), m, length(cX), C(u), pval(u));
end

% Get MATLAB's cross-correlation implementation (normalized to coeff), to
% confirm
[xCm,lagsm] = xcorr(X,Y,M-1,'coeff');

%% Plot results

figure;
hold on;
ph1 = plot(lagsm,xCm,'k--');
ph2 = plot(lags,Cm,'b--');
ph3 = plot(lags,C,'k-','linewidth',2);

if causal
  ph = plot(lag.*[1 1],get(gca,'ylim'),'r-','linewidth',1); uistack(ph,'bottom');
end
phm = plot(lags,quantsm(:,1),'r--','linewidth',1); uistack(phm,'bottom');
phm = plot(lags,quantsm(:,2),'r--','linewidth',1); uistack(phm,'bottom');
ph = plot(lags,quants(:,1),'r-','linewidth',1); uistack(ph,'bottom');
ph = plot(lags,quants(:,2),'r-','linewidth',1); uistack(ph,'bottom');
axis tight;

xlabel('Lag');
ylabel('Corr');
title('Cross-correlations');
legend([ph1 ph2 ph3 ph phm],'Matlab xcorr','Matlab corr','Our corr','threshold exact','treshold standard');

figure;
hold on;
ph1 = plot(lags,pvalm,'b--');
ph2 = plot(lags,pval,'k-','linewidth',2);

set(gca,'yscale','log');

if causal
  ph = plot(lag.*[1 1],get(gca,'ylim'),'r-','linewidth',1); uistack(ph,'bottom');
end
ph = plot(get(gca,'xlim'), alpha.*ones(1,2),'r-','linewidth',1); uistack(ph,'bottom');
axis tight;

xlabel('Lag');
ylabel('P-value');
title('P-values');
legend([ph1 ph2 ph],'t-test','exact','threshold');

if ~causal
  fprintf('FPR for exact t-test: %.2f\n', mean(pval<=alpha)); 
  fprintf('FPR for standard t-test: %.2f\n', mean(pvalm<=alpha)); 
end