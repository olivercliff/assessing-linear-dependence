% Runs the numerical evaluation tests from the paper (see below).

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Oliver M. Cliff <oliver.m.cliff@gmail.com>,
%
% If you use this code for your research, please cite the following paper:
%
% Oliver M. Cliff, Leonardo Novelli, Ben D Fulcher, James M. Shine,
% Joseph T. Lizier, "Assessing the significance of directed and multivariate
% measures of linear dependence between time series," Phys. Rev. Research 3,
% 013145 (2021).
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

% Choose settings from the paper (or '' for your own custom configuration)
myfig = '1a';

% This is which experiment (data-point) within the plot/subfigure.
% It increases the autocorrelation or other properties (see paper)
exp = 1;

addpath(genpath('../'));

if isempty(myfig)
  config.T = 2^9; % Dataset length
  config.R = 1000; % Number of runs/trials
  config.S = 1000; % Number of samples for MC distribution
  config.dim_X = 1; % set X dimension
  config.dim_Y = 1; % set Y dimension
  config.dim_W = 0; % set W dimension

  config.to_filter = 1; % None=0, FIR/MA=1, IIR/ARMA=2
  config.filter_order = 24; % FIR/IIR Filter order
  config.ar = true; % Include the AR parameters (or use IID, false)
  config.lm = false;
  config.r = 3.65;

  config.causal = false; % Non-causal (null model, should we be testing true negatives or false positives?)

  config.is_granger = true; % Granger causality (not CMI)
  config.p = 'auto'; % Optimal embedding, set to numeric for a specific embedding
  config.q = 'auto';
  
  config.whiten = false;
  config.whiten_both = true;
  config.arma_p_max = 0;
  config.arma_q_max = 0;

  config.is_pc = false;

  config.alpha = 0.05; % significance level

  config.seed = now; % RNG seed
  
  numericalEvaluation(config);
else
%   numericalEvaluation(myfig,exp,['./results/fig-' myfig '_exp-' num2str(exp) '.mat'],false,false);
  numericalEvaluation(myfig,exp,[],true,false);
end