% JIDT_TEST Confirms that the measure and LR test match the output of JIDT (a well-known
% toolkit for infodynamics).

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

clear
close all

if ~exist('granger_causality.m','file')
  addpath('..');
end

%% Set up some parameters
% Choose which measure we want to check:
test_options = {'mi','cmi','mvmi','gc','mvgc'};
which_test = 'mvgc'; % Or use, e.g., which_test = test_options{4};

test_id = find(strcmp(test_options,which_test));
if isempty(test_id)
  error('Choose one of the following tests: {%s, %s, %s, %s, %s}\n', test_options{:});
end

% Choose the sample length and RNG (if wanted)
T = 1000; % Number of samples
seed = now; % RNG seed

%% Set up the calculators
javaaddpath('../jidt/infodynamics.jar');  
switch test_id
  case {1,3}
    calc = javaObject('infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian');
  case 2
    calc = javaObject('infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian');
  case {4,5}
    calc = javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorMultiVariateGaussian');
end

if test_id < 4
  compute_measure = @(X,Y,W) mvmi(X,Y,W,'none');
else
  compute_measure = @(X,Y,W) mvgc(X,Y,W,[-1 -1],'none');
end

%% Set up variables and generate random samples

switch test_id
  case {1,4}
    dim_X = 1;
    dim_Y = 1;
    dim_W = 0;
  case 2
    dim_X = 1;
    dim_Y = 1;
    dim_W = 1;
  case {3,5}
    dim_X = 2;
    dim_Y = 2;
    dim_W = 0;
end

% Dimension of Z
D = dim_X + dim_Y + dim_W;

% Partitions
p_X = 1:dim_X;
p_Y = dim_X+1:dim_X+dim_Y;
p_W = dim_X+dim_Y+1:D;

% Seed the RNG
rng(seed);

% Generate random samples
Z = randn(D,T);

X = Z(p_X,:)';
Y = Z(p_Y,:)';
W = Z(p_W,:)';

%% Obtain measurements and p-values

% From our code..
[measure,stats] = compute_measure(X,Y,W);
pval = significance(measure,stats,'lr');

% ...and JIDT
switch test_id
  case {1,3}
    calc.initialise(dim_X,dim_Y)
    calc.setObservations(X,Y);
  case 2
    calc.initialise(dim_X,dim_Y,dim_W);
    calc.setObservations(X,Y,W);
  case {4,5}
    calc.setProperty('k_HISTORY', num2str(stats.p));
    calc.setProperty('l_HISTORY', num2str(stats.q));
    calc.initialise(dim_Y,dim_X);
    calc.setObservations(Y,X);
end
  
measure_JIDT = calc.computeAverageLocalOfObservations();
pval_JIDT = calc.computeSignificance.pValue;

% If testing Granger, multiple TE by 2
if test_id > 3
  measure_JIDT = 2*measure_JIDT;
end

%% Print results

fprintf('Measure from JIDT: %.5g (p = %.5g)\n', measure_JIDT, pval_JIDT);
fprintf('Measure from this code: %.5g (p = %.5g)\n', measure, pval);