clear; clc; close all;
% Load input signal
data = load('NP.1875.01.HNE.dat');
% Samplig rate, Fs, as samples-per-second
Fs = 200;
% Detrend input signal
x = detrend(data(:,2));
% Compute whiten signal within frequency rage of 0.1 to 20 Hz.
xnew = whitening(x, Fs, 'freq', [0.1,20]); 
% Plot original and whitened signal 
figure
plot(x,'k'); hold; plot(xnew,'r'); legend('original','prewhitened')