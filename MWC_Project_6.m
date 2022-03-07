clear; clc; close all;

%% Chapter 6 - Non-linear wave transformation

%% 6.1 Observations of skewness and asymmetry  

% Main script for the computation of the skewness and asymmetry wave evolution 
% for for a given time-series of free-surface elevation. 

% Loading the files LowTide.txt and highTide.txt
lowTide = load('lowTide.txt');
midTide = load('midTide.txt');
highTide = load('highTide.txt');
lowP1 = lowTide(:,1);
lowP3 = lowTide(:,2);
lowP6 = lowTide(:,5);

midP1 = midTide(:,1);
midP3 = midTide(:,2);
midP6 = midTide(:,5);

highP1 = highTide(:,1);
highP3 = highTide(:,2);
highP6 = highTide(:,5);


% Compute skewness and asymmetry for all timeseries 
[Sk_low,As_low] = Skewness_asymmetry(lowTide);
[Sk_mid,As_mid] = Skewness_asymmetry(midTide);
[Sk_high,As_high] = Skewness_asymmetry(highTide);