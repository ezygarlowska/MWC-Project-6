clear; clc; close all;

%% Chapter 7 - Non-linear wave transformation

%% 7.1 Observations of skewness and asymmetry 

%% 7.1.1 Computation from free-surface elevation time-series

% Main script for the computation of the skewness and asymmetry wave evolution 
% for for a given time-series of free-surface elevation. 

% Loading the files LowTide.txt and highTide.txt
lowTide = load('lowTide.txt');
midTide = load('midTide.txt');
highTide = load('highTide.txt');

%loading statistics from chapter 2
data=load('StatisticsEgmond_3.mat','Hm_tot','Hrms_tot','H13_tot','position');
Hrms_tot=data.Hrms_tot;
position=data.position;

% Compute skewness and asymmetry for all timeseries 
[Sk_low,As_low] = Skewness_asymmetry(lowTide);
[Sk_mid,As_mid] = Skewness_asymmetry(midTide);
[Sk_high,As_high] = Skewness_asymmetry(highTide);
Sk= [Sk_low; Sk_mid; Sk_high]';
As= [As_low; As_mid; As_high]';

% Loading h and T to compute the Ursell number 
data = load('MeanWaterDepth.txt');
Tt = [7.58, 6.69, 5.54];
hm = data;
k = k_fun(Tt,hm);

% Ursell number
Ur = Ursell(k,hm,Hrms_tot); 

%% 7.1.2 Comparison with the empirical fit defined by Ruessink et al. (2012)

Ur_E = logspace(-2, 2, 100);
[Sk_E,As_E] = empirical_fun(Ur_E); 

figure; 
subplot(2,1,1);
semilogx(Ur_E,Sk_E);
hold on; 
semilogx(Ur,Sk,'o');

subplot(2,1,2);
semilogx(Ur_E,As_E);
hold on; 
semilogx(Ur,As,'o');

%% 7.2 Modelling Sk and As

%Load bed profile 
prof=load('prof1018.txt');

% Definition of the array profile, input argument for BJmodel
x = prof(:,1);  
zb = prof(:,2); 
profile = [x zb];
Nx= length(x); 

% Loading data from the BJ model 
waves=load('waves.mat');
eta=waves.waves.eta; 
ht=waves.waves.ht; 
Hrms=waves.waves.Hrms; 

% Ursell number
k = k_fun(Tt,ht);
Ur_BJ = Ursell(k,ht,Hrms); 

% Empirical Sk and As
[Sk_BJ_E,As_BJ_E] = empirical_fun(Ur_BJ); 

for i=1:3
figure; 
subplot(3,1,1);
plot(x,As_BJ_E(:,i));
hold on; 
scatter(position,As(:,i),'o');
xlim([4000 5000]);

subplot(3,1,2);
plot(x,Sk_BJ_E(:,i));
hold on; 
scatter(position,Sk(:,i),'o');
xlim([4000 5000]);

subplot(3,1,3);
plot(x,zb);
xlim([4000 5000]);
end


