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

% Loading statistics from chapter 2
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

% Create a vector for Ursell number
Ur_E = logspace(-2, 2, 50);
% Compute skewness and asymmetry empirically from the Ursell number
[Sk_E,As_E] = empirical_fun(Ur_E); 

% Compare empirical and observations
figure; 
subplot(2,1,1);
semilogx(Ur_E,Sk_E);
hold on; 
semilogx(Ur,Sk,'o');
err = std(Sk_E)*ones(size(Sk_E));
errorbar(Ur_E,Sk_E,err);
% std_dev = std(Sk_E);
% curve1 = Sk_E + std_dev;
% curve2 = Sk_E - std_dev;
% x2 = [Ur_E, fliplr(Ur_E)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g');
title('Empirical and observed asymmetry and skeweness');
legend('Empirical','Observed low tide','Observed mid tide','Observed high tide','Empirical +/- std');
ylabel('Skewness','FontWeight','bold');

subplot(2,1,2);
semilogx(Ur_E,As_E);
hold on; 
semilogx(Ur,As,'o');
err = std(As_E)*ones(size(Sk_E));
errorbar(Ur_E,As_E,err);
legend('Empirical','Observed low tide','Observed mid tide','Observed high tide','Empirical +/- std');
xlabel('Ursell number (log10 scale)','FontWeight','bold');
ylabel('Asymmetry','FontWeight','bold');

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
title('Cross shore evolution of asymmetry and skeweness'); 
legend('Empirical','Observed low tide'); 
ylabel('Asymmetry','FontWeight','bold');

subplot(3,1,2);
plot(x,Sk_BJ_E(:,i));
hold on; 
scatter(position,Sk(:,i),'o');
xlim([4000 5000]);
legend('Empirical','Observed low tide'); 
xlabel('','FontWeight','bold');
ylabel('Skeweness','FontWeight','bold');

subplot(3,1,3);
plot(x,zb);
xlim([4000 5000]);
xlabel('Position from offshore sensor [m]','FontWeight','bold');
ylabel('Depth [m]','FontWeight','bold');
end

%% 7.3 Cross-shore evolution of the orbital velocity
Uw = 1;
T = 10;

%phi -pi/2 and r 0,0.3,0.6
phi = -pi/2;
r = [0, 0.3, 0.6];
for i = 1:length(r)
    [u_orbital1,t_orbital] = waveshape(r(i),phi,Uw,T);
    u_orbital(i,:) = u_orbital1;
    u_orbital1 = [];
end
figure()
subplot(3,1,1);
plot(t_orbital,u_orbital)
title("Wave shape for phi = -pi/2 and r = [0, 0.3, 0.6]")
ylim([-1.3 1.5])
ylabel('orbital velocity [m/s]','FontWeight','bold')
legend(["r = 0" "r = 0.3" "r = 0.6"])


%phi 0 and r 0 0.3 0.6
phi = 0;
r = [0, 0.3, 0.6];
for i = 1:length(r)

    [u_orbital1 t_orbital] = waveshape(r(i),phi,Uw,T);
    u_orbital(i,:) = u_orbital1;
    u_orbital1 = [];
    
end
subplot(3,1,2);
plot(t_orbital,u_orbital)
title("Wave shape for phi = 0 and r = [0, 0.3, 0.6]")
ylim([-1.3 1.5])
ylabel('orbital velocity [m/s]','FontWeight','bold')
legend(["r = 0" "r = 0.3" "r = 0.6"])

%phi -pi/2 -pi/4 0 and r 0.6
phi = [-pi/2 -pi/4 0];
r = 0.6;
for i = 1:length(phi)

    [u_orbital1, t_orbital] = waveshape(r,phi(i),Uw,T);
    u_orbital(i,:) = u_orbital1;
    u_orbital1 = [];
    
end
subplot(3,1,3);
plot(t_orbital,u_orbital)
title("Wave shape for phi = [-pi/2, -pi/4, 0] and r = 0.6")
ylim([-1.3 1.5])
ylabel('orbital velocity [m/s]','FontWeight','bold')
xlabel('time [s]','FontWeight','bold')
legend(["phi = -pi/2" "phi = -pi/4" "phi = 0"])


%% Uw for given h, Hrms and T

Uw = Uw_fun(hm(:,1),Hrms_tot(:,1),T); %We are just analysing the low tide so that is why we take the first column oh Hrms and hm

%Predicted time-series of orbital velocity
x = [1000, 4400, 4500, 4700, 4920];
h_low = hm(:,1); %low tide h
k_low = k_fun(T,h_low);
Ur_Elow = Ursell(k_low,h_low,Hrms_tot(:,1));
[Sk_Elow,As_Elow] = empirical_fun(Ur_Elow); 

%And now we compute the orbital velocity for the positions x and the r and
%phi from the empirical function that got the Ursell number from the
%hm data in the first column (lowtide).
u_orb_low = [];
for i = 1:length(h_low) %Number of sensors
    r(i) = computation_r(Sk_Elow(i),As_Elow(i));
    phi(i) = computation_phi(Sk_Elow(i),As_Elow(i));

    [u_orbital1, t_orbital] = waveshape(r(i),phi(i),Uw(i),T);
    u_orb_low(:,i) = u_orbital1;
    u_orbital1 = [];    
end

figure()
plot(t_orbital,u_orb_low)
title("Orbital velocity per position in x")
ylim([-0.7 1.1])
ylabel('orbital velocity [m/s]','FontWeight','bold')
xlabel('time [s]','FontWeight','bold')
legend(["P1" "P3" "P4" "P5" "P6"])


