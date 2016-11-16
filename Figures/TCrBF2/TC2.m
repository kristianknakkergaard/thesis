close all; clear all; clc; 
data = load('TCrBF2.txt');

rBF = data(:,1);
T_C = data(:,2);


TCfit = fit(rBF, T_C, 'a*exp(-b*1/x^2)', 'lower', [0 0], 'upper', [1000 1000]);

TCfitted = TCfit.a * exp(- TCfit.b ./rBF.^2);

figure(1)
hold on
plot(rBF, T_C, 'r-')
plot(rBF, TCfitted, 'k-.')




