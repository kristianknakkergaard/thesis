clc; close all; clear all; 
data = load('kdepend');

k       = data(:,1);
Delta12 = data(:,3);

fit1 = fit(k, Delta12, 'a*exp(-b*x^2)','Startpoint',[0.3 0.15]);

Delta12fit = fit1.a * exp(-fit1.b * k.^2);

figure(1)
hold on
plot(k, Delta12,'k-')
plot(k, Delta12fit,'r-')
hold off