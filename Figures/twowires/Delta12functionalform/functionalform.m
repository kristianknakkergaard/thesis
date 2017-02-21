clc; close all; clear all; 
data = load('kdepend');

k       = data(:,1);
Delta12 = data(:,3);

fit1 = fit(k(abs(k) >= 1.2), Delta12(abs(k) >= 1.2), 'a*besselk(0, sqrt(b*x^2 + c))','Startpoint',[0.3 0.15, 0.3]);
fit2 = fit(k(abs(k) < 1.2), Delta12(abs(k) < 1.2), 'a - b*x^2' ,               'Startpoint',[1,  1        ]); 

Delta12fit1 = fit1.a * besselk(0, sqrt(fit1.b*k.^2 + fit1.c));
Delta12fit2 = fit2.a - fit2.b * k.^2;

figure(1)
hold on
plot(k, Delta12,'k-')
plot(k(abs(k) >= 1.5), Delta12fit1(abs(k) >= 1.5),'r-')
plot(k(abs(k) < 1.5), Delta12fit2(abs(k) < 1.5), 'g-')
hold off
axis([-25 25 0 0.4])