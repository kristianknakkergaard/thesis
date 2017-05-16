close all; clear all; clc; 

data = load('kdepend');

k       = data(:, 1); 
Delta   = data(:, 2);
mu      = data(:, 3);

epsilon = k.^2 - mu; 
dk      = diff(k); 
E2      = epsilon.^2 + Delta.^2; 


fit1 = fit(k, Delta,'x/(x^2+a)'); 

dDelta  = [diff(Delta)./dk; 0];

plot(k, Delta, 'r-', k, k./(k.^2+fit1.a) , 'b-')
axis([-25 25 -0.4 0.4])


fit1.a

dk = [dk; 0];

CS1 = 1/( 4 .* pi ) .* sum( (epsilon .* dDelta - 2 .* Delta .* k )./E2 .* dk )

Deltafit    = k./(k.^2+fit1.a);
dDeltafit   = - (k.^2-fit1.a)./(k.^2+fit1.a).^2;
E2fit       = epsilon.^2 + Deltafit.^2; 


CS12 = 1/( 4 .* pi ) .* sum( (epsilon .* dDeltafit - 2 .* Deltafit .* k )./E2fit .* dk )