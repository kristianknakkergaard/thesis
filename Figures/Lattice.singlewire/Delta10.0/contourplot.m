clc; close all; clear all; 
fill = linspace(0.01, 1.0, 20);
xi   = linspace(1.0, 21.0, 40);
CS1  = load('xifilldepend');

map = [1.0, 1.0, 1.0
    1.0, 0, 0.0
    0, 0, 1.0];

figure
hold on
xlabel('filling')
ylabel('xi / a')
title('t2 / t1 = 10.0')
colormap(map)
contourf(fill, xi, CS1, [0 1 2])