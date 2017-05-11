clc; close all; clear all; 
fill = linspace(0.01, 1.0, 40);
xi   = linspace(1.0, 21.0, 80);
CS1  = load('xifilldepend');

map = [1.0, 1.0, 1.0
    1.0, 0, 0.0
    0, 0, 1.0];

figure
hold on
xlabel('filling')
ylabel('xi / a')
title('t2 / t1 = 1.0')
colormap(map)
image(fill, xi, CS1, 'CDataMapping', 'scaled')
axis([0.01, 1.0, 1.0, 21.0])