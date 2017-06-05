clc; close all; clear all; 

N_b = 30; 

x_b = rand(1, N_b); 
y1_b = 0.4 + 0.1 * rand(1, N_b);


N_f = 4;

%occupied sites: 
x_f = 1/16 * [1, 3, 7, 11, 15];

N = 1000;
x = linspace(0, 1, N);
y = 0.1 * cos(pi * 16 * x) + 0.5;

y_f = 0.45;

figure(1)
hold on
plot(x, y, '-', 'color', [0, 0, 0] + 0.6)
plot(x_b, y1_b, 'r.', 'MarkerSize', 10)
plot(x_f, y_f, 'b.', 'MarkerSize', 40)


set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'visible','off');

set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis([0, 1, 0.35, 1.65])
hold off
