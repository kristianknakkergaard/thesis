clc; close all; clear all; 

N_b = 100; 

x_b = rand(1, N_b); 
y_b = rand(1, N_b); 

N_f = 5;

y1_f = linspace(1/(N_f + 2), 1 - 1/(N_f + 2), N_f) + 0.05 .* rand(1, N_f); 
y2_f = linspace(1/(N_f + 2), 1 - 1/(N_f + 2), N_f) + 0.05 .* rand(1, N_f);

N = 200;
y    = linspace(0, 1, N);

x1_f = 1/3; 
x2_f = 2/3;

figure(1)
hold on
plot(x1_f*ones(1, N), y, '-', 'color', [0, 0, 0] + 0.6)
plot(x2_f*ones(1, N), y, '-', 'color', [0, 0, 0] + 0.6)
plot(x_b, y_b, 'b.', 'MarkerSize', 15)
plot(x1_f, y1_f, 'r.', 'MarkerSize', 40)
plot(x2_f, y2_f, 'r.', 'MarkerSize', 40)


set(gca,'xtick',[])
set(gca,'xticklabel',[])
 set(gca,'visible','off');

set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'visible','off');
axis([0, 1, 0, 1])
hold off
