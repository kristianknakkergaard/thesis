clc; close all; clear all; 

x_b = rand(1, 10000); 
y_b = rand(1, 10000); 
z_b = rand(1, 10000);

x1_f = rand(1, 100); 
%x2_f = rand(1, 100); 
y1_f = 1/2; 
%y2_f = 2/3;
z_f  = 1/2; 

figure(1)
hold on
plot3(x_b, y_b, z_b, 'b.', 'MarkerSize', 4)
plot3(x1_f, y1_f, z_f, 'r.', 'MarkerSize', 2 * 40/7)
%plot3(x2_f, y2_f, z_f, 'r.', 'MarkerSize', 2 * 40/7)

box on
ax = gca;
ax.BoxStyle = 'full';

set(gca,'xtick',[])
set(gca,'xticklabel',[])

set(gca,'ytick',[])
set(gca,'yticklabel',[])

set(gca,'ztick',[])
set(gca,'zticklabel',[])

view([2/3, 1/5, 1/7])

hold off
