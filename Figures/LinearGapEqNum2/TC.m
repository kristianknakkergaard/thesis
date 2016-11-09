clc; clear all; close all;
data = load('FindTC');

T_C  = data(:,3);
rBF =  data(:,1);

plot(rBF,T_C,'k-')
