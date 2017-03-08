clear all; close all; clc;
data = load('corfunc1');

x = data(:,1);
corfunc = data(:,2);

%corfit = fit(x, corfunc, 'a*besselj(1,b*x)*exp(-c*sqrt(abs(x)))');

%corfuncfit = corfit.a.*besselj(1,corfit.b*x).*exp(-corfit.c*sqrt(abs(x)));

figure(1)
hold on
plot(x, corfunc, 'r-')
%plot(x, corfuncfit, 'k-')

[peaks, locs] = findpeaks(corfunc);

xpeaks = x(locs);
diff(xpeaks)