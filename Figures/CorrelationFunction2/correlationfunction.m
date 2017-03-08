clear all; close all; clc;
data = load('corfunc1');

x1 = data(:,1);
corfunc1 = data(:,2);

%corfit = fit(x, corfunc, 'a*besselj(1,b*x)*exp(-c*sqrt(abs(x)))');

%corfuncfit = corfit.a.*besselj(1,corfit.b*x).*exp(-corfit.c*sqrt(abs(x)));

%figure
%hold on
%plot(x, corfunc, 'r-')
%plot(x, corfuncfit, 'k-')

[peaks, locs] = findpeaks(abs(corfunc1));

xpeaks = x1(locs);

[peaksfit1, gof1] = fit(xpeaks(xpeaks > 0), peaks(xpeaks > 0), 'a*exp(-(x/b)^c )','StartPoint',[1, 1, 0.5]);
%[logpeaksfit1, gof1] = fit(xpeaks(xpeaks > 0), log(peaks(xpeaks > 0)), 'a + b * x^c','StartPoint',[1, -1, 1]);

peaksfit1function = peaksfit1.a * exp(- (x1(x1 > 0)./peaksfit1.b).^peaksfit1.c);
%logpeaksfit1function = logpeaksfit1.a + logpeaksfit1.b * x(x > 0).^logpeaksfit1.c;

figure
hold on
plot(xpeaks(xpeaks > 0), peaks(xpeaks > 0), 'k.')
%plot(x(x > 0), logpeaksfit1function, 'r-')
plot(x1(x1 > 0), peaksfit1function, 'r-')
axis([0 50 0 0.1])


data3 = load('corfunc3');
x3 = data3(:,1);
corfunc3 = data3(:,2);

%corfit = fit(x, corfunc, 'a*besselj(1,b*x)*exp(-c*sqrt(abs(x)))');

%corfuncfit = corfit.a.*besselj(1,corfit.b*x).*exp(-corfit.c*sqrt(abs(x)));

%figure
%hold on
%plot(x, corfunc, 'r-')
%plot(x, corfuncfit, 'k-')

[peaks3, locs3] = findpeaks(abs(corfunc3));

xpeaks3 = x3(locs3);

[peaksfit3, gof3] = fit(xpeaks3(xpeaks3 > 0), peaks3(xpeaks3 > 0), 'a*exp(-(x/b)^c )','StartPoint',[1, 1, 0.5]);
%[logpeaksfit1, gof1] = fit(xpeaks(xpeaks > 0), log(peaks(xpeaks > 0)), 'a + b * x^c','StartPoint',[1, -1, 1]);

peaksfit3function = peaksfit3.a * exp(- (x3(x3 > 0)./peaksfit3.b).^peaksfit3.c);
%logpeaksfit1function = logpeaksfit1.a + logpeaksfit1.b * x(x > 0).^logpeaksfit1.c;

figure
hold on
plot(xpeaks3(xpeaks3 > 0), peaks3(xpeaks3 > 0), 'k.')
%plot(x(x > 0), logpeaksfit1function, 'r-')
plot(x3(x3 > 0), peaksfit3function, 'r-')
axis([0 50 0 0.1])