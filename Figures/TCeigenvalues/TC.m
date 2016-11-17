clear all; close all; clc
%Parameters:
rBB     = 0.01;
rBF     = 0.1;
mB      = 7/40;
nB      = 100.0;
aB      = pi .* rBB / nB^(1/3);
aBF     = pi * rBF / nB^(1/3);
xi      = pi/sqrt(8.0 * nB * aB );
Wfactor = 4.0 / pi^2 * (1.0/mB + mB + 2.0) * nB * aBF^2;


retneg = mB.^(-2) * nB .* aB; 

% k-values:
k_up = 50.0;
dk   = 0.03;
Nk   = round(k_up/dk);
k    = linspace(0, k_up, Nk)';

L = zeros(Nk, Nk);
W = zeros(Nk, Nk);

for j = 1:Nk
    kj     = (j-1)*dk;
    W(:,j) = - Wfactor .* log( ((k + kj).^2 + 2.0/(xi^2) )./( (k - kj).^2 + 2.0/(xi^2)) );
end

epsilon = 0;

%Temperatures:
T_low  = 0.001;
T_up   = 0.2;
dT     = 0.001;
NT = round( (T_up - T_low)/dT );
T = linspace(T_low, T_up, NT);

%eigenvalues:
number_eigen = 5;
eigen = zeros(number_eigen, NT);

Tcheck = 0

for iT = 1:NT
    mu = 1 + pi^2./12 * T(iT)^2;
    
    for j = 1:Nk
        kj      = (j-1)*dk;
        epsilon = kj^2 - mu;
        L(:, j) = - 1.0./pi .* W(:,j) .* tanh(epsilon./(2*T(iT))) ./ (2 * epsilon) * dk;    
    end
    
    Tcheck = Tcheck + 1
    
    eigen(:, iT) = eigs(L,number_eigen)'; 
end

data = [T; eigen];

fileID = fopen('TCeigen.txt','w');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n',data);
fclose(fileID);
 
figure(1)
hold on

for j = 1:number_eigen
    plot(T, eigen(j,:),'m-')
end
hold off



