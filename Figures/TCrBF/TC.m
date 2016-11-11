clear all; close all; clc
%Parameters:
rBB    = 0.01;

rBF_up   = 0.1;
rBF_low  = 0.001;
drBF     = 0.0001;
NrBF     = round((rBF_up - rBF_low)/drBF);
rBF      = linspace(rBF_low, rBF_up, NrBF);

mB     = 7/40;
nB     = 100.0;
aB     = pi * rBB / nB^(1/3);
xi     = pi/sqrt(8.0 * nB * aB );
retneg = mB.^(-2) * nB * aB; 

% k-values:
k_up = 100.0;
dk   = 0.11;
Nk   = round(k_up/dk);
k    = linspace(0, k_up, Nk)';

L = zeros(Nk, Nk);
W = zeros(Nk, 1);
epsilon = 0;

%Temperatures:
Tguess = 0.13;
dT     = 0.00005;
T_C     = zeros(1, NrBF);

rBFcheck = 0;
iTmax = 400;
for irBF = 1:NrBF
    rBFnow  = rBF(NrBF - irBF + 1);
    aBF     = pi * rBFnow / nB^(1/3);
    T_low   = Tguess - iTmax * dT;
    T_up    = Tguess + dT;
    NT      = (T_up - T_low)/dT;
    
    Wfactor = 4.0 / pi^2 * (1.0/mB + mB + 2.0) * nB * aBF^2;
    
    Tcheck   = 0;
    for iT = 1:NT
        T  = T_up - iT * dT;
        mu = 1 + pi^2./12 * T.^2; 
        
        kcheck = 0;
        for j = 1:Nk
            kj      = (j-1)*dk;
            epsilon = kj^2 - mu;
            
            W       = - Wfactor .* log( ((k + kj).^2 + 2.0/(xi^2) )./( (k - kj).^2 + 2.0/(xi^2)) );
            L(:, j) = -1.0./pi .* W .* tanh(epsilon./(2*T)) ./ (2 * epsilon) * dk;
            kcheck = kcheck + 1;
            check = [rBFcheck, Tcheck, kcheck]
        end
        Tcheck = Tcheck + 1;
        
        if eigs(L,1) > 1
            T_C(irBF) = T;
            break;
        end
        
    end
    
    if Tcheck == iTmax
        break;
    end
    
    rBFcheck = rBFcheck + 1;
    check = [rBFcheck, Tcheck, kcheck]
    
    Tguess = T_C(irBF); 
end

rBF = fliplr(rBF);
data = [rBF; T_C];

fileID = fopen('TCrBF.txt','w');
fprintf(fileID,'%12s %12s\n','rBF','TC');
fprintf(fileID,'%12.8f %12.8f\n',data);
fclose(fileID);
 




