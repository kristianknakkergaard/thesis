clear all; close all; clc
%Parameters:
rBB_up   = 0.03;
rBB_low  = 0.004;
drBB     = 0.00002;
NrBB     = round((rBB_up - rBB_low)/drBB);
rBB      = linspace(rBB_low, rBB_up, NrBB);

rBF     = 0.1;
mB      = 7/40;
nB      = 100.0;
aB      = pi .* rBB / nB^(1/3);
aBF     = pi * rBF / nB^(1/3);
Wfactor = 4.0 / pi^2 * (1.0/mB + mB + 2.0) * nB * aBF^2;


retneg = mB.^(-2) * nB .* aB; 

% k-values:
k_up = 100.0;
dk   = 0.11;
Nk   = round(k_up/dk);
k    = linspace(0, k_up, Nk)';

L = zeros(Nk, Nk);
W = zeros(Nk, 1);
epsilon = 0;

%Temperatures:
Tguess = 0.276;
dT     = 0.00008;
T_C     = zeros(1, NrBB);

rBBcheck = 0;
iTmax = 400;
for irBB = 1:NrBB
    rBBnow  = rBB(irBB);
    aBnow   = pi .* rBBnow / nB^(1/3);
    xi      = pi/sqrt(8.0 * nB * aBnow );
    
    T_low   = Tguess - iTmax * dT;
    T_up    = Tguess + dT;
    NT      = (T_up - T_low)/dT;
    Tcheck  = 0;
    
    for iT = 1:NT
        T  = T_up - iT * dT;
        mu = 1 + pi^2./12 * T.^2; 
        
        for j = 1:Nk
            kj      = (j-1)*dk;
            epsilon = kj^2 - mu;
         
            W       = - Wfactor .* log( ((k + kj).^2 + 2.0/(xi^2) )./( (k - kj).^2 + 2.0/(xi^2)) );
            L(:, j) = - 1.0./pi .* W .* tanh(epsilon./(2*T)) ./ (2 * epsilon) * dk;
           
        end
        Tcheck = Tcheck + 1;
        check = [rBBcheck, Tcheck, eigs(L,1)]
        
        if eigs(L,1) > 1
            T_C(irBB) = T;
            break;
        end
        
    end
    
    if Tcheck == iTmax
        break;
    end
    
    rBBcheck = rBBcheck + 1;
    check = [rBBcheck, Tcheck, eigs(L,1)]
    
    Tguess = T_C(irBB); 
end

data = [rBB; T_C; retneg];

fileID = fopen('TCrBB.txt','w');
fprintf(fileID,'%12s %12s %12s\n','rBB','TC', 'retneg');
fprintf(fileID,'%12.8f %12.8f %12.8f \n',data);
fclose(fileID);
 




