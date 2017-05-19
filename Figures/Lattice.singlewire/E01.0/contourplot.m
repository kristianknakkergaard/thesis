clc; close all; clear all; 
xifill  = load('xifilldepend');
xifill2 = load('xifilldepend2');

E0diff   = xifill(:,5);
E02      = xifill2(:,4);

E0diffpercentage = -E0diff ./ E02;

Nxi   = 61;
Nfill = 101;

xi_low   = 1.0;
dxi      = 0.1;
fill_low = 0.0;
dfill    = 0.01;

fileID = fopen('E0diffpercentagexifilldepend','w');
for j = 1 : Nxi
    xi = xi_low + dxi * (j-1);
    for i = 1: Nfill
        fill                = fill_low + dfill * (i - 1);
        CS11                = xifill((j-1)*Nfill + i, 3);
        CS12                = xifill2((j-1)*Nfill + i, 3);
        if E0diff((j-1)*Nfill + i) < 0.001
            fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f \n', fill, xi, CS11, CS12, 0);
        else
            fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f \n', fill, xi, CS11, CS12, E0diffpercentage((j-1)*Nfill + i));
        end
            
        
    end
    fprintf(fileID, '\n');
end
   
E0diffcal = load('E0diffpercentagexifilldepend');
max(E0diffcal(:,5))
