clc; close all; clear all; 
xifill1  = load('xifilldepend1');
xifill2  = load('xifilldepend2');

E01      = xifill1(:,4);
E02      = xifill2(:,4);
E0diff   = E02 - E01;

Nxi   = 61;
Nfill = 101;

xi_low   = 1.0;
dxi      = 0.1;
fill_low = 0.0;
dfill    = 0.01;

fileID = fopen('E0diff','w');
for j = 1 : Nxi
    xi = xi_low + dxi * (j-1);
    for i = 1: Nfill
        fill        = fill_low + dfill * (i - 1);
        CS11        = xifill1((j-1)*Nfill + i, 3);
        CS12        = xifill2((j-1)*Nfill + i, 3);
        E0diffvalue = E0diff((j-1)*Nfill + i);
        fprintf(fileID,'%12.8f \t %12.8f \t %12.8f \t %12.8f \t %12.8f \n', fill, xi, CS11, CS12, E0diffvalue);
    end
    fprintf(fileID, '\n');
end
    


