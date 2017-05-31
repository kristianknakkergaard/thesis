clc; close all; clear all; 
data = load('xifilldependkopi');

E01k = data(:, 5);
E02k = data(:, 6);

Deltamax2k = data(:, 8);

E0diffvector = E02k - E01k; 

winding1k = data(:,3);
winding2k = data(:,4);

winding = zeros(length(E0diffvector),1);

tolerance = 8e-2;

for ixi = 1:length(E0diffvector)
    if Deltamax2k(ixi) < tolerance
        winding(ixi) = winding1k(ixi);
        E0diffvector(ixi) = 0;
    else
        winding(ixi) = winding2k(ixi);
    end  
end

fid = fopen('Deltamaxcorrected.txt', 'w+');

fill_low = 0.0;
dfill    = 0.01;
Nfill    = 101;

xi_low   = 1.0;
dxi      = 0.1;
Nxi      = 61;

for ixi = 1:Nxi
    
    xi = xi_low + (ixi  -1) * dxi;
    
    for jfill = 1:Nfill
        
        fill = fill_low + (jfill - 1) * dfill;
        fprintf(fid, '%4.3f \t %4.3f \t %4.3f \t %4.3f \n', fill, xi, winding(jfill + (ixi - 1)*Nfill ), E0diffvector(jfill + (ixi - 1)*Nfill ) );
    end
    fprintf(fid, '\n');
end
fclose(fid);