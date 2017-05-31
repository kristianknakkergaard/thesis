clc; close all; clear all; 
data = load('xifilldependkopi');

E01k = data(:, 5);
E02k = data(:, 6);

E0diffvector = E02k - E01k; 

winding1k = data(:,3);
winding2k = data(:,4);

winding = zeros(length(E0diffvector),1);

tolerance = 5e-3;
for ixi = 1:length(E0diffvector)
    if E0diffvector(ixi) < 0
        winding(ixi) = winding2k(ixi);
    elseif abs(E0diffvector(ixi)) / abs(E01k(ixi)) < tolerance
        winding(ixi) = winding2k(ixi);
    else
        winding(ixi) = winding1k(ixi);
    end  
end

fid = fopen('E0diffcorrected.txt', 'w+');

fill_low = 0.0;
dfill    = 0.01;
Nfill    = 101;

xi_low   = 4.0;
dxi      = 0.1;
Nxi      = 61;

for ixi = 1:Nxi
    
    xi = xi_low + (ixi  -1) * dxi;
    
    for jfill = 1:Nfill
        
        fill = fill_low + (jfill - 1) * dfill;
        fprintf(fid, '%4.3f \t %4.3f \t %4.3f \n', fill, xi, winding(jfill + (ixi - 1)*Nfill ) );
    end
    fprintf(fid, '\n');
end
fclose(fid);