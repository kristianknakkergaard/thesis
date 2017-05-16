clc; close all; clear all; 
E01k = load('E0depend1k');
E02k = load('E0depend2k');

E0diffmatrix = E02k - E01k;

fid = fopen('E0diff.txt', 'w+');
for i=1:size(E0diffmatrix, 1)
    fprintf(fid, '%8.3f \t', E0diffmatrix(i,:));
    fprintf(fid, '\n');
end
fclose(fid);
[MAX, I_MAX] = max(E0diffmatrix(:));
[MIN, I_MIN] = min(E0diffmatrix(:));

[I_MAX_row, I_MAX_col] = ind2sub(size(E0diffmatrix),I_MAX);
[I_MIN_row, I_MIN_col] = ind2sub(size(E0diffmatrix),I_MIN);

xi_low = 1.0;
dxi = 0.125;
xi_min = xi_low + (I_MIN_row - 1) * 0.125
n_min = (I_MIN_col - 1) * 0.0125


