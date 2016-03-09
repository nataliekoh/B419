function [out] = getPCvariables(header, coeff, N)
% This function returns a table of strings containing the variables that
% are most heavily weighted in each PC, along with their PC loadings. 

% Inputs: header is a row vector containing variable names as strings;
% coeff is a matrix output from PCA that contains PC loadings;
% N is the number of PCs you want to get variable weights from 
    trans_header = header';
    tmp = zeros(N, 1);
    tmp_coeff = zeros(N, 1);
    for i = 1:N,
        [row, ~] = find(coeff(:,i) == max(coeff(:,i)));
        tmp(i,1) = row;
        tmp_coeff(i,1) = max(coeff(:,i));
        clear row;
    end
    header_names = {};
    for j = 1:N,
        header_names(j) = header(tmp(j,1));
    end
    out = table; 
    out.Variable = header_names';
    out.Variable_Col = tmp;
    out.PC_Loading = tmp_coeff;
end
