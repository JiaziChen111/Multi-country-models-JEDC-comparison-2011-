% Polynomial_2d.m is a routine that constructs the bases of a complete
% second-degree ordinary polynomial
% -------------------------------------------------------------------------
% Inputs:  "z" is the data points on which the polynomial bases must be 
%          constructed; n_rows-by-dimen; 
%
% Output:  "pol_bases" is the bases of the complete polynomial of the given 
%          degree 
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function pol_bases = Polynomial_2d(z)

[n_rows,dimen] = size(z); % Compute the number of rows, n_rows, and number 
                          % of variables (columns), dimen, in matrix z for  
                          % which the polynomial bases must be constructed

    % 1. Bases of the first-degree polynomial (the default option)
    % ------------------------------------------------------------
    pol_bases = [ones(n_rows,1) z];  % these bases include a constant 
                                     % and linear bases

    
for r = 1:n_rows % for each row in z ...
    i = dimen+1; % Index number of a polynomial bases; a constant and linear
                 % bases are indexed from 1 to dimen+1, and higher order 
                 % polynomial bases will be indexed from dimen+2 and on
    
    % 2. Bases of the second-degree polynomial
    % ----------------------------------------            
        for j1 = 1:dimen   
            for j2 = j1:dimen
                i = i+1;
                pol_bases(r,i) = z(r,j1)*z(r,j2);
            end
        end
end                             % end of the loop of the rows
