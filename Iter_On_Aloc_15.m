% Iter_On_Aloc_15.m is a routine that solves for consumption in Models 1   
% and 5 by performing fixed-point iteration-on-allocation on consumption of 
% country 1; see Maliar, Maliar and Judd, (2011), "Solving the Multi-Country 
% Real Business Cycle Model Using Ergodic Set Methods", Journal of Economic 
% Dynamics and Control 35, 207-228 (henceforth, MMJ, 2011). 
% -------------------------------------------------------------------------
% Inputs:  "c_old_1" is a guess on consumption of country 1; 
%          "k", "a" and "kprime" are current-period capital, current-period
%          productivity levels and next-period capital, respectively;  
%          "N" is the number of countries; 
%          "alpha", "gam", "phi", "A", "tau" are the parameters of the
%          model;
%          "cdamp" and "citer" are the parameters in iteration-on-allocation 
%
% Output:  "c" is individual consumption of N countries
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [c]  = Iter_On_Aloc_15(c_old_1,k,a,kprime,N,alpha,gam,phi,A,tau,cdamp,citer)

n_rows=size(k,1);             % n_rows is the number of rows in k and a
c = zeros(n_rows,N);          % allocate memory for consumption allocations

% 1. Aggregate consumption, C
%----------------------------
C = (A*k.^alpha.*a - kprime+k-phi/2*(kprime./k-1).^2.*k)*ones(N,1);
% C is computed by summing up individual consumption, which in turn, is 
% found from the individual budget constraints; n_rows-by-1

% 2. Individual consumption if the countries are identical in "gamma" 
%--------------------------------------------------------------------
if (gam - ones(1,N)*mean(gam)) == 0;   % If the countries are identical in 
                                       % "gamma", ... 
    c = C*ones(1,N)/N;                 % Individual consumption is equal to 
                                       % average consumption

% 3. Individual consumption if the countries differ in "gamma": iteration-
% on-allocation
%--------------------------------------------------------------------------
else     
    c(:,1) = c_old_1(:,1);             % A guess on consumption of country 1;
                                       % n_rows-by-1
    
    for n = 1:citer;                   % Number of iterations for iteration-
                                       % on-allocation 

    % Given consumption of country 1, find consumption of the other
    % countries using condition (12) in MMJ (2011)
    for j = 2:N   
        c(:,j) = (tau(1)/tau(j))^(-gam(j))*c(:,1).^(gam(j)/gam(1));     
    end % n_rows-by-N
    
    % Recompute consumption of country 1 using condition (13) in MMJ (2011) 
    c_tilda_1 = C-sum(c(:,2:N),2);     % n_rows-by-1
    
    % Restrict consumption of country 1 to be in the interval [0.5A, 1.5A] 
    % where A is steady-state consumption
    c_tilda_1 = c_tilda_1.*(c_tilda_1>(0.5*A)).*(c_tilda_1<(1.5*A))+0.5*A*(c_tilda_1<=(0.5*A))+1.5*A*(c_tilda_1>=(1.5*A));

    % Update consumption of country 1 for the subsequent iteration using 
    % damping 
    c(:,1) = c(:,1)*(1-cdamp)+c_tilda_1*cdamp;
    
    % Compute the difference between consumption of country 1 at the end and
    % the beginning of the iteration (the convergence criterion is adjusted 
    % to the damping parameter as described in MMJ, 2011)
    dif = mean(mean(abs(1-c_tilda_1./c(:,1))))/cdamp;

    % Stop iterating if consumption of country 1 is computed with an average  
    % absolute error of 1e-10    
    if dif<1e-10; break; end
    end
end