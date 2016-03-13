% Simulation_15.m is a routine that simulates a time-series solution for  
% Models 1 and 5: it first computes the capital series using the given 
% capital policy functions, and it then computes the consumption series 
% using fixed-point iteration-on-allocation; see Maliar, Maliar and Judd, 
% (2011), "Solving the Multi-Country Real Business Cycle Model Using Ergodic 
% Set Methods", Journal of Economic Dynamics and Control 35, 207-228 
% (henceforth, MMJ, 2011). 
% -------------------------------------------------------------------------
% Inputs:    "a" is the given vector of current-period productivity levels;
%            "vk" and "vc" are, respectively, the coefficients of the 
%            capital policy functions of N countries and the consumption 
%            policy function of country 1;
%            "alpha", "gam", "phi", "A", "tau" are the parameters of the
%            model
%
% Output:    "c", "k" are the capital and consumption series of N countries
%            that have the same length as "a"
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [c k]  = Simulation_15(vk,vc,a,k0_init,alpha,gam,phi,A,tau)

% 1. Compute the capital series using the given capital policy functions
%-----------------------------------------------------------------------
[n_rows,N] = size(a);             % Infer the simulation length and number of 
                                  % countries 

k = ones(n_rows+1,N);             % Allocate memory for the capital series                                
k(1,1:N) = k0_init;               % Initial condition for capital
pol_bases = ones(n_rows,1+2*N+N*(2*N+1)); 
                                  % Allocate memory for the polynomial bases

for t = 1:n_rows
    z = [k(t,:) a(t,:)];          % Current-period state variables 
    pol_bases(t,:) = Polynomial_2d(z); 
                                  % Construct bases of the complete second-
                                  % degree polynomial 
    k(t+1,1:N) = pol_bases(t,:)*vk;     % Next-period capital 
end

% 2. Compute the consumption series using iteration-on-allocation
% ---------------------------------------------------------------
% 2.1. Parameters in iteration-on-allocation
% ------------------------------------------
cdamp    = 0.01;       % Damping parameter in iteration-on-allocation
citer    = 1e+6;       % Number of iterations on consumption allocations in
                       % iteration-on-allocation; under a very large citer 
                       % (like here), iteration-on-allocation is in effect 
                       % an infinite loop, which breaks when the convergence 
                       % in the consumption allocations with a given accuracy 
                       % of 1e-10 average error is achieved 

% 2.2 Inputs for iteration-on-allocation
%---------------------------------------
kprime = k(2:n_rows+1,1:N);   % Series of next-period period capital
k      = k(1:n_rows,1:N);     % Series of current-period capital

% 2.3 Aggregate consumption, C
%-----------------------------
C = (A*k.^alpha.*a - kprime+k-phi/2*(kprime./k-1).^2.*k)*ones(N,1);
% C is computed by summing up individual consumption, which in turn, is 
% found from the individual budget constraints; n_rows-by-1

% 2.4 Individual consumption if the countries are identical in "gamma" 
%--------------------------------------------------------------------
if (gam - ones(1,N)*mean(gam)) == 0;   % If the countries are identical in 
                                       % "gamma", ... 
    c = C*ones(1,N)/N;                 % Individual consumption is equal to 
                                       % average consumption

% 2.5 Individual consumption if the countries differ in "gamma": iteration-
% on-allocation
%--------------------------------------------------------------------------
else 
    
    c(:,1) = pol_bases*vc;             % A guess on consumption of country 1;
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