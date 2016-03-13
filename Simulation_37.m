% Simulation_37.m is a routine that simulates a time-series solution for  
% Models 3 and 7: it first computes the capital series using the given 
% capital policy functions, and it then computes the consumption and labor 
% series using fixed-point iteration-on-allocation; see Maliar, Maliar and 
% Judd, (2011), "Solving the Multi-Country Real Business Cycle Model Using 
% Ergodic Set Methods", Journal of Economic Dynamics and Control 35, 207-228 
% (henceforth, MMJ, 2011). 
% -------------------------------------------------------------------------
% Inputs:    "a" is the given vector of current-period productivity levels;
%            "vk" and "vh" are, respectively, the coefficients of the 
%            capital and labor policy functions of N countries;
%            "alpha", "gam", "psi", "Le", "phi", "A", "tau" are the parameters 
%            of the model
%
% Output:    "c", "h", "k" are the capital and consumption series of N 
%            countries that have the same length as "a"
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [c h k]  = Simulation_37(vk,vh,a,k0_init,alpha,gam,psi,Le,phi,A,tau)

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

% 2. Compute the consumption and labor series using iteration-on-allocation
% -------------------------------------------------------------------------
% 2.1. Parameters in iteration-on-allocation
%-------------------------------------------
hdamp    = 0.01;      % Damping parameter in iteration-on-allocation
hiter    = 1e+6;      % Number of iterations on labor allocations in
                      % iteration-on-allocation; under a very large hiter 
                      % (like here), iteration-on-allocation is in effect 
                      % an infinite loop, which breaks when the convergence 
                      % in the labor allocations with a given accuracy 
                      % of 1e-10 average error is achieved 

% 2.2 Inputs for iteration-on-allocation
%---------------------------------------
kprime = k(2:n_rows+1,1:N);   % Series of next-period period capital
k      = k(1:n_rows,1:N);     % Series of current-period capital

% 2.3 Individual labor: iteration-on-allocation
%----------------------------------------------

    h =  pol_bases*vh;        % A guess on labor of N countries; n_rows-by-N
      
    for n = 1:hiter;          % Number of iterations for iteration-
                              % on-allocation 
    
    % Given labor of a country, find consumption of the country using 
    % condition (31) in MMJ (2011)
    c=psi/(1-psi)*A*(1-alpha)*a.*k.^alpha.*h.^-alpha.*(Le-h); % n_rows-by-N
    
    % Recompute labor of country N using condition (30) in MMJ (2011) 
    h_tilda(:,N)=(((c + kprime-k+phi/2*(kprime./k-1).^2.*k)*ones(N,1) -(A*k(:,1:N-1).^alpha.*h(:,1:N-1).^(1-alpha).*a(:,1:N-1))*ones(N-1,1))./(A*k(:,N).^alpha.*a(:,N))).^(1/(1-alpha));
    % n_rows-by-1

    % Given labor of country N and consumption of N countries, find labor 
    % of countries 1,...,N-1 using condition (32) in MMJ (2011)
    for j = 1:N-1   
        h_tilda(:,j)=Le-(tau(N)/tau(j)*c(:,N).^(psi*(1-1/gam(N))-1)./c(:,j).^(psi*(1-1/gam(j))-1).*(Le-h(:,N)).^((1-psi)*(1-1/gam(N)))).^(1/(1-psi)/(1-1/gam(j)));
    end % n_rows-by-(N-1)

    % Restrict labor of N countries to be in the interval [0.5, 1.5] 
    h_tilda = h_tilda.*(h_tilda>0.5).*(h_tilda<1.5)+0.5*(h_tilda<=0.5)+1.5*(h_tilda>=1.5);

    % Update labor of N countries for the subsequent iteration using 
    % damping 
    h = h*(1-hdamp)+h_tilda*hdamp;
    
    % Compute the difference between labor of N countries at the end and
    % the beginning of the iteration (the convergence criterion is adjusted 
    % to the damping parameter as described in MMJ, 2011)
    dif = mean(mean(abs(1-h_tilda./h)))/hdamp;
    
    % Stop iterating if labor of N countries is computed with an average  
    % absolute error of 1e-10    
    if dif<1e-10; break; end
end