% Productivity.m is a routine that draws random series of the productivity 
% shocks and simulates the corresponding series of the productivity levels; 
% see Maliar, Maliar and Judd, (2011), "Solving the Multi-Country Real 
% Business Cycle Model Using Ergodic Set Methods", Journal of Economic 
% Dynamics and Control 35, 207-228 (henceforth, MMJ, 2011). 
% -------------------------------------------------------------------------
% Inputs:  "T" is the simulation length; T>=1;
%          "N" is the number of countries; N>=1;
%          "a_init" is the initial condition for the productivity levels of
%          N countries; 1-by-N;
%          "rho" and "sigma" are the parameters of the model

% Output:  "a" are the time series of the productivity levels of N countries; 
%          T-by-N
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function a = Productivity(T,N,a_init,sigma,rho)

EPSI = randn(T,1);  % A random draw of common-for-all-countries productivity 
                    % shocks for T periods; T-by-1 
                    
epsi = randn(T,N);  % A random draw of country-specific productivity shocks 
                    % for T periods and N countries; T-by-N

epsi = (epsi+EPSI*ones(1,N))*sigma; 
                    % Compute the error terms in the process for productivity 
                    % level using condition (3) in MMJ (2011); T-by-N

a(1,1:N) = a_init;  % Initial condition for the productivity levels; 1-by-N

for t = 1:T-1; 
    a(t+1,:) = a(t,:).^rho.*exp(epsi(t+1,:)); 
                    % Compute the next-period productivity levels using 
                    % condition (3) in MMJ (2011); 1-by-N 
end;