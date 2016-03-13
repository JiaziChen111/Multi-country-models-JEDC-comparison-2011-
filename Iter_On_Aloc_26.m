% Iter_On_Aloc_26.m is a routine that solves for consumption and labor in  
% Models 2 and 6 by performing fixed-point iteration-on-allocation on labor 
% of country 1; see Maliar, Maliar and Judd, (2011), "Solving the Multi-
% Country Real Business Cycle Model Using Ergodic Set Methods", Journal of 
% Economic Dynamics and Control 35, 207-228 (henceforth, MMJ, 2011). 
% -------------------------------------------------------------------------
% Inputs:  "h_old_1" is a guess on labor of country 1; 
%          "k", "a" and "kprime" are current-period capital, current-period
%          productivity levels and next-period capital, respectively;  
%          "N" is the number of countries; 
%          "alpha", "gam", "nu", "b", "phi", "A", "tau" are the parameters
%          of the model;
%          "hdamp" and "hiter" are the parameters in iteration-on-allocation 
%
% Output:  "c" and "h" are individual consumption and labor of N countries,
%          respectively
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [c h]  = Iter_On_Aloc_26(h_old_1,k,a,kprime,N,alpha,gam,nu,b,phi,A,tau,hdamp,hiter)

n_rows=size(k,1);             % n_rows is the number of rows in k and a
c = zeros(n_rows,N);          % allocate memory for consumption allocations

% Individual labor: iteration-on-allocation
%------------------------------------------
h(:,1) = h_old_1(:,1);      % A guess on labor of country 1; n_rows-by-1 
    
for n = 1:hiter;            % Number of iterations for iteration-
                            % on-allocation 

% Given labor of country 1, find labor of the other countries using 
% condition (28) in MMJ (2011)
for j = 2:N   
    h(:,j)=(a(:,j).*k(:,j).^alpha./a(:,1)./k(:,1).^alpha*tau(1)*b(1)/tau(j)/b(j)).^(nu(j)/(1+alpha*nu(j))).*h(:,1).^((nu(j)+alpha*nu(j)*nu(1))/(nu(1)+alpha*nu(j)*nu(1)));
end % n_rows-by-N
    
% Given labor of a country, find consumption of the country using condition    
% (29) in MMJ (2011)
for j=1:N
    c(:,j)=(b(j)*h(:,j).^(1/nu(j))/A/(1-alpha)./a(:,j)./k(:,j).^alpha./h(:,j).^-alpha).^-gam(j);     
end % n_rows-by-N

% Recompute labor of country 1 using condition (30) in MMJ (2011) 
h_tilda_1=(((c + kprime-k+phi/2*(kprime./k-1).^2.*k)*ones(N,1) -(A*k(:,2:N).^alpha.*h(:,2:N).^(1-alpha).*a(:,2:N))*ones(N-1,1))./(A*k(:,1).^alpha.*a(:,1))).^(1/(1-alpha));
% n_rows-by-1

% Restrict labor of country 1 to be in the interval [0.5, 1.5] 
h_tilda_1 = h_tilda_1.*(h_tilda_1>0.5).*(h_tilda_1<1.5)+0.5*(h_tilda_1<=0.5)+1.5*(h_tilda_1>=1.5);

% Update labor of country 1 for the subsequent iteration using 
% damping 
h(:,1) = h(:,1)*(1-hdamp)+h_tilda_1*hdamp;
    
% Compute the difference between labor of country 1 at the end and the
% beginning of the iteration (the convergence criterion is adjusted 
% to the damping parameter as described in MMJ, 2011)
dif = mean(mean(abs(1-h_tilda_1./h(:,1))))/hdamp; 

% Stop iterating if labor of country 1 is computed with an average  
% absolute error of 1e-10    
if dif<1e-10; break; end
end