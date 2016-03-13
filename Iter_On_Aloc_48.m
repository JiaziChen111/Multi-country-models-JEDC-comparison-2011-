% Iter_On_Aloc_48.m is a routine that solves for consumption and labor in  
% Models 4 and 8 by performing fixed-point iteration-on-allocation on labor 
% of N countries; see Maliar, Maliar and Judd, (2011), "Solving the Multi-
% Country Real Business Cycle Model Using Ergodic Set Methods", Journal of   
% Economic Dynamics and Control 35, 207-228 (henceforth, MMJ, 2011). 
% -------------------------------------------------------------------------
% Inputs:  "h_old" is a guess on labor of N countries; 
%          "k", "a" and "kprime" are current-period capital, current-period
%          productivity levels and next-period capital, respectively;  
%          "N" is the number of countries; 
%          "alpha", "gam", "xi", "mu", "b", "phi", "A", "Le", "tau" are the 
%          parameters of the model;
%          "hdamp" and "hiter" are the parameters in iteration-on-allocation 
%
% Output:  "c" and "h" are individual consumption and labor of N countries,
%          respectively
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [c h]  = Iter_On_Aloc_48(h_old,k,a,kprime,N,alpha,gam,xi,mu,b,phi,A,Le,tau,hdamp,hiter)

    n_rows=size(k,1);             % n_rows is the number of rows in k and a
    c = zeros(n_rows,N);          % allocate memory for consumption allocations
    h_tilda = zeros(n_rows,N);    % allocate memory for labor allocations
    
% Individual labor: iteration-on-allocation
%------------------------------------------

    h = h_old;                    % A guess on labor of N countries; 
                                  % n_rows-by-N 
                                  
    for n = 1:hiter;              % Number of iterations for iteration-
                                  % on-allocation 
    
    % Given labor of a country, find consumption of the country using 
    % condition (33) in MMJ (2011)
    for j=1:N
        c(:,j)=((1-alpha)*A*a(:,j).*h(:,j).^(mu(j)-1).*(alpha*k(:,j).^mu(j)+(1-alpha)*h(:,j).^mu(j)).^(1/mu(j)-1)/b(j)).^xi(j).*(Le-h(:,j));     
    end   % n_rows-by-N
    
    % Recompute labor of country 1 using condition (35) in MMJ (2011) 
    h_tilda(:,1)=(((((c + kprime-k+phi/2*(kprime./k-1).^2.*k)*ones(N,1) -(A*a(:,2:N).*(alpha*k(:,2:N).^(ones(n_rows,1)*mu(2:N))+(1-alpha)*h(:,2:N).^(ones(n_rows,1)*mu(2:N))).^(1./(ones(n_rows,1)*mu(2:N))))*ones(N-1,1))./a(:,1)/A).^mu(1)-alpha*k(:,1).^mu(1))/(1-alpha)).^(1/mu(1));      
    % n_rows-by-1

    % Given labor of country N and consumption of N countries, find labor 
    % of countries 2,...,N using condition (34) in MMJ (2011)
    for j = 2:N   
        h_tilda(:,j) = Le-(((tau(1)/tau(j)*(c(:,1).^(1-1/xi(1))+b(1)*(Le-h(:,1)).^(1-1/xi(1))).^((1/xi(1)-1/gam(1))/(1-1/xi(1))).*c(:,1).^(-1/xi(1))./c(:,j).^(-1/xi(j))).^((1-1/xi(j))/(1/xi(j)-1/gam(j)))-c(:,j).^(1-1/xi(j)))/b(j)).^(1/(1-1/xi(j)));      
    end % n_rows-by-(N-1)

    % Restrict labor of N countries to be in the interval [0.5, 1.5] 
    h_tilda = h_tilda.*(h_tilda>0.5).*(h_tilda<1.5)+0.5*(h_tilda<=0.5)+1.5*(h_tilda>=1.5);

    % Update labor of N countries for the subsequent iteration using 
    % damping 
    h = h*(1-hdamp)+h_tilda*hdamp;
    
    % Compute the difference between labor of N countries at the end the
    % the beginning of the iteration (the convergence criterion is adjusted 
    % to the damping parameter as described in MMJ, 2011)
    dif = mean(mean(abs(1-h_tilda./h)))/hdamp;
    
    % Stop iterating if labor of N countries is computed with an average  
    % absolute error of 1e-10    
    if dif<1e-10; break; end
    end
end
