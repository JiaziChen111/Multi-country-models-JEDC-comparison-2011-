% Accuracy_Test_15.m is a routine for evaluating accuracy of the solutions    
% to Models 1 and 5: it computes approximation errors in the optimality  
% conditions on a given set of points in the state space; see Maliar, Maliar 
% and Judd, (2011), "Solving the Multi-Country Real Business Cycle Model Using 
% Ergodic Set Methods", Journal of Economic Dynamics and Control 35, 207-228 
% (henceforth, MMJ, 2011), and Juillard and Villemot, (2011), "Multi-Country 
% Real Business Cycle Models: Accuracy Tests and Test Bench", Journal of 
% Economic Dynamics and Control 35, 178-185.
% -------------------------------------------------------------------------
% Inputs:    "k" and "a" are, respectively, current-period capital and 
%            productivity levels, in the given set of points on which the 
%            accuracy is tested; 
%            "vk" and "vc" are, respectively, the coefficients of the 
%            capital policy functions of N countries and the consumption 
%            policy function of country 1;
%            "alpha", "gam", "phi", "beta", "A", "tau", "rho" and "vcv"
%            are the parameters of the model;
%            "discard" is the number of data points to discard 
%
% Outputs:   "Errors_mean" and "Errors_max" are, respectively, the mean and
%            maximum approximation errors across all optimality conditions;
%            "Errors_max_EE", "Errors_max_MUC", "Errors_max_MUL", and 
%            "Errors_max_RC" are the maximum approximation errors disaggre- 
%            gated by optimality conditions 
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [Errors_mean Errors_max Errors_max_EE Errors_max_MUC Errors_max_MUL Errors_max_RC time_test Errors] = Accuracy_Test_15(k,a,vk,vc,alpha,gam,phi,beta,A,tau,rho,vcv,discard)

tic                   % Start counting time for running the test        

[P,N] = size(a);      % Infer the number of points on which accuracy is 
                      % evaluated P and the number of countries N
                      
% 1. Parameters in iteration-on-allocation
%-----------------------------------------
cdamp    = 0.01;      % Damping parameter in iteration-on-allocation
citer    = 1e+6;      % Number of iterations on consumption allocations in
                      % iteration-on-allocation; under a very large citer 
                      % (like here), iteration-on-allocation is in effect 
                      % an infinite loop, which breaks when the convergence 
                      % in the consumption allocations with a given accuracy 
                      % of 1e-10 average error is achieved 

% 2. Integration method for evaluating accuracy 
% ---------------------------------------------
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(N,vcv);
                             % Monomial integration rule with 2N^2+1 nodes

% 3. A complete second-degree polynomial for the test
%----------------------------------------------------
X = Polynomial_2d([k a]);    % Form a complete second-degree polynomial on 
                             % the given set of points 

% 4. Given the solution for capital, compute consumption on the given set 
% of points  
%------------------------------------------------------------------------
for p = 1:P;                 % For each given point, ...     
    
       p                     % Display the point (with the purpose of  
                             % monitoring the progress) 
                                     
        % 4.1 Variables in point p 
        % ------------------------        
        k0 = k(p,1:N);       % N capital stocks of period t
        a0 = a(p,1:N);       % N productivity levels of period t
        X0 = X(p,:);         % Complete (second-degree) polynomial 
                             % bases at t
                                    
        % 4.2 Capital and consumption choices at t
        % ----------------------------------------
        k1 = X0*vk; 
        % Compute a row-vector of capital of period t+1 (chosen at t) using
        % the corresponding capital policy functions; 1-by-N
        
        c0(1,1) = X0*vc; 
        % Compute a guess on current-period consumption of country 1 using 
        % the corresponding consumption policy function of country 1; 1-by-1

        c0(1,1:N)  = Iter_On_Aloc_15(c0(1,1),k0,a0,k1,N,alpha,gam,phi,A,tau,cdamp,citer); 
        % Compute current-period consumption using iteration-on-allocation 
        % (consumption of country 1 is computed by fixed-point iteration,
        % and consumption of the other countries is found analytically, 
        % given consumption of country 1); 1-by-N

        % 4.3 Capital and consumption choices at t+1
        %--------------------------------------------
        a1 = (ones(n_nodes,1)*a0).^rho.*exp(epsi_nodes);    
        % Compute the next-period productivity levels in each integration node
        % using condition (3) in MMJ (2011); n_nodes-by-N

        k1_dupl = ones(n_nodes,1)*k1; 
        % Duplicate k1 n_nodes times to create a matrix with n_nodes identical
        % rows; n_nodes-by-N 

        X1 = Polynomial_2d([k1_dupl a1]);
        % Form a complete second-degree polynomial (at t+1) in the given point 
        
        k2 = X1*vk; 
        % Compute capital of period t+2 (chosen at t+1) using the second-
        % degree capital policy functions; n_nodes-by-N 

        c1(1:n_nodes,1) = X1*vc;
        % Compute a guess on next-period consumption of country 1 using 
        % the corresponding consumption policy function of country 1; 
        % n_nodes-by-1

        c1(1:n_nodes,1:N) = Iter_On_Aloc_15(c1(1:n_nodes,1),k1_dupl,a1,k2,N,alpha,gam,phi,A,tau,cdamp,citer); 
        % Compute next-period consumption using iteration-on-allocation 
        % (consumption of country 1 is computed by fixed-point iteration,
        % and consumption of the other countries is found analytically, 
        % given consumption of country 1); n_nodes-by-N

        % 5. Approximation errors in point p
        %-----------------------------------
        
        % 5.1 Lagrange multiplier associated with the aggregate resource
        % constraint
        %---------------------------------------------------------------
        for j = 1:N
            MUC0j(1,j) = tau(j).*c0(1,j).^(-1./gam(j)); 
            % Compute a country's marginal utility of consumption multiplied 
            % by its welfare weight
        end
        lambda0 = mean(MUC0j,2);
        % An optimality condition w.r.t. consumption of period t equates 
        % the Lagrange multiplier of the aggregate resource constraint of 
        % period t and each country's marginal utility of consumption 
        % multiplied by its welfare weight; to infer the Lagrange multiplier,  
        % we average across N countries; 1-by-1
        
        for j = 1:N
            MUC1j(1:n_nodes,j) = tau(j).*c1(1:n_nodes,j).^(-1./gam(j));
            % Compute a country's marginal utility of consumption multiplied 
            % by its welfare weight
        end
        lambda1 = mean(MUC1j,2);
        % Similarly, the Lagrange multiplier of the aggregate resource 
        % constraint of period t+1 is equal to a country's marginal utility 
        % of consumption multiplied by its welfare weight; to infer the 
        % Lagrange multiplier, we average across N countries; 1-by-n_nodes
        
        % 5.2 Unit-free Euler-equation errors
        %------------------------------------
        for j = 1:N
            Errors(p,j) = 1-weight_nodes'*(beta*lambda1/lambda0/(1+phi*(k1(1,j)/k0(1,j)-1)).*(1+alpha*A*k1(1,j)^(alpha-1)*a1(1:n_nodes,j)+phi/2*(k2(1:n_nodes,j)/k1(1,j)-1).*(k2(1:n_nodes,j)/k1(1,j)+1)));
        % A unit-free Euler-equation approximation error of country j
        end
        
        % 5.2 Unit-free errors in the optimality conditions w.r.t. consumption
        %---------------------------------------------------------------------
        for j = 1:N;
            Errors(p,N+j) = 1-lambda0./(tau(j)*c0(1,j)^(-1/gam(j))); 
        % A unit-free approximation error in the optimality condition w.r.t. 
        % consumption of country j (this condition equates marginal utility 
        % of consumption, multiplied by the welfare weight, and the 
        % Lagrange multiplier of the aggregate resource constraint)
        end
        
        % 5.3 Unit-free errors in the optimality conditions w.r.t. labor 
        %---------------------------------------------------------------
        Errors(p,2*N+1:3*N) = zeros(N,1);
        % These errors  are zero by construction in Models 1 and 5 
        
        % 5.4 Unit-free approximation error in the aggregate resource constraint
        %-----------------------------------------------------------------------
        Errors(p,3*N+1) = 1-(c0(1,1:N) + k1(1,1:N)-k0(1,1:N))*ones(N,1)/((A*k0(1,1:N).^alpha.*a0(1,1:N)-phi/2*(k1(1,1:N)./k0(1,1:N)-1).^2.*k0(1,1:N))*ones(N,1));
        % This error is a unit-free expression of resource constraint (2) 
        % in MMJ (2011)
        
        % 5.5 Approximation errors in the capital-accumulation equation
        %--------------------------------------------------------------
        Errors(p,3*N+2:4*N+1) = zeros(N,1);
        % These errors are always zero by construction

end

% 6. Mean and maximum approximation errors computed after discarding the first 
% "discard" observations
%----------------------------------------------------------------------------

 % 6.1 Approximation errors across all the optimality conditions
 %--------------------------------------------------------------
 Errors_mean = log10(mean(mean(abs(Errors(1+discard:end,:))))); 
 % Average absolute approximation error 
 
 Errors_max = log10(max(max(abs(Errors(1+discard:end,:)))));    
 % Maximum absolute approximation error
 
 % 6.2 Maximum approximation errors disaggregated by the optimality conditions
 %----------------------------------------------------------------------------
 Errors_max_EE = log10(max(max(abs(Errors(1+discard:end,1:N)))));    
 % Across N Euler equations

 Errors_max_MUC = log10(max(max(abs(Errors(1+discard:end,N+1:2*N)))));    
 % Across N optimality conditions w.r.t. consumption (conditions on marginal  
 % utility of consumption, MUC)

 Errors_max_MUL = log10(max(max(abs(Errors(1+discard:end,2*N+1:3*N)))));    
 % Across N optimality conditions w.r.t. labor (conditions on marginal 
 % utility of labor, MUL)

 Errors_max_RC = log10(max(max(abs(Errors(1+discard:end,3*N+1)))));    
 % In the aggregate resource constraint 
 
time_test = toc;     % Time needed to run the test     