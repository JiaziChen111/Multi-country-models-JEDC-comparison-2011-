% MATLAB software that solves Models 1 and 5 using the stochastic-simulation 
% algorithm (SSA) and the cluster-grid algorithm (CGA) as described in the 
% article "Solving the Multi-Country Real Business Cycle Model Using Ergodic 
% Set Methods" by Serguei Maliar, Lilia Maliar and Kenneth L. Judd, (2011), 
% Journal of Economic Dynamics and Control 35, 207-228 (henceforth, MMJ, 2011).   
%
% This version: March 31, 2011 (earlier versions: 2004, 2007, 2009)
% 
% ------------------------------------------------------------------------
% The software that solves Models 1 and 5 uses the following files: 
% ------------------------------------------------------------------------
% 1. "Main_SSA_CGA_15.m"  computes a first-degree polynomial solution under 
%                         SSA and second-degree polynomial solution under  
%                         CGA for Models 1 and 5
% 2. "Iter_On_Aloc_15.m"  performs (fixed-point) iteration-on-allocation on
%                         consumption of country 1
% 3. "Simulation_15.m"    simulates a time-series solution given the 
%                         computed policy functions 
% 4. "Accuracy_Test_15.m" computes approximation errors in the optimality 
%                         conditions on a given set of points in the state 
%                         space 
% 5. "Productivity.m"     generates random draws of the productivity shocks  
%                         and simulates the corresponding series of the  
%                         productivity levels
% 6. "Clusters.m"         constructs clusters from simulated series and 
%                         computes clusters' centers to be used as a grid
% 7. "Polynomial_2d.m"    constructs a complete first- and second-degree 
%                         polynomials 
% 8. "Monomials_1.m"      constructs integration nodes and weights for an N-
%                         dimensional monomial (non-product) integration rule 
%                         with 2N nodes 
% 9. "Monomials_2.m"      constructs integration nodes and weights for an N-
%                         dimensional monomial (non-product) integration rule 
%                         with 2N^2+1 nodes
%10. "Monomials_3.m"      constructs integration nodes and weights for an N-
%                         dimensional monomial integration rule with 2^N 
%                         nodes
%11. "Draws_Solut.mat"    contains the series of the productivity levels of 
%                         length 10,000 for 10 countries that are used for 
%                         computing solutions
%12. "Test_Points.mat"    contains the series of the productivity levels of 
%                         length 10,200 for 10 countries that are used for 
%                         implementing the accuracy tests
% -------------------------------------------------------------------------
% Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
% The code may be used, modified and redistributed under the terms provided 
% in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

% 1. Choose the model, number of countries and simulation length
% --------------------------------------------------------------
Model = 5;       % Choose the model (either Model=1 or Model=5)
N     = 2;       % Choose the number of countries, 1<=N<=10 (note that the 
                 % code also works for the one-country case, N=1)
T     = 10000;   % Choose the simulation length for the solution procedure,
                 % T<=10,000   
                 
% To solve models with N>10 or T>10,000, one needs to simulate new series
% of the productivity levels by enabling the code in paragraph 6  

% 2. Common-for-all-countries parameters
% --------------------------------------
alpha   = 0.36;   % Capital share in output
beta    = 0.99;   % Discount factor
phi     = 0.5;    % Adjustment-cost parameter "phi" 
rho     = 0.95;   % Persistence of the log of the productivity level
sigma   = 0.01;   % Standard deviation of shocks to the log of the 
                  % productivity level
vcv = sigma^2*(eye(N)+ones(N,N)); 
                  % Variance-covariance matrix of the countries' productivity 
                  % shocks in which diagonal terms are equal to 2*sigma^2   
                  % and in which off-diagonal terms are equal to sigma^2; 
                  % this vcv follows from condition (3) of MMJ (2011), 
                  % according to which a country's shock has both common-for-
                  % all-countries and country-specific components; N-by-N

% 3. The utility-function parameter, "gamma"
% -----------------------------------------
gam = ones(1,N);                     % "gamma" is identical for all countries 
                                     % (Model 1 is the default)
if (Model == 5)&&(N>1);              % In Model 5 with N>1 countries, ...
    gam = (0.25:((1-0.25)/(N-1)):1); % "gamma" is uniformly distributed across 
                                     % countries in the interval [0.25, 1]
end

% 4. The normalizing constant, A, and welfare weights, tau
% --------------------------------------------------------
A       = (1-beta)/alpha/beta;    % The constant in the production function  
for j=1:N
   tau(j) = 1/A^(-1/gam(j));      % The welfare weight of country j
end

% The above normalization ensures that steady state of capital of all 
% countries is equal to one 

% 5. Initial condition
% --------------------
k(1,1:N) = 1;  % Initial condition for capital (is equal to steady state)
a(1,1:N) = 1;  % Initial condition for the productivity level (is equal to 
               % steady state)

% 6. Construct the productivity levels, a, for the solution procedure 
% -------------------------------------------------------------------
% a = Productivity(T,N,a(1,1:N),sigma,rho); 
                               % Generate a random draw of the productivity 
                               % shocks and simulate the corresponding  
                               % series of the productivity levels of length   
                               % T periods for N countries 
% save Draws_Solut a;          % Save the series of the productivity levels  
                               % into a file "Draws_Solut.mat" 
load Draws_Solut;              % Load the previously saved series of the 
                               % productivity levels of length 10,000 for 
                               % 10 countries
a = a(1:T,1:N);                % Restrict the series of the productivity 
                               % levels to the given T<=10,000 and N<=10
                               
% 7. The matrices for the SSA and CGA polynomial coefficients 
% ------------------------------------------------------------                             
npol_2d = 1+2*N+N*(2*N+1);  % Number of terms in a complete second-degree 
                            % polynomial for the case of N countries
VK = zeros(npol_2d,N,2);    % Matrix of polynomial coefficients of the 
                            % capital policy functions of N countries for two  
                            % solution methods, SSA and CGA; npol_2d-by-N-by-2 
VC = zeros(npol_2d,1,2);    % Matrix of polynomial coefficients of the 
                            % consumption policy function of country 1 for  
                            % two solution methods, SSA and CGA; 
                            % npol_2d-by-1-by-2
% _________________________________________________________________________
                               
% The first-degree polynomial solution under SSA
% _________________________________________________________________________
%
tic;                  % Start counting time for computing the SSA solution
                            
% 8. The SSA parameters  
% ---------------------
kdamp     = 0.05;     % Damping parameter for (fixed-point) iteration on 
                      % the coefficients of the capital policy functions
cdamp     = 0.01;     % Damping parameter for (fixed-point) iteration on
                      % consumption allocations
citer     = 10;       % Number of iterations on consumption allocations in
                      % iteration-on-allocation per iteration on the
                      % coefficients of the capital policy functions
dif_SSA	  = 1e+10;    % Convergence criterion (initially is not satisfied)

% To achieve convergence under N>10, one may need to modify the values of 
% the algorithm's parameters (such as cdamp, kdamp, citer), to refine initial 
% guess and to change the interval in which consumption of country 1 is forced 
% to stay in iteration-on-allocation. 

% 9. Initialize the first-degree capital policy functions of N countries 
% for SSA 
%-----------------------------------------------------------------------                          
vk_1d  = [zeros(1,N); diag(ones(1,N)*0.9);diag(ones(1,N)*0.1)]; 
% Matrix of polynomial coefficients of size (1+2N)-by-N: for each country  
% (each column), 1+2N rows correspond to a constant, N coefficients on the
% capital stocks, k(t,1),...,k(t,N), and N coefficients on the productivity 
% levels, a(t,1),...,a(t,N)

% As an initial guess, assume that a country's j capital depends only on 
% its own capital and productivity level as k(t+1,j)=0.9*k(t,j)+0.1*a(t,j); 
% (note that if we are in the steady state, k(t+1,j)=0.9*k(t,j)+0.1*a(t,j)=1)

% Note that diag(ones(1,N)*q) delivers an N-by-N matrix with  diagonal 
% entries equal to q. 

% 10. Initialize the consumption and capital series
% -------------------------------------------------
c_old_1 = ones(T,1)*A; % A guess on consumption series of country 1; this 
                       % series is used both to initialize iteration-on-
                       % allocation (initially, consumption of country 1 is 
                       % set at the steady-state level equal to A) and to 
                       % check the convergence on the subsequent iteration; 
                       % T-by-1
k_old = ones(T+1,N);   % Initialize the series of next-period capital of N
                       % countries; these series are used to check the
                       % convergence on the subsequent iteration (initially, 
                       % capital can take any value); (T+1)-by-N
                       
% 11. The main iterative cycle of SSA
% -----------------------------------              
while dif_SSA > 1e-7;   % The convergence criterion (which is unit free 
                        % because dif_SSA is unit free)
%
% 11.1 Generate time series for capital
% -------------------------------------
for t = 1:T 
    x(t,:) = [1 k(t,:) a(t,:)];  % The first-degree polynomial bases at t
    k(t+1,:) = x(t,:)*vk_1d;     % Compute next-period capital using vk_1d
end

% 11.2 Compute consumption series of all countries using iteration-on-
% allocation (in Models 1 and 5, consumption of country 1 is computed by 
% fixed-point iteration, and consumption of the other countries is found 
% analytically, given consumption of country 1)
%--------------------------------------------------------------------------
[c]  = Iter_On_Aloc_15(c_old_1,k(1:T,:),a(1:T,:),k(2:T+1,:),N,alpha,gam,phi,A,tau,cdamp,citer); 
% T-by-N

% 11.3 Evaluate the percentage (unit-free) difference between the series  
% from the previous and current iterations
% ----------------------------------------------------------------------
dif_SSA = (mean(mean(abs(1-k./k_old)))/kdamp+mean(mean(abs(1-c(:,1)./c_old_1)))/cdamp)/2
                   % The convergence criterion is adjusted to the damping 
                   % parameters as described in MMJ (2011)

% 11.4 Monte Carlo realizations of the right side of the Euler equation, y, 
% in condition (41) in MMJ (2011)
%--------------------------------------------------------------------------
for j = 1:N
   y(1:T-1,j) = beta*c(2:T,j).^(-1/gam(j))./c(1:T-1,j).^(-1/gam(j))./(1+phi*(k(2:T,j)./k(1:T-1,j)-1)).*(1+alpha*A*k(2:T,j).^(alpha-1).*a(2:T,j)+phi/2*(k(3:T+1,j)./k(2:T,j)-1).*(k(3:T+1,j)./k(2:T,j)+1)).*k(2:T,j);
end  % (T-1)-by-N

% 11.5 Compute and update the coefficients of the capital policy functions 
% ------------------------------------------------------------------------
vk_hat_1d  = x(1:T-1,:)\y;    % Compute new coefficients of the capital 
                              % policy functions using the linear least-squares 
                              % method by QR factorization; this method is
                              % implemented with "backslash operator" (in 
                              % which explanatory and independent variables 
                              % are given by x and y, respectively)
vk_1d = kdamp*vk_hat_1d + (1-kdamp)*vk_1d; 
                              % Update the coefficients of the capital  
                              % policy functions using damping 
                                     
% 11.6 Store the series to be used on the subsequent iteration 
%-------------------------------------------------------------
c_old_1 = c(:,1);  % Store the series for consumption of country 1 to be used 
                   % as an input for performing iteration-on-allocation and   
                   % for checking the convergence on the subsequent iteration
k_old = k;         % Store capital series to be used for checking the 
                   % convergence on the subsequent iteration
end;

% 12. The SSA solution output 
% ---------------------------
vc_1d = x\c(:,1);           % Compute the coefficients of the (first-degree) 
                            % consumption policy function of country 1 
VK(1:1+2*N,1:N,1) = vk_1d;  % Store the coefficients of the (first-degree)  
                            % capital policy functions of N countries 
                            % computed by SSA
VC(1:1+2*N,1,1) = vc_1d;    % Store the coefficients of the consumption 
                            % policy function of country 1 computed by SSA
time_SSA        = toc;      % Time needed to compute the SSA solution

% _________________________________________________________________________

% The second-degree polynomial solution under CGA
% _________________________________________________________________________

tic;                   % Start counting time for computing the CGA solution

% 13. Create a grid by clustering the simulated data and construct a complete 
% second-degree polynomial on the grid
% -------------------------------------------------------------------------
n_G = 500;                    % Number of clusters (grid points) to be created
Data = [k(1:T,:) a(1:T,:)];   % Matrix of the data to be clustered: it 
                              % includes N capital stocks and N productivity 
                              % levels for each t=1,...,T; T-by-2N
[Grid] = Clusters(Data,n_G);  % Construct clusters and compute their 
                              % centers; the centers of the clusters, Grid,   
                              % are our grid points (each cluster's center  
                              % is described by N capital stocks and N 
                              % productivity levels); n_G-by-2N                              
X0_G = Polynomial_2d(Grid);   % Form a complete second-degree polynomial on
                              % the grid points (centers of the clusters)  
                              % for future regressions: the matrix X0_G
                              % consists of a column of ones, the linear 
                              % and quadratic polynomial bases on the grid
                              % points

% 14. Select the integration method 
% ---------------------------------
[n_nodes,epsi_nodes,weight_nodes] = Monomials_1(N,vcv);
                             % Monomial integration rule with 2N nodes
                             
%[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(N,vcv);
                             % Monomial integration rule with 2N^2+1 nodes

%[n_nodes,epsi_nodes,weight_nodes] = Monomials_3(N,vcv);
                             % Monomial integration rule with 2^N nodes;  
                             % this rule coincides with Gauss-Hermite 
                             % quadrature (product) integration rule with 
                             % 2 nodes in each dimension 

%n_nodes = 1; epsi_nodes = zeros(1,N); weight_nodes = 1;                             
                             % Gauss-Hermite quadrature integration rule 
                             % with one node 
                              
% 15. Initialize the second-degree capital policy functions for CGA
% -----------------------------------------------------------------
vk_2d = VK(:,1:N,1);         % Start from the first-degree polynomial solution
                             % obtained under SSA

% 16. Initialize the values for next-period capital and the country's 1 
% current- and next-period consumption on the grid
%----------------------------------------------------------------------
for i = 1:n_G;
    
   k1_old_G(1,1:N,i) = [1 Grid(i,:)]*vk_1d; 
                      % Compute an initial guess on next-period capital 
                      % (chosen at t) on the grid using the SSA capital  
                      % policy functions; 1-by-N-by-n_G
   a1_G(1:n_nodes,1:N,i) = (ones(n_nodes,1)*Grid(i,N+1:2*N)).^rho.*exp(epsi_nodes); 
                      % Compute the next-period productivity levels on the  
                      % grid for all integration node using condition (3) 
                      % in MMJ (2011); n_nodes-by-N-by-n_G
   c0_old_1_G(1,1,i) = [1 Grid(i,:)]*vc_1d;
                      % Compute an initial guess on current-period consump-
                      % tion of country 1 on the grid using the SSA solution; 
                      % 1-by-1-by-n_G
   c1_old_1_G(1:n_nodes,1,i) = [ones(n_nodes,1) ones(n_nodes,1)*k1_old_G(:,:,i) a1_G(1:n_nodes,1:N,i)]*vc_1d;
                      % Compute an initial guess on next-period consumption 
                      % of country 1 on the grid for all integration nodes 
                      % using the SSA solution; n_nodes-by-1-by-n_G 
  
   % k1_old_G is used for checking the convergence on the subsequent iteration; 
   % c0_old_1_G and c1_old_1_G are used both for performing iteration-on-
   % allocation and for checking the convergence on the subsequent iteration
   
   % All the variables in paragraph 16 are defined to have a similar three-
   % dimensional structure for consistency
end

% 17. The CGA parameters
% ----------------------
kdamp     = 0.1;      % Damping parameter for (fixed-point) iteration on 
                      % the coefficients of the capital policy functions
cdamp     = 0.05;     % Damping parameter for (fixed-point) iteration on
                      % consumption allocations
citer     = 3;        % Number of iterations on consumption allocations in
                      % iteration-on-allocation per iteration on the
                      % coefficients of the capital policy functions
dif_CGA	  = 1e+10;    % Convergence criterion (initially is not satisfied)

% To achieve convergence under N>10, one may need to modify the values of 
% the algorithm's parameters (such as cdamp, kdamp, citer), to refine initial 
% guess and to change the interval in which consumption of country 1 is  
% forced to stay in iteration-on-allocation 

% 18. The main iterative cycle of CGA
% -----------------------------------              
while dif_CGA > 1e-7; % The convergence criterion (which is unit free 
                      % because dif_CGA is unit free)

    for i=1:n_G;      % For each grid point (cluster's center) ... 
        
        % 18.1 Variables in grid point i 
        % ------------------------------        
        k0 = Grid(i,1:N);           % N capital stocks of period t 
        a0 = Grid(i,N+1:2*N);       % N productivity levels of period t 
        X0 = X0_G(i,:);             % Complete (second-degree) polynomial 
                                    % bases (at t) formed on k0 and a0; 
                                    % 1-by-npol_2d 
        c0(1,1) = c0_old_1_G(1,1,i);% A guess on current-period consumption
                                    % of country 1 
        c1(1:n_nodes,1) = c1_old_1_G(1:n_nodes,1,i);                          
                                    % A guess on next-period consumption of 
                                    % country 1 for n_nodes states 
        a1 = a1_G(1:n_nodes,1:N,i); % Next-period productivity levels for 
                                    % n_nodes states 
        
        % 18.2 Capital and consumption choices at t
        % -----------------------------------------
        k1 = X0*vk_2d; 
        % Compute a row-vector of capital of period t+1 (chosen at t) using
        % the second-degree capital policy functions; 1-by-N 

        c0(1,1:N)  = Iter_On_Aloc_15(c0(1,1),k0,a0,k1,N,alpha,gam,phi,A,tau,cdamp,citer); 
        % Compute current-period consumption of all countries (consumption 
        % of country 1 is computed by fixed-point iteration, and consumption 
        % of the other countries is found analytically, given consumption 
        % of country 1); 1-by-N

        % 18.3 Capital and consumption choices at t+1
        %--------------------------------------------
        k1_dupl = ones(n_nodes,1)*k1; 
        % Duplicate k1 n_nodes times to create a matrix with n_nodes identical
        % rows; n_nodes-by-N 

        X1 = Polynomial_2d([k1_dupl a1]);
        % Form a complete second-degree polynomial (at t+1) on k1_dupl and 
        % a1; n_nodes-by-npol_2d 
        
        k2 = X1*vk_2d; 
        % Compute capital of period t+2 (chosen at t+1) using the second- 
        % degree capital policy functions; n_nodes-by-N 

        c1(1:n_nodes,1:N)  = Iter_On_Aloc_15(c1(1:n_nodes,1),k1_dupl,a1,k2,N,alpha,gam,phi,A,tau,cdamp,citer); 
        % Compute next-period consumption of all countries (consumption 
        % of country 1 is computed by fixed-point iteration, and consumption 
        % of the other countries is found analytically, given consumption 
        % of country 1); n_nodes-by-N

        % 18.4. The integral in the right side of the Euler equation (that 
        % is equal to next-period capital) in condition (41) of MMJ (2011)
        %------------------------------------------------------------------        
        for j = 1:N
            k1_hat_G(i,j) = weight_nodes'*(beta*c1(1:n_nodes,j).^(-1/gam(j))/c0(1,j)^(-1/gam(j))/(1+phi*(k1(1,j)/k0(1,j)-1)).*(1+alpha*A*k1(1,j)^(alpha-1)*a1(1:n_nodes,j)+phi/2*(k2(1:n_nodes,j)/k1(1,j)-1).*(k2(1:n_nodes,j)/k1(1,j)+1))*k1(1,j));
            % Approximate the integral as a weighted sum of the integrand
            % in n_nodes with the integration weights weight_nodes
        end 

        % 18.5 Variables in all grid points 
        %----------------------------------
        c0_new_1_G(1,1,i) = c0(1,1);
        % Store current-period consumption of country 1 to be used as an 
        % input for performing iteration-on-allocation and for checking 
        % the convergence on the subsequent iteration; 1-by-1-by-n_G 

        c1_new_1_G(1:n_nodes,1,i) = c1(1:n_nodes,1);   
        % Store next-period consumption of country 1 to be used as an 
        % input for performing iteration-on-allocation and for checking the  
        % convergence on the subsequent iteration; n_nodes-by-1-by-n_G 

        k1_new_G(1,1:N,i) = k1; 
        % Store next-period capital to be used for checking the convergence  
        % on the subsequent iteration; 1-by-N-by-n_G 
    end  
    
% 18.6 Evaluate the percentage (unit-free) difference between the values on  
% the grid from the previous and current iterations
% -------------------------------------------------------------------------
dif_CGA = (mean(mean(abs(1-k1_new_G./k1_old_G)))/kdamp+mean(mean(abs(1-c0_new_1_G./c0_old_1_G)))/cdamp)/2
                   % The convergence criterion is adjusted to the damping 
                   % parameters as described in MMJ (2011)    

% 18.7 Compute and update the coefficients of the capital policy functions 
% ------------------------------------------------------------------------
vk_hat_2d = X0_G\k1_hat_G;                 % Compute the new coefficients of
                                           % the capital policy functions                          
vk_2d = kdamp*vk_hat_2d + (1-kdamp)*vk_2d; % Update the coefficients using
                                           % damping
                                          
% 18.8 Store the values on the grid to be used on the subsequent iteration 
%-------------------------------------------------------------------------
c0_old_1_G = c0_new_1_G; 
                   % Store current-period consumption of country 1 on the 
                   % grid to be used as an input for performing iteration- 
                   % on-allocation and for checking the convergence on 
                   % the subsequent iteration
c1_old_1_G = c1_new_1_G; 
                   % Store next-period consumption of country 1 on the grid
                   % to be used as an input for iteration-on-allocation 
k1_old_G = k1_new_G;     
                   % Store next-period capital on the grid for checking the
                   % convergence on the subsequent iteration          
end

% 19. The CGA solution output
% ---------------------------
vc_2d = X0_G\squeeze(c0_new_1_G); 
                     % Compute the coefficients of the (second-degree) 
                     % consumption policy function of country 1 
VK(:,1:N,2) = vk_2d; % Store the coefficients of the (second-degree) capital   
                     % policy functions computed by CGA
VC(:,1,2) = vc_2d;   % Store the coefficients of the consumption policy 
                     % function of country 1 computed by CGA
time_CGA  = toc;   % Time needed to compute the CGA solution

% 20. Accuracy test of the SSA and CGA solutions: errors on spheres of 
% radia r = 0.01, 0.1, 0.3 and on a stochastic simulation 
% --------------------------------------------------------------------

% 20.1 Specify a set of points on which the accuracy is evaluated
%----------------------------------------------------------------
load Test_Points;  
         % This data file contains a set of state variables for testing 
         % accuracy for N=10 countries from Juillard and Villemot (2011);
         % k_r001 and a_r001 are 1,000 points on a sphere of radius r=0.01;
         % k_r01 and a_r01 are 1,000 points on a sphere of radius r=0.1;
         % k_r03 and a_r03 are 1,000 points on a sphere of radius r=0.3;
         % a_sim is a series of 10,200 observations for the productivity 
         % levels constructed using a random draw 

k_r001 = k_r001(:,1:N);        % Restrict the series of capital for the test  
                               % on a sphere of radius r=0.01 to the given 
                               % N<=10 
k_r01  = k_r01(:,1:N);         % The same for r=0.1
k_r03  = k_r03(:,1:N);         % The same for r=0.3

a_r001 = a_r001(:,1:N);        % Restrict the series of the productivity 
                               % levels for the test on a sphere of radius   
                               % r=0.01 to the given N<=10 
a_r01  = a_r01(:,1:N);         % The same for r=0.1
a_r03  = a_r03(:,1:N);         % The same for r=0.3
          
T_sim = 10200;                 % Choose the simulation length for the test 
                               % on a stochastic simulation, T_sim<=10,200 

a_sim = a_sim(1:T_sim,1:N);    % Restrict the series of the productivity 
                               % levels for the test on a stochastic 
                               % simulation to the given T_sim<=10,200 and 
                               % N<=10                             
          
k0_sim(1,1:N) = 1;  % Initial condition for capital (equal to steady state)
                 
% To implement the tests on spheres for models with N>10, one needs to 
% construct new sets of points as described in Juillard and Villemot (2011); 
% to implement the test on a stochastic simulation with N>10 or T_sim>10,200, 
% one needs to simulate new series of the productivity levels with larger N 
% or T_sim by enabling the code in paragraph 6

% 20.2 Compute errors on a stochastic simulation for SSA (m=1) and CGA (m=2)
% ------------------------------------------------------------------------
for m = 1:2 
    
    % Simulate the time series solution under the given capital-policy- 
    % function coefficients, VK(:,:,m) with m=1 for SSA and m=2 for CGA 
    % (for m=1, the coefficients on the quadratic polynomial bases are 
    % all zeros)
    %--------------------------------------------------------------------
    [c_sim k_sim] = Simulation_15(VK(:,:,m),VC(:,:,m),a_sim,k0_sim,alpha,gam,phi,A,tau);
    
    % Errors across 1,000 points on radius r=0.01
    % -------------------------------------------
    discard = 0;    % Do not discard observations
    [Errors_mean(m,1),Errors_max(m,1),Errors_max_EE(m,1),Errors_max_MUC(m,1),Errors_max_MUL(m,1),Errors_max_RC(m,1),time_test(m,1)] = Accuracy_Test_15(k_r001,a_r001,VK(:,:,m),VC(:,1,m),alpha,gam,phi,beta,A,tau,rho,vcv,discard);
    
    % Errors across 1,000 points on radius r=0.1
    % ------------------------------------------
    discard = 0; 
    [Errors_mean(m,2),Errors_max(m,2),Errors_max_EE(m,2),Errors_max_MUC(m,2),Errors_max_MUL(m,2),Errors_max_RC(m,2),time_test(m,2)] = Accuracy_Test_15(k_r01,a_r01,VK(:,:,m),VC(:,1,m),alpha,gam,phi,beta,A,tau,rho,vcv,discard);

    % Errors across 1,000 points on radius r=0.3
    % ------------------------------------------
    discard = 0; 
    [Errors_mean(m,3),Errors_max(m,3),Errors_max_EE(m,3),Errors_max_MUC(m,3),Errors_max_MUL(m,3),Errors_max_RC(m,3),time_test(m,3)] = Accuracy_Test_15(k_r03,a_r03,VK(:,:,m),VC(:,1,m),alpha,gam,phi,beta,A,tau,rho,vcv,discard);
    
    % Errors across 10,000 points on a stochastic simulation
    % ------------------------------------------------------
    discard = 200; % Discard the first 200 observations to remove the effect
                   % of the initial conditions 
    [Errors_mean(m,4),Errors_max(m,4),Errors_max_EE(m,4),Errors_max_MUC(m,4),Errors_max_MUL(m,4),Errors_max_RC(m,4),time_test(m,4)] = Accuracy_Test_15(k_sim,a_sim,VK(:,:,m),VC(:,1,m),alpha,gam,phi,beta,A,tau,rho,vcv,discard);

    % Errors_mean    is the unit-free average absolute approximation error  
    %                across 4N+1 optimality conditions (in log10) 
    % Errors_max     is the unit-free maximum absolute approximation error   
    %                across 4N+1 optimality conditions (in log10) 
    % Errors_max_EE  is the unit-free maximum absolute approximation error  
    %                across N Euler equations (in log10)
    % Errors_max_MUC is the unit-free maximum absolute approximation error   
    %                across N conditions on marginal utility of consumption 
    %                (optimality conditions w.r.t. consumption) (in log10)
    % Errors_max_MUL is the unit-free maximum absolute approximations error   
    %                across N conditions on marginal utility of labor 
    %                (optimality conditions w.r.t. labor) (in log10)
    % Errors_max_RC  is the unit-free maximum absolute approximation error   
    %                in the aggregate resource constraint (in log10)
    % The errors in capital-accumulation equations are zero by construction  
end

% 21. Display the results for SSA 
% -------------------------------
disp(' '); disp('SSA OUTPUT'); disp(' '); 
disp('RUNNING TIME (in seconds):'); disp('');
disp('a) for computing the solution'); disp(time_SSA);
disp('b) for implementing the accuracy test'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim'); disp(time_test(1,:));
disp('APPROXIMATION ERRORS (log10):'); disp(''); 
disp('a) mean error across 4N+1 optimality conditions'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_mean(1,:))
disp('b) max error across 4N+1 optimality conditions'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max(1,:))
disp('c) max error across N Euler equations');   
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_EE(1,:))
disp('d) max error across N conditions on MUC'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_MUC(1,:))
disp('e) max error across N conditions on MUL'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_MUL(1,:))
disp('f) max error in the resource constraint'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_RC(1,:))

% 22. Display the results for CGA 
% -------------------------------
disp(' '); disp('CGA OUTPUT'); disp(' '); 
disp('RUNNING TIME (in seconds):'); disp('');
disp('a) for computing the solution'); disp(time_CGA);
disp('b) for running the accuracy test'); 
disp('   r=0.01    r=0.1     r=0.3    st.sim'); disp(time_test(2,:));
disp('APPROXIMATION ERRORS (log10):'); disp(''); 
disp('a) mean error across 4N+1 optimality conditions'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_mean(2,:))
disp('b) max error across 4N+1 optimality conditions'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max(2,:))
disp('c) max error across N Euler equations');   
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_EE(2,:))
disp('d) max error across N conditions on MUC'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_MUC(2,:))
disp('e) max error across N conditions on MUL'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_MUL(2,:))
disp('f) max error in the resource constraint'); 
disp('   r=0.01    r=0.1     r=0.3     st.sim');disp(Errors_max_RC(2,:))

save M5_N2