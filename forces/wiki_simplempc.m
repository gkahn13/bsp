% FORCES - Fast interior point code generation for multistage problems.
% Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
% Automatic Control Laboratory, ETH Zurich.
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% simple MPC - double integrator
% 
%  min   xN'*P*xN + sum_{i=1}^{N-1} xi'*Q*xi + ui'*R*ui
% xi,ui
%       s.t. x1 = x
%            x_i+1 = A*xi + B*ui  for i = 1...N-1
%            xmin <= xi <= xmax   for i = 1...N
%            umin <= ui <= umax   for i = 1...N
%
% and P is solution of Ricatti eqn. from LQR problem

% make sure you have FORCES on your path
addpath('FORCES');

%% system
nx = 2;
nu = 1;
A = [1.1 1; 0 1];
B = [1; 0.5];


%% MPC setup
N = 10;
Q = eye(nx);
R = eye(nu);
[~,P] = dlqr(A,B,Q,R);
umin = -0.5;       umax = 0.5;
xmin = [-5; -5]; xmax = [5; 5];   

%% FORCES multistage form - zi = [xi, ui] for i=1...N-1 and zN = xN

stages = MultistageProblem(N);

for i = 1:N
    % initial stage
    if( i==1 )
        
        % dimension
        stages(i).dims.n = nx+nu; % number of stage variables
        stages(i).dims.r = 2*nx;  % number of equality constraints        
        stages(i).dims.l = nx+nu; % number of lower bounds
        stages(i).dims.u = nx+nu; % number of upper bounds
        stages(i).dims.p = 0;     % number of polytopic constraints
        stages(i).dims.q = 0;     % number of quadratic constraints
        
        % cost
        stages(i).cost.H = blkdiag(Q,R);
        stages(i).cost.f = zeros(stages(i).dims.n,1);
        
        % lower bounds
        stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
        stages(i).ineq.b.lb = [xmin; umin]; % lower bound for this stage variable
        
        % upper bounds
        stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
        stages(i).ineq.b.ub = [xmax; umax]; % upper bound for this stage variable
        
        % equality constraints
        stages(i).eq.C = [eye(nx), zeros(nx,nu); A, B];
        params(1) = newParam('z1',1,'eq.c'); % RHS of first eq. constr. is a parameter: [x0, 0]
        
    end
    
    % stages along horizon
    if( i>1 && i<N )       
        
        % dimension
        stages(i).dims.n = nx+nu; % number of stage variables
        stages(i).dims.r = nx;    % number of equality constraints        
        stages(i).dims.l = nx+nu; % number of lower bounds
        stages(i).dims.u = nx+nu; % number of upper bounds
        stages(i).dims.p = 0;     % number of polytopic constraints
        stages(i).dims.q = 0;     % number of quadratic constraints
        
        % cost
        stages(i).cost.H = blkdiag(Q,R);
        stages(i).cost.f = zeros(stages(i).dims.n,1);
        
        % lower bounds
        stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
        stages(i).ineq.b.lb = [xmin; umin]; % lower bound for this stage variable
        
        % upper bounds
        stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
        stages(i).ineq.b.ub = [xmax; umax]; % upper bound for this stage variable
        
        % equality constraints
        stages(i).eq.C = [A, B];
        stages(i).eq.c = zeros(nx,1);
        if( i==2 )
            stages(i).eq.D = [zeros(nx,nx+nu); -eye(nx), zeros(nx,nu)];
        else
            stages(i).eq.D = [-eye(nx), zeros(nx,nu)];
        end
        
    end
    
    % final stage
    if( i==N )
        
        % dimension
        stages(i).dims.n = nx;    % number of stage variables
        stages(i).dims.r = 0;     % number of equality constraints        
        stages(i).dims.l = nx;    % number of lower bounds
        stages(i).dims.u = nx;    % number of upper bounds
        stages(i).dims.p = 0;     % number of polytopic constraints
        stages(i).dims.q = 0;     % number of quadratic constraints
        
        % cost
        stages(i).cost.H = P;
        stages(i).cost.f = zeros(stages(i).dims.n,1);
        
        % lower bounds
        stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
        stages(i).ineq.b.lb = xmin; % lower bound for this stage variable
        
        % upper bounds
        stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
        stages(i).ineq.b.ub = xmax; % upper bound for this stage variable
        
        % equality constraints        
        stages(i).eq.D = -eye(nx);
        
    end
end

%% define outputs of the solver
outputs = newOutput('u1',1,nx+1:nx+nu);

%% solver settings
codeoptions = getOptions('myMPC');

%% generate code
generateCode(stages,params,codeoptions,outputs);


%% simulate
x1 = [-4; 2];
kmax = 30;
X = zeros(nx,kmax+1); X(:,1) = x1;
U = zeros(nu,kmax);
problem.z1 = zeros(2*nx,1);
for k = 1:kmax
    problem.z1(1:nx) = X(:,k);
    [solverout,exitflag,info] = myMPC(problem);
    if( exitflag == 1 )
        U(:,k) = solverout.u1;
    else
        info
        error('Some problem in solver');
    end
    X(:,k+1) = A*X(:,k) + B*U(:,k);
end

%% plot
figure(1); clf;
subplot(2,1,1); grid on; title('states'); hold on;
plot([1 kmax], [xmax xmax]', 'r--'); plot([1 kmax], [xmin xmin]', 'r--');
ylim(1.1*[min(xmin),max(xmax)]); stairs(1:kmax,X(:,1:kmax)');
subplot(2,1,2);  grid on; title('input'); hold on;
plot([1 kmax], [umax umax]', 'r--'); plot([1 kmax], [umin umin]', 'r--');
ylim(1.1*[min(umin),max(umax)]); stairs(1:kmax,U(:,1:kmax)');