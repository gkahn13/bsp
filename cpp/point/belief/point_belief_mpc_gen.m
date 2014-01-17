function beliefmpcgen()

% FORCES - Fast interior point code generation for multistage problems.
% Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
% Automatic Control Laboratory, ETH Zurich.

rootDir = pwd(1:strfind(pwd,'bsp')-1);
addpath strcat(rootDir,'/forces');

close all;
clear all;

% problem setup
N = 14;
nx = 2;
ns = 3;
nb = 5;
nu = 2;
stages = MultistageProblem(N+1);

Q = 10*eye(ns);
Qfinal = 10*eye(ns);
R = 1*eye(nu);

% first stage
i=1;
istr = sprintf('%02d',i);

% dimensions
stages(i).dims.n = nb+nu;           % number of stage variables
stages(i).dims.l = nb+nu;           % number of lower bounds
stages(i).dims.u = nb+nu;           % number of upper bounds
stages(i).dims.r = 2*nb;            % number of equality constraints
stages(i).dims.p = 0;               % number of affine constraints
stages(i).dims.q = 0;               % number of quadratic constraints

% cost
stages(i).cost.H = 2*blkdiag(zeros(nx,nx), Q, R);
stages(i).cost.f = zeros(stages(i).dims.n,1);

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
params(1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable

params(end+1) = newParam(['C',istr], i, 'eq.C');
params(end+1) = newParam(['e',istr], i, 'eq.c');

for i = 2:N
    istr = sprintf('%02d',i);
    
    % dimension
    stages(i).dims.n = nb+nu;    % number of stage variables
    stages(i).dims.l = nb+nu;    % number of lower bounds
    stages(i).dims.u = nb+nu;    % number of upper bounds
    stages(i).dims.r = nb;       % number of equality constraints
    stages(i).dims.p = 0;        % number of polytopic constraints
    stages(i).dims.q = 0;        % number of quadratic constraints
    
    % cost
    stages(i).cost.H = 2*blkdiag(zeros(nx,nx), Q, R);
    stages(i).cost.f = zeros(stages(i).dims.n,1);
    
    % lower bounds
    stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
    params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable
    
    % upper bounds
    stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
    params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable
        
    % equality constraints
    params(end+1) = newParam(['C',istr], i, 'eq.C');
    params(end+1) = newParam(['e',istr], i, 'eq.c');
    if( i==2 )
        stages(i).eq.D = [zeros(nb,nb+nu); -eye(nb), zeros(nb,nu)];
    else
        stages(i).eq.D = [-eye(nb), zeros(nb,nu)];
    end
    
end

% final stage
i = N+1;
istr = sprintf('%02d',i);

% dimension
stages(i).dims.n = nb;    % number of stage variables
stages(i).dims.r = 0;    % number of equality constraints
stages(i).dims.l = nb;    % number of lower bounds
stages(i).dims.u = nb;    % number of upper bounds
stages(i).dims.p = 0;     % number of polytopic constraints
stages(i).dims.q = 0;     % number of quadratic constraints

% cost
stages(i).cost.H = 2*blkdiag(zeros(nx,nx), Qfinal);
stages(i).cost.f = zeros(stages(i).dims.n,1);

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb');

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub');

% equality constraints
%stages(i).eq.C = blkdiag(eye(nx), zeros(ns,ns));
%params(end+1) = newParam(['e',istr], i, 'eq.c');
stages(i).eq.D = -eye(nb);

%--------------------------------------------------------------------------
% define outputs of the solver
for i=1:N
    var = sprintf('z%d',i);
    outputs(i) = newOutput(var,i,1:nb+nu);
end
i=N+1;
var = sprintf('z%d',i);
outputs(i) = newOutput(var,i,1:nb);

% solver settings
codeoptions = getOptions('beliefMPC');
codeoptions.printlevel = 0;
codeoptions.timing=0;

% generate code
generateCode(stages,params,codeoptions,outputs);

end