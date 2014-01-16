function statebspgen()

% FORCES - Fast interior point code generation for multistage problems.
% Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
% Automatic Control Laboratory, ETH Zurich.

close all;
clear all;

% problem setup
N = 14;
nx = 2;
nu = 2;
stages = MultistageProblem(N+1);

% first stage
i=1;
istr = sprintf('%02d',i);

% dimensions
stages(i).dims.n = nx+nu;           % number of stage variables
stages(i).dims.l = nx+nu;           % number of lower bounds
stages(i).dims.u = nx+nu;           % number of upper bounds
stages(i).dims.r = nx;              % number of equality constraints
stages(i).dims.p = 0;               % number of affine constraints
stages(i).dims.q = 0;               % number of quadratic constraints

% cost
params(1) = newParam(['H',istr], i, 'cost.H', 'diag'); % diagonal hessian
params(end+1) = newParam(['f',istr], i, 'cost.f'); % gradient

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable

params(end+1) = newParam(['C',istr], i, 'eq.C');
params(end+1) = newParam(['e',istr], i, 'eq.c');
stages(i).eq.D = [eye(nx), zeros(nx,nu)];

for i = 2:N
    istr = sprintf('%02d',i);
    
    % dimension
    stages(i).dims.n = nx+nu;    % number of stage variables
    stages(i).dims.l = nx+nu;    % number of lower bounds
    stages(i).dims.u = nx+nu;    % number of upper bounds
    stages(i).dims.r = nx;       % number of equality constraints
    stages(i).dims.p = 0;        % number of polytopic constraints
    stages(i).dims.q = 0;        % number of quadratic constraints
    
    % cost
    params(end+1) = newParam(['H',istr], i, 'cost.H', 'diag'); % diagonal Hessian
    params(end+1) = newParam(['f',istr], i, 'cost.f'); % gradient
    
    % lower bounds
    stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
    params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable
    
    % upper bounds
    stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
    params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable
        
    % equality constraints
    params(end+1) = newParam(['C',istr], i, 'eq.C');
    params(end+1) = newParam(['e',istr], i, 'eq.c');
    stages(i).eq.D = [-eye(nx), zeros(nx,nu)]; 
end

% final stage
i = N+1;
istr = sprintf('%02d',i);

% dimension
stages(i).dims.n = nx;    % number of stage variables
stages(i).dims.r = nx;    % number of equality constraints
stages(i).dims.l = nx;    % number of lower bounds
stages(i).dims.u = nx;    % number of upper bounds
stages(i).dims.p = 0;     % number of polytopic constraints
stages(i).dims.q = 0;     % number of quadratic constraints

% cost
params(end+1) = newParam(['H',istr], i, 'cost.H', 'diag'); % diagonal Hessian
params(end+1) = newParam(['f',istr], i, 'cost.f'); % gradient

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb');

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub');

% equality constraints
params(end+1) = newParam(['e',istr], i, 'eq.c');
stages(i).eq.D = -eye(nx);

%--------------------------------------------------------------------------
% define outputs of the solver
for i=1:N
    var = sprintf('z%d',i);
    outputs(i) = newOutput(var,i,1:nx+nu);
end
i=N+1;
var = sprintf('z%d',i);
outputs(i) = newOutput(var,i,1:nx);

% solver settings
codeoptions = getOptions('stateMPC');
codeoptions.printlevel = 0;
codeoptions.timing = 0;

% generate code
generateCode(stages,params,codeoptions,outputs);

end