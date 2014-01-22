function arm_state_mpc_gen(timesteps)

% currently goal is a vector of joints
% need to change if we want the goal to be a position

% FORCES - Fast interior point code generation for multistage problems.
% Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
% Automatic Control Laboratory, ETH Zurich.

currDir = pwd;
disp('currDir');
disp(currDir);
endPwdIndex = strfind(currDir,'bsp') - 1;
rootDir = currDir(1:endPwdIndex);
forcesDir = strcat(rootDir,'bsp/forces');
addpath(forcesDir);
disp(strcat(rootDir,'bsp/forces'));

% problem setup
N = timesteps - 1;
nx = 6;
nu = 6;
stages = MultistageProblem(N+1);

% first stage
i=1;
istr = sprintf('%d',i);

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
    istr = sprintf('%d',i);
    
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
istr = sprintf('%d',i);

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
mpcname = 'armStateMPC';
codeoptions = getOptions(mpcname);
codeoptions.printlevel = 0;
codeoptions.timing=0;

% generate code
generateCode(stages,params,codeoptions,outputs);


disp('Unzipping into mpc...');
outdir = 'mpc/';
system(['mkdir -p ',outdir]);
header_file = [mpcname,num2str(timesteps),'.h'];
src_file = [mpcname,num2str(timesteps),'.c'];
system(['unzip -p ',mpcname,'.zip include/',mpcname,'.h > ',outdir,header_file]);
system(['unzip -p ',mpcname,'.zip src/',mpcname,'.c > ',outdir,src_file]);
system(['rm -rf ',mpcname,'.zip @CodeGen']);

disp('Replacing incorrect #include in .c file...');
str_to_delete = ['#include "../include/',mpcname,'.h"'];
str_to_insert = ['#include "',mpcname,'.h"'];
mpc_src = fileread([outdir,src_file]);
mpc_src_new = strrep(mpc_src,str_to_delete,str_to_insert);

fid = fopen([outdir,src_file],'w');
fwrite(fid,mpc_src_new);
fclose(fid);

end