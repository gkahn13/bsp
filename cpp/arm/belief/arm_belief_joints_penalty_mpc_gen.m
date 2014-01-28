function arm_belief_joints_penalty_mpc_gen(timesteps)

% FORCES - Fast interior point code generation for multistage problems.
% Copyright (C) 2011-12 Alexander Domahidi [domahidi@control.ee.ethz.ch],
% Automatic Control Laboratory, ETH Zurich.

%currDir = pwd;
%disp('currDir');
%disp(currDir);
%endPwdIndex = strfind(currDir,'bsp') - 1;
%rootDir = currDir(1:endPwdIndex);
%forcesDir = strcat(rootDir,'bsp/forces');
%addpath(forcesDir);
%disp(strcat(rootDir,'bsp/forces'));

% problem setup
N = timesteps - 1;
nx = 6;
ns = ((nx+1)*nx)/2;
nb = nx+ns;
nu = nx;
stages = MultistageProblem(N+1);

% first stage
i=1;
istr = sprintf('%d',i);

% dimensions
stages(i).dims.n = 3*nb+nu;           % number of stage variables
stages(i).dims.l = 3*nb+nu;           % number of lower bounds
stages(i).dims.u = nb+nu;             % number of upper bounds
stages(i).dims.r = nb;                % number of equality constraints
stages(i).dims.p = 0;                 % number of affine constraints
stages(i).dims.q = 0;                 % number of quadratic constraints

% cost
params(1) = newParam(['Q',istr], i, 'cost.H');
params(end+1) = newParam(['f',istr], i, 'cost.f');

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.l; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.u; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable

params(end+1) = newParam(['C',istr], i, 'eq.C');
params(end+1) = newParam(['e',istr], i, 'eq.c');
stages(i).eq.D = [eye(nb), zeros(nb,2*nb+nu)];

for i = 2:N
    istr = sprintf('%d',i);
    
    % dimension
    stages(i).dims.n = 3*nb+nu;    % number of stage variables
    stages(i).dims.l = 3*nb+nu;    % number of lower bounds
    stages(i).dims.u = nb+nu;      % number of upper bounds
    stages(i).dims.p = 0;          % number of polytopic constraints
    stages(i).dims.q = 0;          % number of quadratic constraints
    stages(i).dims.r = nb;         % number of equality constraints
    
    % cost
    params(end+1) = newParam(['Q',istr], i, 'cost.H');
    params(end+1) = newParam(['f',istr], i, 'cost.f');
    
    % lower bounds
    stages(i).ineq.b.lbidx = 1:stages(i).dims.l; % lower bound acts on these indices
    params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable
    
    % upper bounds
    stages(i).ineq.b.ubidx = 1:stages(i).dims.u; % upper bound acts on these indices
    params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable
        
    % equality constraints
    params(end+1) = newParam(['C',istr], i, 'eq.C');
    params(end+1) = newParam(['e',istr], i, 'eq.c');
    stages(i).eq.D = [-eye(nb), zeros(nb,2*nb+nu)];
    
end

% final stage
i = N+1;
istr = sprintf('%d',i);

% dimension
stages(i).dims.n = nb;    % number of stage variables
stages(i).dims.r = nb;    % number of equality constraints
stages(i).dims.l = nb;    % number of lower bounds
stages(i).dims.u = nb;    % number of upper bounds
stages(i).dims.p = 2*nx;     % number of polytopic constraints
stages(i).dims.q = 0;     % number of quadratic constraints

% cost
params(end+1) = newParam(['Q',istr], i, 'cost.H');
params(end+1) = newParam(['f',istr], i, 'cost.f');

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.l; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb');

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.u; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub');

% equality constraints
%stages(i).eq.C = blkdiag(eye(nx), zeros(ns,ns));
params(end+1) = newParam(['e',istr], i, 'eq.c');
stages(i).eq.D = -eye(nb);

stages(i).ineq.p.A = [eye(nx), zeros(nx,ns);-eye(nx), zeros(nx,ns)];
params(end+1) = newParam(['b',istr], i, 'ineq.p.b');

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
mpcname = 'beliefPenaltyMPC';
codeoptions = getOptions(mpcname);
codeoptions.printlevel = 0;
codeoptions.timing = 0;
codeoptions.maxit = 40;

% generate code
generateCode(stages,params,codeoptions,outputs);

%disp('Unzipping into mpc...');
%outdir = 'mpc/';
%system(['mkdir -p ',outdir]);
%header_file = [mpcname,num2str(timesteps),'.h'];
%src_file = [mpcname,num2str(timesteps),'.c'];
%system(['unzip -p ',mpcname,'.zip include/',mpcname,'.h > ',outdir,header_file]);
%system(['unzip -p ',mpcname,'.zip src/',mpcname,'.c > ',outdir,src_file]);
%system(['rm -rf ',mpcname,'.zip @CodeGen']);

%disp('Replacing incorrect #include in .c file...');
%str_to_delete = ['#include "../include/',mpcname,'.h"'];
%str_to_insert = ['#include "',mpcname,'.h"'];
%mpc_src = fileread([outdir,src_file]);
%mpc_src_new = strrep(mpc_src,str_to_delete,str_to_insert);

%fid = fopen([outdir,src_file],'w');
%fwrite(fid,mpc_src_new);
%fclose(fid);

end