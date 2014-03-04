function parameter_state_penalty_mpc_gen(timesteps)

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
alpha_belief = 10;
alpha_final_belief = 10;
alpha_control = 1;
N = timesteps - 1;
nx = 8;

nu = 2;
nj = 4; % number of joints
np = 4; % number of parameters to be estimated
stages = MultistageProblem(N+1);


% first stage
i=1;
istr = sprintf('%d',i);

% dimensions
stages(i).dims.n = 3*nx+nu;           % number of stage variables
stages(i).dims.l = 3*nx+nu;           % number of lower bounds
stages(i).dims.u = nx+nu;             % number of upper bounds
stages(i).dims.r = nx;              % number of equality constraints
stages(i).dims.p = 0;                 % number of affine constraints
stages(i).dims.q = 0;                 % number of quadratic constraints

% cost

params(1) = newParam(['H',istr], i, 'cost.H','diag');
params(end+1) = newParam(['f',istr], i, 'cost.f');

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.l; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.u; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable

params(end+1) = newParam(['C',istr], i, 'eq.C');
params(end+1) = newParam(['e',istr], i, 'eq.c');
stages(i).eq.D = [eye(nx), zeros(nx,2*nx+nu)];

for i = 2:N
    istr = sprintf('%d',i);
    
    % dimension
    stages(i).dims.n = 3*nx+nu;    % number of stage variables
    stages(i).dims.l = 3*nx+nu;    % number of lower bounds
    stages(i).dims.u = nx+nu;      % number of upper bounds
    stages(i).dims.r = nx;         % number of equality constraints
    stages(i).dims.p = 0;          % number of polytopic constraints
    stages(i).dims.q = 0;          % number of quadratic constraints
    
    % cost
    %stages(i).cost.H = 2*blkdiag(zeros(nx,nx), Q, R, zeros(nb,nb), zeros(nb, nb));
    params(end+1) = newParam(['H',istr], i, 'cost.H','diag');
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
    stages(i).eq.D = [-eye(nx), zeros(nx,2*nx+nu)];
    
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
params(end+1) = newParam(['H',istr], i, 'cost.H', 'diag');
params(end+1) = newParam(['f',istr], i, 'cost.f');
%stages(i).cost.f = zeros(stages(i).dims.n,1);

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.l; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb');

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.u; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub');

% equality constraints
%stages(i).eq.C = blkdiag(eye(nx), zeros(ns,ns));
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
mpcname = 'statePenaltyMPC';
codeoptions = getOptions(mpcname);
codeoptions.printlevel = 0;
codeoptions.timing=0;
codeoptions.maxit=60;

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

function [H] = hess_structure(alpha_belief, nj, np)

Q = blkdiag(zeros(nj,nj), alpha_belief*eye(np,np));

QTQ = Q'*Q;
n = size(QTQ,1);
sdim = (n*(n+1))/2;

% construct the full hessian
Qfull = [];
for i=1:n
Qfull = blkdiag(Qfull,QTQ);
end

%fprintf('Full trace: %f\n',vec(SS)'*Qfull*vec(SS));

%s = SS(find(tril(ones(n,n))));

A = zeros(n*n,sdim);

idx = 0;
% Written w.r.t c++ indices, added 1 for matlab 1-based indexing
for i = 0:n-1
    for j=0:n-1
        if (i <= j)
            A(idx+1,(2*n-i+1)*i/2+j-i+1) = 1;
        else
            A(idx+1,(2*n-j+1)*j/2+i-j+1) = 1;
        end
        idx = idx+1;
    end
end

H = A'*Qfull*A;
%disp(H)
%fprintf('Reduced trace: %f\n',s'*H*s);

end