function point_control_mpc_gen(timesteps)

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
nu = 2;
stages = MultistageProblem(N);

% first stage
i=1;
istr = sprintf('%d',i);

% dimensions
stages(i).dims.n = nu;           % number of stage variables
stages(i).dims.l = nu;           % number of lower bounds
stages(i).dims.u = nu;           % number of upper bounds
stages(i).dims.r = nu;            % number of equality constraints
stages(i).dims.p = 0;            % number of affine constraints
stages(i).dims.q = 0;            % number of quadratic constraints

% cost
params(1) = newParam(['H',istr], i, 'cost.H', 'diag'); % diagonal hessian
params(end+1) = newParam(['f',istr], i, 'cost.f'); % gradient

% lower bounds
stages(i).ineq.b.lbidx = 1:stages(i).dims.n; % lower bound acts on these indices
params(end+1) = newParam(['lb',istr], i, 'ineq.b.lb'); % lower bound for this stage variable

% upper bounds
stages(i).ineq.b.ubidx = 1:stages(i).dims.n; % upper bound acts on these indices
params(end+1) = newParam(['ub',istr], i, 'ineq.b.ub'); % upper bound for this stage variable

stages(i).eq.C = zeros(nu,nu);
stages(i).eq.c = zeros(nu,1);

for i = 2:N
    istr = sprintf('%d',i);
    
    % dimension
    stages(i).dims.n = nu;    % number of stage variables
    stages(i).dims.l = nu;    % number of lower bounds
    stages(i).dims.u = nu;    % number of upper bounds
    stages(i).dims.r = nu;       % number of equality constraints
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
    
    stages(i).eq.C = zeros(nu,nu);
    stages(i).eq.c = zeros(nu,1);
    stages(i).eq.D = zeros(nu,nu);
end

%--------------------------------------------------------------------------
% define outputs of the solver
for i=1:N
    var = sprintf('z%d',i);
    outputs(i) = newOutput(var,i,1:nu);
end

% solver settings
mpcname = 'controlMPC';
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