function opts = getOptions(varargin)
% Returns options structure for generating code with FORCES.
% 
%    OPTS = GETOPTIONS returns a default options struct.
%
%    OPTS = GETOPTIONS(NAME) returns a default options struct with the 
%    solver named NAME.
%
% See also GENERATECODE FORCES_LICENSE

% default name
if (nargin == 1 && isa(varargin{1},'char') )
    opts.name = varargin{1};
else
    opts.name = 'FORCES';
end

% maximum number of iterations
opts.maxit = 30;

% default line search options
opts.linesearch.factor_aff = 0.9;
opts.linesearch.factor_cc = 0.95;
opts.linesearch.minstep = 1e-8;
opts.linesearch.maxstep = 0.995;

% default accuracy / convergence options
opts.accuracy.mu = 1e-6;
opts.accuracy.ineq = 1e-6;
opts.accuracy.eq = 1e-6;
opts.accuracy.rdgap = 1e-4;

% printlevel
opts.printlevel = 2;

% initialization
opts.init = 0;

% optimizationlevel
opts.optlevel = 1;

% overwrite
opts.overwrite = 2;		% 0 never overwrite existing code
						% 1 always overwrites existing code
						% 2 ask to overwrite existing code

% timing
opts.timing = 1;

% datatype
opts.floattype = 'double';

% parallel
opts.parallel = 0;		% 0 design code for single core
						% 1 design code for multiple cores
                        
% regularization: if Mii < epsilon then Mii = delta, where
%                 Mii is the ith diagonal during Cholesky fact.
opts.regularize.epsilon = 1e-13;
opts.regularize.delta = 4e-4;

