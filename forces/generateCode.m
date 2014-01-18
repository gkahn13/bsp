function success = generateCode(stages,params,settings,outvars)
% Generate a custom interior point solver for multistage problems using
% FORCES.
%
%    SUCCESS = GENERATECODE(STAGES) generates, downloads and compiles your
%    custom solver for the multistage problem STAGES. Default settings are
%    used, and the default parameter is 'stages(1).eq.c', i.e. the offset
%    term in the first equality constraint (as typically used in MPC). The
%    default output are simply all stage variables.
%
%    SUCCESS = GENERATECODE(STAGES,PARAMS) does the above but with user
%    defined parameters.
%
%    SUCCESS = GENERATECODE(STAGES,PARAMS,SETTINGS) does the above but with
%    user defined parameters and settings.
%
%    SUCCESS = GENERATECODE(STAGES,PARAMS,SETTINGS,OUTVARS) does the above 
%    but with user defined parameters, settings and outputs.
%
% SEE ALSO MULTISTAGEPROBLEM NEWPARAM NEWOUTPUT GETOPTIONS FORCES_LICENSE

if( nargin < 4 )
    % default output variables
    outvars = getAllOutputs(stages);
end
if( nargin < 3 )
    % default codegen options
    settings = getOptions;
end
if( nargin < 2 )
    % default parameter setting
    params = newParam('z0', 1, 'eq.c');
end

% copyright
fprintf('\n');
fprintf('FORCES - Fast interior point code generation for multistage problems.\n');
fprintf('Copyright (C) 2011-13 Alexander Domahidi [domahidi@control.ee.ethz.ch]\n');
fprintf('Automatic Control Laboratory, ETH Zurich.\n\n');

      
forcesurl = 'http://forces.ethz.ch';
createClassFromWsdl([forcesurl,'/CodeGen.asmx?Wsdl']);
obj = CodeGen;

% save structs in matfile
disp('Preparing data to be sent...');
[~,tmpname] = fileparts(tempname);
matfile = [tmpname,'.mat'];
save([tmpname,'.mat'],'stages','params','settings','outvars','-v6');

% read mat file as string stream
fid = fopen(matfile,'r');
byteStream = fread(fid);
fclose(fid);

% generate code
disp('Sending request to server...');
tic;
try
    zipurl = [forcesurl,generateCodeFromMatlab(obj, '135bd27a-27a6-4a8a-81c8-77a605c798d3', byteStream)];
catch err
    delete(matfile);
    throw(err);
end
t = toc;
fprintf('Code successfully generated in %4.2f sec. ',t);

% delete temporary file
delete(matfile);

% download code
disp('Downloading generated code...');
[~,file,ext] = fileparts(zipurl);
outdir = file;
filename = [file,ext];
urlwrite(zipurl,filename);
disp(['Code has been downloaded to ''',filename,'''. Extracting...']);

% check overwrite
if(isfield(settings, 'overwrite'))
    if(settings.overwrite == 0)         % never overwrite
        overwrite = 0;
    elseif(settings.overwrite == 1)     % always overwrite
        overwrite = 1; 
    else                                % ask to overwrite
        overwrite = 2;
    end
else                                    % ask to overwrite  
   overwrite = 2;    
end

% unzip
% dcontent = dir(outdir);
% if( size(dcontent,1) > 0 )
%     if(overwrite == 1)
%         rmdir(outdir,'s');
%         unzip(filename,file);
%         disp('Code successfully unzipped. MEXing your solver...');
%     elseif(overwrite == 2)
%         confirm = input(['Directory ''',outdir,''' will be overwritten. Proceed? [y]/n '],'s');
%         if( strcmp(confirm,'y') || isempty(confirm) )
%             rmdir(outdir,'s');
%             unzip(filename,file);
%             disp('Code successfully unzipped. MEXing your solver...');
%         end
%     end
% else
%     unzip(filename,file);
%     disp('Code successfully unzipped. MEXing your solver...');
% end

% make mex file
% clear(file);
% cd([file,'/interface']);
% makemex;
% cd ..
% cd ..
% disp(' ');
% disp('Code generation and MEXing successfully completed.');
% disp(['Type ''help ',file,''' for more documentation. Happy solving!']);

% return success
success = 1;