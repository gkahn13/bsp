function hess_structure()
clear all;

load SS.mat;
load J.mat;

JTJ = J'*J;
n = size(JTJ,1);
sdim = (n*(n+1))/2;

% construct the full hessian
Hfull = [];
for i=1:n
Hfull = blkdiag(Hfull,JTJ);
end

fprintf('Full trace: %f\n',vec(SS)'*Hfull*vec(SS));

s = SS(find(tril(ones(n,n))));

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

H = A'*Hfull*A;
disp(H)
fprintf('Reduced trace: %f\n',s'*H*s);

end