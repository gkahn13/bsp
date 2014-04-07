function pf()

%rng('default');
format shortg;

xdim = 1; wdim = 1;
zdim = 1; vdim = 1;

T = 200;
N = 1000;

x0 = zeros(xdim,1);
P0 = 5*eye(xdim,1);

Q = 1*eye(wdim);
R = 10*eye(vdim);

X = zeros(xdim,T);
Z = zeros(zdim,T);

X(:,1) = x0 + chol(P0)'*randn(xdim,1);
Z(:,1) = h(X(:,1), chol(R)'*randn(vdim,1));

for t=1:T
  X(:,t+1) = f(X(:,t),t,chol(Q)'*randn(wdim,1));
  Z(:,t+1) = h(X(:,t+1),chol(R)'*randn(vdim,1));
end

mean_pf = zeros(xdim,T+1);
mean_pf(:,1) = x0;

Xp = repmat(x0,1,N) + chol(P0)'*randn(xdim,N);
pw = zeros(1,N);

figure('units','pixel','outerposition',  [0 0 1200 800]);
hor = 1:T+1;
clf; hold on

for t=1:T
    fprintf('Timestep: %i\n',t);
    
    %if (t==1 || rem(t,5)==0)
    %    plot(t+1, Xp(1,:),'ko');
    %end
    
    for i = 1:N
        Xp(i) = f(Xp(i),t,chol(Q)'*randn(wdim,1));
        e = Z(:,t+1) - h(Xp(i),zeros(vdim,1));
        pw(i) = gauss_likelihood(e,R);
    end
    pw = pw/sum(pw);
    mean_pf(:,t+1) = pw*Xp';
    
    Xp = resample(Xp,pw);
end
    
plot(hor, mean_pf, 'b--','linewidth',5);
plot(hor, X, 'r--','linewidth',5);

grid on;
set(gcf,'PaperSize', [10 5]);
set(gcf,'PaperPosition',[0.1 0.1 10 5]);
xlabel('time steps');
ylabel('hidden states');
axis([0 T+1 -30 30]);
%legend('PF','ground truth','location','southwest');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnew = f(x, t, w)
xnew = x*0.5 + (25/(1 + x'*x))*x + 8*cos(1.2*t)*ones(size(x,1),1) + w;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x, n)
z = 5*sin(x) + n;
%z = (x'*x)/20 + n;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = resample(Xp,pw)
cdf = cumsum(pw);
diff = cdf'*ones(1,length(pw)) - ones(length(pw),1)*rand(1,length(pw));
diff = (diff <= 0) * 2 + diff;
[~, idx] = min(diff);
P = Xp(idx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = gauss_likelihood(v,S)
D = size(v,1);
Sf = chol(S)';
M = Sf\v; % equivalent to writing inv(Sf)*v
% M is the normalised innovation of v, and M(:,i)'*M(:,i) gives the Mahalanobis 
% distance for each v(:,i). 
E = -0.5 * sum(M.*M, 1); % term inside exponential of Gaussian formula
% Note: writing sum(x.*x, 1) is a fast way to compute sets of inner-products.
C = (2*pi)^(D/2) * prod(diag(Sf)); % normalising term (makes Gaussian hyper-volume equal 1)
w = exp(E) / C; % likelihood
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
