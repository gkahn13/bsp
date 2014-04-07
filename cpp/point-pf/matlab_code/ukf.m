function ukf()

rng('default');
format shortg;

xdim = 1; wdim = 1;
zdim = 1; vdim = 1;

T = 100;

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

mean_ukf = zeros(xdim,T+1);
cov_ukf = cell(1,T+1);

mean_ukf(:,1) = x0;
cov_ukf{1} = P0;
x = x0;
P = P0;

for t=1:T
    [xSigmaPts, wSigmaPts, nsp] = UT(x, P);
        
    wSigmaPts_xmat = repmat(wSigmaPts(:,2:nsp),xdim,1);
    
    xPredSigmaPts = zeros(xdim,nsp);
    for i=1:nsp
        xPredSigmaPts(:,i) = f(xSigmaPts(:,i),t,zeros(wdim,1));
    end
    xPred = sum(wSigmaPts_xmat .* (xPredSigmaPts(:,2:nsp) - repmat(xPredSigmaPts(:,1),1,nsp-1)),2);
    xPred = xPred+xPredSigmaPts(:,1);
    
    exSigmaPt = xPredSigmaPts(:,1)-xPred;
    PPred   = wSigmaPts(nsp+1)*(exSigmaPt*exSigmaPt');
    
    exSigmaPt = xPredSigmaPts(:,2:nsp) - repmat(xPred,1,nsp-1);
    PPred   = PPred + (wSigmaPts_xmat .* exSigmaPt) * exSigmaPt' + Q;
    
    [xPredSigmaPts, wSigmaPts, nsp] = UT(xPred, PPred);
    wSigmaPts_zmat = repmat(wSigmaPts(:,2:nsp),zdim,1);
    
    zPredSigmaPts = zeros(zdim,nsp);
    for i=1:nsp
        zPredSigmaPts(:,i) = h(xPredSigmaPts(:,i),zeros(vdim,1));
    end
    zPred = sum(wSigmaPts_zmat .* (zPredSigmaPts(:,2:nsp) - repmat(zPredSigmaPts(:,1),1,nsp-1)),2);
    zPred = zPred+zPredSigmaPts(:,1);
    
    exSigmaPt = xPredSigmaPts(:,1)-xPred;
    ezSigmaPt = zPredSigmaPts(:,1)-zPred;
    PxzPred   = wSigmaPts(nsp+1)*exSigmaPt*ezSigmaPt';
    S         = wSigmaPts(nsp+1)*(ezSigmaPt*ezSigmaPt');
    
    exSigmaPt = xPredSigmaPts(:,2:nsp) - repmat(xPred,1,nsp-1);
    ezSigmaPt = zPredSigmaPts(:,2:nsp) - repmat(zPred,1,nsp-1);
    S         = S + (wSigmaPts_zmat .* ezSigmaPt) * ezSigmaPt' + R;
    PxzPred   = PxzPred + exSigmaPt * (wSigmaPts_zmat .* ezSigmaPt)';

    % measurement update
    K  = PxzPred / S;
    innovation = Z(:,t+1) - zPred;
    x = xPred + K*innovation;
    P = PPred - K*S*K';
    
    mean_ukf(:,t+1) = x;
    cov_ukf{t+1} = P;
end

figure('units','pixel','outerposition',  [0 0 1200 800]);
hor = 1:T+1;
clf; hold on
ff2 = [mean_ukf' + 2*sqrt(cell2mat(cov_ukf')); flipdim(mean_ukf' - 2*sqrt(cell2mat(cov_ukf')),1)];
fill([hor';  flipdim(hor',1)], ff2, [7 7 8]/8, 'EdgeColor', [7 7 8]/8);
plot(hor, X, 'r--','linewidth',5);
%plot(hor, mean_ukf, 'b--','linewidth',5);
plot(hor, mean_ukf + 2*sqrt(cell2mat(cov_ukf)), 'color', [0 4 0]/8);
plot(hor, mean_ukf - 2*sqrt(cell2mat(cov_ukf)), 'color', [0 4 0]/8);
  
grid on;
set(gcf,'PaperSize', [10 5]);
set(gcf,'PaperPosition',[0.1 0.1 10 5]);
xlabel('time steps');
ylabel('hidden states');
axis([0 T+1 -30 30]);
legend('UKF','ground truth','location','southwest');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnew = f(x, t, w)
xnew = x*0.5 + (25/(1 + x'*x))*x + 8*cos(1.2*t)*ones(size(x,1),1) + w;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = h(x, n)
z = 5*sin(x*3) + n;
%z = (x'*x)/20 + n;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xPts, wPts, nPts] = UT(x,P)
alpha=1;
beta=0;
kappa=2;
n = size(x(:),1);
nPts = 2*n+1;
kappa = alpha^2*(n+kappa)-n;
Psqrtm=(chol((n+kappa)*P))';  
xPts=[zeros(size(P,1),1) -Psqrtm Psqrtm];
xPts = xPts + repmat(x,1,nPts);  
wPts=[kappa 0.5*ones(1,nPts-1) 0]/(n+kappa);
wPts(nPts+1) = wPts(1) + (1-alpha^2) + beta;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%