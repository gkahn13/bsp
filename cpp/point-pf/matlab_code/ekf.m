function ekf()

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

mean_ekf = zeros(xdim,T+1);
cov_ekf = cell(1,T+1);

mean_ekf(:,1) = x0;
cov_ekf{1} = P0;
x = x0;
P = P0;

for t=1:T
    % predict
    F = numerical_jacobian(@f,1,x,t,zeros(wdim,1));
    G = numerical_jacobian(@f,3,x,t,zeros(wdim,1));
    x = f(x,t,zeros(xdim,1));
    P = F*P*F' + G*Q*G';
    
    % update
    zpred = h(x,zeros(vdim,1));
    H = numerical_jacobian(@h,1,x,zeros(vdim,1));
    V = numerical_jacobian(@h,2,x,zeros(vdim,1));
    P12 = P*H';
    Pz = H*P12 + V*R*V';
    R1 = chol(Pz);
    U = P12/R1;
    x = x + U*(R1'\(Z(:,t+1)-zpred));
    P = P - U*U';
    
    mean_ekf(:,t+1) = x;
    cov_ekf{t+1} = P;
    
    %fprintf('Mean: %f cov: %f\n',mean_ekf(:,t+1),cov_ekf{t+1});
 end

figure('units','pixel','outerposition',  [0 0 1200 800]);
hor = 1:T+1;
clf; hold on
ff2 = [mean_ekf' + 2*sqrt(cell2mat(cov_ekf')); flipdim(mean_ekf' - 2*sqrt(cell2mat(cov_ekf')),1)];
fill([hor';  flipdim(hor',1)], ff2, [7 7 8]/8, 'EdgeColor', [7 7 8]/8);
plot(hor, X, 'r--','linewidth',5);
%plot(hor, mean_ekf, 'b--','linewidth',5);
plot(hor, mean_ekf + 2*sqrt(cell2mat(cov_ekf)), 'color', [0 4 0]/8);
plot(hor, mean_ekf - 2*sqrt(cell2mat(cov_ekf)), 'color', [0 4 0]/8);
  
grid on;
set(gcf,'PaperSize', [10 5]);
set(gcf,'PaperPosition',[0.1 0.1 10 5]);
xlabel('time steps');
ylabel('hidden states');
axis([0 T+1 -30 30]);
legend('EKF','ground truth','location','southwest');

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
function J = numerical_jacobian(model, idx, varargin)

step = 1e-6; %0.00390625

x = varargin{idx};
y = feval(model, varargin{:});
lenx = length(x);
leny = length(y);
J = zeros(leny, lenx);

for i=1:lenx
    xu = x(i) + step;
    xl = x(i) - step;

    varargin{idx}(i) = xu;
    yu = feval(model, varargin{:}); 
    varargin{idx}(i) = xl;
    yl = feval(model, varargin{:});    
    varargin{idx}(i) = x(i);
    
    J(:,i) = (yu - yl)/(xu - xl);  
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%